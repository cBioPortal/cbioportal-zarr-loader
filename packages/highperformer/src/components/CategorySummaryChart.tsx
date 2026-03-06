import { useState } from 'react'
import { Group } from '@visx/group'
import { Pie } from '@visx/shape'
import { scaleBand, scaleLinear } from '@visx/scale'
import { AxisBottom, AxisLeft } from '@visx/axis'
import { Popover, Segmented, Typography } from 'antd'
import { ExpandOutlined } from '@ant-design/icons'
import type { SelectionGroup } from '../store/useAppStore'
import type { RGB } from '../utils/colors'
import { ALL_CELLS_GROUP_ID } from '../hooks/useAllCellsSummary'
import ChartModal from './ChartModal'

function groupLabel(id: number): string {
  return id === ALL_CELLS_GROUP_ID ? 'All Cells' : `Group ${id}`
}

interface CategorySummaryChartProps {
  name: string
  categoryMap: { label: string; color: RGB }[]
  countsByGroup: Map<number, Uint32Array>
  groups: SelectionGroup[]
}

interface CategoryRow {
  label: string
  counts: Record<number, number>
  total: Record<number, number>
}

type StackMode = 'count' | 'percent'

const MARGIN = { top: 4, right: 12, bottom: 4, left: 100 }
const BAR_HEIGHT = 16
const MAX_CATS_INLINE = 5

function buildData(
  categoryMap: { label: string; color: RGB }[],
  countsByGroup: Map<number, Uint32Array>,
  activeGroups: SelectionGroup[],
): CategoryRow[] {
  const data: CategoryRow[] = []
  for (let ci = 0; ci < categoryMap.length; ci++) {
    const counts: Record<number, number> = {}
    let hasAny = false
    const totals: Record<number, number> = {}
    for (const g of activeGroups) {
      const c = countsByGroup.get(g.id)
      const count = c ? c[ci] : 0
      counts[g.id] = count
      const total = c ? Array.from(c).reduce((a, b) => a + b, 0) : 0
      totals[g.id] = total
      if (count > 0) hasAny = true
    }
    if (hasAny) {
      data.push({ label: categoryMap[ci].label, counts, total: totals })
    }
  }
  data.sort((a, b) => {
    const sumA = Object.values(a.counts).reduce((x, y) => x + y, 0)
    const sumB = Object.values(b.counts).reduce((x, y) => x + y, 0)
    return sumB - sumA
  })
  return data
}

interface GroupInfo {
  groupId: number
  color: string
  count: number
  pct: string
}

interface PopoverData {
  label: string
  highlightGroupId: number
  groups: GroupInfo[]
}

function StackedBarChart({ data, activeGroups, width, margin, mode }: {
  data: CategoryRow[]
  activeGroups: SelectionGroup[]
  width: number
  margin: typeof MARGIN
  mode: StackMode
}) {
  const [hoverKey, setHoverKey] = useState<string | null>(null)
  const [popover, setPopover] = useState<{ key: string; data: PopoverData } | null>(null)
  const rowHeight = BAR_HEIGHT + 4
  const height = data.length * rowHeight + margin.top + margin.bottom
  const innerWidth = width - margin.left - margin.right

  // Grand total across all categories and groups (for % of whole)
  let grandTotal = 0
  for (const d of data) {
    for (const g of activeGroups) grandTotal += d.counts[g.id] ?? 0
  }

  let maxValue = 1
  if (mode === 'count') {
    for (const d of data) {
      const sum = activeGroups.reduce((acc, g) => acc + (d.counts[g.id] ?? 0), 0)
      if (sum > maxValue) maxValue = sum
    }
  } else {
    // Find the largest category % to scale the axis
    for (const d of data) {
      const sum = activeGroups.reduce((acc, g) => acc + (d.counts[g.id] ?? 0), 0)
      const pct = grandTotal > 0 ? (sum / grandTotal) * 100 : 0
      if (pct > maxValue) maxValue = pct
    }
  }

  const yScale = scaleBand({
    domain: data.map((d) => d.label),
    range: [margin.top, height - margin.bottom],
    padding: 0.15,
  })

  const xScale = scaleLinear({
    domain: [0, maxValue],
    range: [0, innerWidth],
  })

  return (
    <svg
      width={width}
      height={height}
      onMouseLeave={() => setHoverKey(null)}
    >
      <Group left={margin.left}>
        {data.map((d) => {
          const y = yScale(d.label) ?? 0
          const barH = yScale.bandwidth()
          let xOffset = 0

          return activeGroups.map((g) => {
            const count = d.counts[g.id] ?? 0
            const value = mode === 'percent'
              ? (grandTotal > 0 ? (count / grandTotal) * 100 : 0)
              : count
            const barWidth = xScale(value)
            const x = xScale(xOffset)
            xOffset += value
            const key = `${d.label}-${g.id}`
            const rowKey = d.label
            const isHovered = hoverKey === key
            const isPopoverOpen = popover?.key === rowKey && popover.data.highlightGroupId === g.id
            const color = `rgb(${g.color.join(',')})`

            const handleClick = () => {
              if (popover?.key === rowKey && popover.data.highlightGroupId === g.id) {
                setPopover(null)
              } else {
                const groups: GroupInfo[] = activeGroups.map((ag) => {
                  const c = d.counts[ag.id] ?? 0
                  const p = grandTotal > 0 ? ((c / grandTotal) * 100).toFixed(1) : '0.0'
                  return { groupId: ag.id, color: `rgb(${ag.color.join(',')})`, count: c, pct: p }
                })
                setPopover({ key: rowKey, data: { label: d.label, highlightGroupId: g.id, groups } })
              }
            }

            const rectEl = (
              <rect
                key={key}
                fill={color}
                stroke={isHovered || isPopoverOpen ? '#333' : 'none'}
                strokeWidth={isHovered || isPopoverOpen ? 1.5 : 0}
                style={{
                  x, y, width: barWidth, height: barH,
                  transition: 'x 300ms ease, width 300ms ease',
                  cursor: 'pointer',
                }}
                onMouseEnter={() => setHoverKey(key)}
                onMouseLeave={() => setHoverKey(null)}
                onClick={handleClick}
              />
            )

            if (isPopoverOpen) {
              return (
                <Popover
                  key={key}
                  open
                  placement="left"
                  onOpenChange={(open) => { if (!open) setPopover(null) }}
                  content={
                    <div style={{ fontSize: 12 }}>
                      <div style={{ fontWeight: 600, marginBottom: 4 }}>{popover.data.label}</div>
                      {popover.data.groups.map((gi) => (
                        <div key={gi.groupId} style={{ fontWeight: gi.groupId === g.id ? 600 : 400 }}>
                          <span style={{ color: gi.color }}>{groupLabel(gi.groupId)}</span>: {gi.count.toLocaleString()} ({gi.pct}%)
                        </div>
                      ))}
                    </div>
                  }
                >
                  {rectEl}
                </Popover>
              )
            }

            return rectEl
          })
        })}
        <AxisLeft
          scale={yScale}
          tickValues={data.map((d) => d.label)}
          stroke="#e8e8e8"
          tickStroke="none"
          tickLabelProps={{
            fill: '#666',
            fontSize: 10,
            textAnchor: 'end',
            dy: '0.33em',
          }}
          tickFormat={(v) => {
            const s = String(v)
            return s.length > 18 ? s.slice(0, 16) + '\u2026' : s
          }}
        />
      </Group>
    </svg>
  )
}

const VERTICAL_MARGIN_BASE = { top: 12, right: 16, left: 48 }
const LABEL_MAX_CHARS = 30
const BAND_WIDTH = 32
const MIN_CHART_WIDTH = 400

function VerticalStackedBarChart({ data, activeGroups, mode }: {
  data: CategoryRow[]
  activeGroups: SelectionGroup[]
  mode: StackMode
}) {
  const [hoverKey, setHoverKey] = useState<string | null>(null)
  const [popover, setPopover] = useState<{ key: string; data: PopoverData } | null>(null)
  const maxLabelLen = Math.min(
    data.reduce((max, d) => Math.max(max, d.label.length), 0),
    LABEL_MAX_CHARS,
  )
  // fontSize 10, ~6px per char. Rotated 45° with textAnchor 'end':
  // vertical extent = len * 6 * sin(45°), horizontal extent = len * 6 * cos(45°)
  const COS45 = 0.707
  const charPx = 6
  const bottomMargin = Math.max(48, maxLabelLen * charPx * COS45 + 30)
  const leftMargin = Math.max(VERTICAL_MARGIN_BASE.left, maxLabelLen * charPx * COS45 + 8)
  const margin = { ...VERTICAL_MARGIN_BASE, bottom: bottomMargin, left: leftMargin }
  const innerWidth = Math.max(data.length * BAND_WIDTH, MIN_CHART_WIDTH)
  const width = innerWidth + margin.left + margin.right
  const innerHeight = 300

  let grandTotal = 0
  for (const d of data) {
    for (const g of activeGroups) grandTotal += d.counts[g.id] ?? 0
  }

  let maxValue = 1
  if (mode === 'count') {
    for (const d of data) {
      const sum = activeGroups.reduce((acc, g) => acc + (d.counts[g.id] ?? 0), 0)
      if (sum > maxValue) maxValue = sum
    }
  } else {
    for (const d of data) {
      const sum = activeGroups.reduce((acc, g) => acc + (d.counts[g.id] ?? 0), 0)
      const pct = grandTotal > 0 ? (sum / grandTotal) * 100 : 0
      if (pct > maxValue) maxValue = pct
    }
  }

  const xScale = scaleBand({
    domain: data.map((d) => d.label),
    range: [0, innerWidth],
    padding: 0.2,
  })

  const yScale = scaleLinear({
    domain: [0, maxValue],
    range: [innerHeight, 0],
    nice: true,
  })

  return (
    <svg
      width={width}
      height={innerHeight + margin.top + margin.bottom}
      onMouseLeave={() => setHoverKey(null)}
    >
      <Group left={margin.left} top={margin.top}>
        {data.map((d) => {
          const x = xScale(d.label) ?? 0
          const barW = Math.min(xScale.bandwidth(), BAND_WIDTH)
          const barOffset = (xScale.bandwidth() - barW) / 2
          let yOffset = 0

          return activeGroups.map((g) => {
            const count = d.counts[g.id] ?? 0
            const value = mode === 'percent'
              ? (grandTotal > 0 ? (count / grandTotal) * 100 : 0)
              : count
            const barHeight = innerHeight - yScale(value)
            const barY = yScale(yOffset + value)
            yOffset += value
            const key = `${d.label}-${g.id}`
            const rowKey = d.label
            const isHovered = hoverKey === key
            const isPopoverOpen = popover?.key === rowKey && popover.data.highlightGroupId === g.id
            const color = `rgb(${g.color.join(',')})`

            const handleClick = () => {
              if (popover?.key === rowKey && popover.data.highlightGroupId === g.id) {
                setPopover(null)
              } else {
                const groups: GroupInfo[] = activeGroups.map((ag) => {
                  const c = d.counts[ag.id] ?? 0
                  const p = grandTotal > 0 ? ((c / grandTotal) * 100).toFixed(1) : '0.0'
                  return { groupId: ag.id, color: `rgb(${ag.color.join(',')})`, count: c, pct: p }
                })
                setPopover({ key: rowKey, data: { label: d.label, highlightGroupId: g.id, groups } })
              }
            }

            const rectEl = (
              <rect
                key={key}
                fill={color}
                stroke={isHovered || isPopoverOpen ? '#333' : 'none'}
                strokeWidth={isHovered || isPopoverOpen ? 1.5 : 0}
                x={x + barOffset}
                y={barY}
                width={barW}
                height={barHeight}
                style={{ cursor: 'pointer', transition: 'y 300ms ease, height 300ms ease' }}
                onMouseEnter={() => setHoverKey(key)}
                onMouseLeave={() => setHoverKey(null)}
                onClick={handleClick}
              />
            )

            if (isPopoverOpen) {
              return (
                <Popover
                  key={key}
                  open
                  placement="top"
                  onOpenChange={(open) => { if (!open) setPopover(null) }}
                  content={
                    <div style={{ fontSize: 12 }}>
                      <div style={{ fontWeight: 600, marginBottom: 4 }}>{popover.data.label}</div>
                      {popover.data.groups.map((gi) => (
                        <div key={gi.groupId} style={{ fontWeight: gi.groupId === g.id ? 600 : 400 }}>
                          <span style={{ color: gi.color }}>{groupLabel(gi.groupId)}</span>: {gi.count.toLocaleString()} ({gi.pct}%)
                        </div>
                      ))}
                    </div>
                  }
                >
                  {rectEl}
                </Popover>
              )
            }

            return rectEl
          })
        })}
        <AxisLeft
          scale={yScale}
          numTicks={5}
          stroke="#e8e8e8"
          tickStroke="#e8e8e8"
          tickLabelProps={{ fill: '#666', fontSize: 10 }}
        />
        <AxisBottom
          scale={xScale}
          top={innerHeight}
          tickValues={data.map((d) => d.label)}
          stroke="#e8e8e8"
          tickStroke="none"
          tickLabelProps={{
            fill: '#666',
            fontSize: 10,
            textAnchor: 'end',
            angle: -45,
            dx: -4,
            dy: -2,
          }}
          tickFormat={(v) => {
            const s = String(v)
            return s.length > LABEL_MAX_CHARS ? s.slice(0, LABEL_MAX_CHARS - 2) + '\u2026' : s
          }}
        />
      </Group>
    </svg>
  )
}

const PIE_RADIUS = 50
const MAX_LEGEND_ITEMS = 8

interface PieDatum {
  label: string
  count: number
  color: string
}

function PieChart({ data, categoryMap, groupId, countsByGroup, radius = PIE_RADIUS, onShowMore }: {
  data: CategoryRow[]
  categoryMap: { label: string; color: RGB }[]
  groupId: number
  countsByGroup: Map<number, Uint32Array>
  radius?: number
  onShowMore?: () => void
}) {
  const [hoverLabel, setHoverLabel] = useState<string | null>(null)
  const [popoverLabel, setPopoverLabel] = useState<string | null>(null)

  const colorMap = new Map<string, string>()
  for (const cat of categoryMap) {
    colorMap.set(cat.label, `rgb(${cat.color.join(',')})`)
  }

  const counts = countsByGroup.get(groupId)
  const grandTotal = counts ? Array.from(counts).reduce((a, b) => a + b, 0) : 0

  const pieData: PieDatum[] = data.map((d) => ({
    label: d.label,
    count: d.counts[groupId] ?? 0,
    color: colorMap.get(d.label) ?? '#ccc',
  }))

  const size = radius * 2
  const legendItems = pieData.slice(0, MAX_LEGEND_ITEMS)
  const hasMore = pieData.length > MAX_LEGEND_ITEMS

  return (
    <div style={{ display: 'flex', gap: 12, alignItems: 'center' }}>
      <svg width={size} height={size} style={{ flexShrink: 0 }} onMouseLeave={() => setHoverLabel(null)}>
        <Group top={radius} left={radius}>
          <Pie
            data={pieData}
            pieValue={(d) => d.count}
            outerRadius={radius - 2}
            innerRadius={0}
            cornerRadius={1}
            padAngle={0.02}
          >
            {(pie) => pie.arcs.map((arc, i) => {
              const d = arc.data
              const isHovered = hoverLabel === d.label
              const isPopoverOpen = popoverLabel === d.label
              const pct = grandTotal > 0 ? ((d.count / grandTotal) * 100).toFixed(1) : '0.0'

              const pathEl = (
                <path
                  key={`arc-${i}`}
                  d={pie.path(arc) ?? ''}
                  fill={d.color}
                  stroke={isHovered || isPopoverOpen ? '#333' : 'white'}
                  strokeWidth={isHovered || isPopoverOpen ? 2 : 0.5}
                  style={{ cursor: 'pointer' }}
                  onMouseEnter={() => setHoverLabel(d.label)}
                  onMouseLeave={() => setHoverLabel(null)}
                  onClick={() => setPopoverLabel(isPopoverOpen ? null : d.label)}
                />
              )

              if (isPopoverOpen) {
                return (
                  <Popover
                    key={`arc-${i}`}
                    open
                    placement="right"
                    onOpenChange={(open) => { if (!open) setPopoverLabel(null) }}
                    content={
                      <div style={{ fontSize: 12 }}>
                        <div style={{ fontWeight: 600, marginBottom: 4 }}>{d.label}</div>
                        <div>{d.count.toLocaleString()} cells ({pct}%)</div>
                      </div>
                    }
                  >
                    {pathEl}
                  </Popover>
                )
              }

              return pathEl
            })}
          </Pie>
        </Group>
      </svg>
      <div style={{ fontSize: 10, lineHeight: '16px', minWidth: 0, overflow: 'hidden' }}>
        {legendItems.map((d) => {
          const pct = grandTotal > 0 ? ((d.count / grandTotal) * 100).toFixed(1) : '0.0'
          return (
            <div
              key={d.label}
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: 4,
                opacity: hoverLabel && hoverLabel !== d.label ? 0.4 : 1,
                fontWeight: hoverLabel === d.label ? 600 : 400,
                whiteSpace: 'nowrap',
                overflow: 'hidden',
                textOverflow: 'ellipsis',
              }}
              onMouseEnter={() => setHoverLabel(d.label)}
              onMouseLeave={() => setHoverLabel(null)}
            >
              <span style={{ width: 8, height: 8, borderRadius: 2, background: d.color, flexShrink: 0 }} />
              <span style={{ overflow: 'hidden', textOverflow: 'ellipsis' }}>{d.label}</span>
              <span style={{ color: '#999', flexShrink: 0 }}>{pct}%</span>
            </div>
          )
        })}
        {hasMore && (
          <div
            style={{ color: '#1677ff', cursor: 'pointer' }}
            onClick={onShowMore}
          >
            +{pieData.length - MAX_LEGEND_ITEMS} more
          </div>
        )}
      </div>
    </div>
  )
}

function CategoryTable({ data, activeGroups }: { data: CategoryRow[]; activeGroups: SelectionGroup[] }) {
  return (
    <table style={{ width: '100%', fontSize: 12, borderCollapse: 'collapse' }}>
      <thead>
        <tr>
          <th style={{ textAlign: 'left', fontWeight: 600, padding: '4px 8px', borderBottom: '1px solid #f0f0f0' }}>Value</th>
          {activeGroups.map((g) => (
            <th key={`count-${g.id}`} style={{ textAlign: 'right', fontWeight: 600, padding: '4px 8px', borderBottom: '1px solid #f0f0f0', color: `rgb(${g.color.join(',')})` }}>
              {groupLabel(g.id)} count
            </th>
          ))}
          {activeGroups.map((g) => (
            <th key={`pct-${g.id}`} style={{ textAlign: 'right', fontWeight: 600, padding: '4px 8px', borderBottom: '1px solid #f0f0f0', color: `rgb(${g.color.join(',')})` }}>
              {groupLabel(g.id)} %
            </th>
          ))}
        </tr>
      </thead>
      <tbody>
        {data.map((row) => (
          <tr key={row.label}>
            <td style={{ padding: '2px 8px' }}>{row.label}</td>
            {activeGroups.map((g) => (
              <td key={`count-${g.id}`} style={{ textAlign: 'right', padding: '2px 8px' }}>
                {(row.counts[g.id] ?? 0).toLocaleString()}
              </td>
            ))}
            {activeGroups.map((g) => {
              const count = row.counts[g.id] ?? 0
              const total = row.total[g.id] ?? 1
              const pct = total > 0 ? ((count / total) * 100).toFixed(1) : '0.0'
              return (
                <td key={`pct-${g.id}`} style={{ textAlign: 'right', padding: '2px 8px' }}>{pct}%</td>
              )
            })}
          </tr>
        ))}
      </tbody>
    </table>
  )
}

type InlineView = 'chart' | 'count' | 'percent' | 'table'

export default function CategorySummaryChart({ name, categoryMap, countsByGroup, groups }: CategorySummaryChartProps) {
  const isSingleGroup = groups.filter((g) => countsByGroup.has(g.id)).length === 1
  const [inlineView, setInlineView] = useState<InlineView>(isSingleGroup ? 'chart' : 'count')
  const [modalOpen, setModalOpen] = useState(false)
  const [modalTab, setModalTab] = useState<'chart' | 'table'>('chart')

  const activeGroups = groups.filter((g) => countsByGroup.has(g.id))

  if (activeGroups.length === 0) return null

  const data = buildData(categoryMap, countsByGroup, activeGroups)
  if (data.length === 0) return null

  // Filter out near-zero categories for the bar chart only
  const chartData = data.filter((d) =>
    activeGroups.some((g) => {
      const total = d.total[g.id]
      return total > 0 && (d.counts[g.id] / total) >= 0.005
    })
  )

  const inlineChartData = chartData.slice(0, MAX_CATS_INLINE)
  const inlineTableData = data.slice(0, MAX_CATS_INLINE)
  const truncatedChart = chartData.length > MAX_CATS_INLINE
  const truncatedTable = data.length > MAX_CATS_INLINE

  const toggleOptions = isSingleGroup
    ? [
        { label: 'Chart', value: 'chart' as const },
        { label: 'Table', value: 'table' as const },
      ]
    : [
        { label: 'Count', value: 'count' as const },
        { label: '%', value: 'percent' as const },
        { label: 'Table', value: 'table' as const },
      ]

  const openModalTable = () => { setModalTab('table'); setModalOpen(true) }

  const renderChart = () => {
    if (isSingleGroup) {
      return <PieChart data={chartData} categoryMap={categoryMap} groupId={activeGroups[0].id} countsByGroup={countsByGroup} onShowMore={openModalTable} />
    }
    const mode: StackMode = inlineView === 'percent' ? 'percent' : 'count'
    return (
      <>
        <StackedBarChart data={inlineChartData} activeGroups={activeGroups} width={280} margin={MARGIN} mode={mode} />
        {truncatedChart && (
          <Typography.Link style={{ fontSize: 10 }} onClick={openModalTable}>
            ...and {chartData.length - MAX_CATS_INLINE} more
          </Typography.Link>
        )}
      </>
    )
  }

  const renderModalChart = () => {
    if (isSingleGroup) {
      return <PieChart data={chartData} categoryMap={categoryMap} groupId={activeGroups[0].id} countsByGroup={countsByGroup} radius={140} />
    }
    const mode: StackMode = (inlineView === 'percent') ? 'percent' : 'count'
    return <VerticalStackedBarChart data={chartData} activeGroups={activeGroups} mode={mode} />
  }

  return (
    <div style={{ marginBottom: 12 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 4 }}>
        <Typography.Text strong style={{ fontSize: 12 }}>{name}</Typography.Text>
        <span style={{ display: 'flex', gap: 8, alignItems: 'center' }}>
          <Segmented
            size="small"
            value={inlineView}
            onChange={(v) => setInlineView(v as InlineView)}
            options={toggleOptions}
            style={{ fontSize: 11 }}
          />
          <ExpandOutlined
            style={{ fontSize: 11, cursor: 'pointer', color: '#1677ff' }}
            onClick={() => { setModalTab('chart'); setModalOpen(true) }}
          />
        </span>
      </div>

      {inlineView === 'table' ? (
        <>
          <CategoryTable data={inlineTableData} activeGroups={activeGroups} />
          {truncatedTable && (
            <Typography.Link style={{ fontSize: 10 }} onClick={openModalTable}>
              ...and {data.length - MAX_CATS_INLINE} more
            </Typography.Link>
          )}
        </>
      ) : (
        renderChart()
      )}

      <ChartModal
        title={name}
        open={modalOpen}
        onClose={() => setModalOpen(false)}
        chart={renderModalChart()}
        table={<CategoryTable data={data} activeGroups={activeGroups} />}
        defaultTab={modalTab}
      />
    </div>
  )
}
