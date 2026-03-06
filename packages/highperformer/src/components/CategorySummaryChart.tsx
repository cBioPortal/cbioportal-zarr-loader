import { useState } from 'react'
import { Group } from '@visx/group'
import { scaleBand, scaleLinear } from '@visx/scale'
import { Bar } from '@visx/shape'
import { AxisLeft } from '@visx/axis'
import { Typography } from 'antd'
import { ExpandOutlined } from '@ant-design/icons'
import type { SelectionGroup } from '../store/useAppStore'
import type { RGB } from '../utils/colors'
import ChartModal from './ChartModal'

interface CategorySummaryChartProps {
  name: string
  categoryMap: { label: string; color: RGB }[]
  countsByGroup: Map<number, Uint32Array>
  groups: SelectionGroup[]
  singleGroupId?: number
}

interface CategoryRow {
  label: string
  counts: Record<number, number>
  total: Record<number, number>
}

const MARGIN = { top: 4, right: 12, bottom: 4, left: 100 }
const MODAL_MARGIN = { top: 4, right: 16, bottom: 4, left: 140 }
const BAR_HEIGHT = 14
const BAR_GAP = 2
const MAX_CATS_INLINE = 15

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

function CategoryBarChart({ data, activeGroups, width, margin }: {
  data: CategoryRow[]
  activeGroups: SelectionGroup[]
  width: number
  margin: typeof MARGIN
}) {
  const numGroups = activeGroups.length
  const rowHeight = numGroups * BAR_HEIGHT + (numGroups - 1) * BAR_GAP + 4
  const height = data.length * rowHeight + margin.top + margin.bottom

  const maxCount = Math.max(
    ...data.flatMap((d) => activeGroups.map((g) => d.counts[g.id] ?? 0)),
    1,
  )

  const yScale = scaleBand({
    domain: data.map((d) => d.label),
    range: [margin.top, height - margin.bottom],
    padding: 0.15,
  })

  const xScale = scaleLinear({
    domain: [0, maxCount],
    range: [0, width - margin.left - margin.right],
  })

  const groupScale = scaleBand({
    domain: activeGroups.map((g) => String(g.id)),
    range: [0, yScale.bandwidth()],
    padding: 0.1,
  })

  return (
    <svg width={width} height={height}>
      <Group left={margin.left}>
        {data.map((d) => {
          const y0 = yScale(d.label) ?? 0
          return activeGroups.map((g) => {
            const barY = y0 + (groupScale(String(g.id)) ?? 0)
            const barWidth = xScale(d.counts[g.id] ?? 0)
            return (
              <Bar
                key={`${d.label}-${g.id}`}
                x={0}
                y={barY}
                width={barWidth}
                height={groupScale.bandwidth()}
                fill={`rgb(${g.color.join(',')})`}
                rx={2}
              />
            )
          })
        })}
        <AxisLeft
          scale={yScale}
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

function CategoryTable({ data, activeGroups }: { data: CategoryRow[]; activeGroups: SelectionGroup[] }) {
  return (
    <table style={{ width: '100%', fontSize: 12, borderCollapse: 'collapse' }}>
      <thead>
        <tr>
          <th style={{ textAlign: 'left', fontWeight: 600, padding: '4px 8px', borderBottom: '1px solid #f0f0f0' }}>Category</th>
          {activeGroups.map((g) => (
            <th key={`count-${g.id}`} style={{ textAlign: 'right', fontWeight: 600, padding: '4px 8px', borderBottom: '1px solid #f0f0f0', color: `rgb(${g.color.join(',')})` }}>
              G{g.id} count
            </th>
          ))}
          {activeGroups.map((g) => (
            <th key={`pct-${g.id}`} style={{ textAlign: 'right', fontWeight: 600, padding: '4px 8px', borderBottom: '1px solid #f0f0f0', color: `rgb(${g.color.join(',')})` }}>
              G{g.id} %
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

export default function CategorySummaryChart({ name, categoryMap, countsByGroup, groups, singleGroupId }: CategorySummaryChartProps) {
  const [showTable, setShowTable] = useState(false)
  const [modalOpen, setModalOpen] = useState(false)

  const activeGroups = singleGroupId != null
    ? groups.filter((g) => g.id === singleGroupId)
    : groups.filter((g) => countsByGroup.has(g.id))

  if (activeGroups.length === 0) return null

  const data = buildData(categoryMap, countsByGroup, activeGroups)
  if (data.length === 0) return null

  const inlineData = data.slice(0, MAX_CATS_INLINE)
  const truncated = data.length > MAX_CATS_INLINE

  const headerLinks = (
    <span style={{ display: 'flex', gap: 8, alignItems: 'center' }}>
      <Typography.Link style={{ fontSize: 11 }} onClick={() => setShowTable(!showTable)}>
        {showTable ? 'Chart' : 'Table'}
      </Typography.Link>
      <ExpandOutlined
        style={{ fontSize: 11, cursor: 'pointer', color: '#1677ff' }}
        onClick={() => setModalOpen(true)}
      />
    </span>
  )

  return (
    <div style={{ marginBottom: 12 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 4 }}>
        <Typography.Text strong style={{ fontSize: 12 }}>{name}</Typography.Text>
        {headerLinks}
      </div>

      {showTable ? (
        <>
          <CategoryTable data={inlineData} activeGroups={activeGroups} />
          {truncated && <Typography.Text type="secondary" style={{ fontSize: 10 }}>...and {data.length - MAX_CATS_INLINE} more</Typography.Text>}
        </>
      ) : (
        <>
          <CategoryBarChart data={inlineData} activeGroups={activeGroups} width={280} margin={MARGIN} />
          {truncated && <Typography.Text type="secondary" style={{ fontSize: 10 }}>...and {data.length - MAX_CATS_INLINE} more</Typography.Text>}
        </>
      )}

      <ChartModal
        title={name}
        open={modalOpen}
        onClose={() => setModalOpen(false)}
        chart={<CategoryBarChart data={data} activeGroups={activeGroups} width={640} margin={MODAL_MARGIN} />}
        table={<CategoryTable data={data} activeGroups={activeGroups} />}
      />
    </div>
  )
}
