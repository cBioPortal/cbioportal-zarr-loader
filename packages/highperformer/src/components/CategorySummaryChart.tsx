import { useState } from 'react'
import { Group } from '@visx/group'
import { scaleBand, scaleLinear } from '@visx/scale'
import { Bar } from '@visx/shape'
import { AxisLeft } from '@visx/axis'
import { Typography } from 'antd'
import type { SelectionGroup } from '../store/useAppStore'
import type { RGB } from '../utils/colors'

interface CategorySummaryChartProps {
  name: string
  categoryMap: { label: string; color: RGB }[]
  countsByGroup: Map<number, Uint32Array>
  groups: SelectionGroup[]
  singleGroupId?: number
}

const MARGIN = { top: 4, right: 12, bottom: 4, left: 100 }
const BAR_HEIGHT = 14
const BAR_GAP = 2

export default function CategorySummaryChart({ name, categoryMap, countsByGroup, groups, singleGroupId }: CategorySummaryChartProps) {
  const [showTable, setShowTable] = useState(false)

  const activeGroups = singleGroupId != null
    ? groups.filter((g) => g.id === singleGroupId)
    : groups.filter((g) => countsByGroup.has(g.id))

  if (activeGroups.length === 0) return null

  // Build data: only categories with non-zero counts across any group
  const data: { label: string; counts: Record<number, number>; total: Record<number, number> }[] = []
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

  if (data.length === 0) return null

  // Sort by total count descending (sum across groups)
  data.sort((a, b) => {
    const sumA = Object.values(a.counts).reduce((x, y) => x + y, 0)
    const sumB = Object.values(b.counts).reduce((x, y) => x + y, 0)
    return sumB - sumA
  })

  const MAX_CATS = 15
  const displayData = data.slice(0, MAX_CATS)
  const truncated = data.length > MAX_CATS

  if (showTable) {
    return (
      <div style={{ marginBottom: 12 }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 4 }}>
          <Typography.Text strong style={{ fontSize: 12 }}>{name}</Typography.Text>
          <Typography.Link style={{ fontSize: 11 }} onClick={() => setShowTable(false)}>Chart</Typography.Link>
        </div>
        <table style={{ width: '100%', fontSize: 11, borderCollapse: 'collapse' }}>
          <thead>
            <tr>
              <th style={{ textAlign: 'left', fontWeight: 600, padding: '2px 4px', borderBottom: '1px solid #f0f0f0' }}>Category</th>
              {activeGroups.map((g) => (
                <th key={g.id} style={{ textAlign: 'right', fontWeight: 600, padding: '2px 4px', borderBottom: '1px solid #f0f0f0', color: `rgb(${g.color.join(',')})` }}>
                  G{g.id}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {displayData.map((row) => (
              <tr key={row.label}>
                <td style={{ padding: '1px 4px', maxWidth: 120, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>{row.label}</td>
                {activeGroups.map((g) => {
                  const count = row.counts[g.id] ?? 0
                  const total = row.total[g.id] ?? 1
                  const pct = total > 0 ? ((count / total) * 100).toFixed(1) : '0.0'
                  return (
                    <td key={g.id} style={{ textAlign: 'right', padding: '1px 4px' }}>{pct}%</td>
                  )
                })}
              </tr>
            ))}
          </tbody>
        </table>
        {truncated && <Typography.Text type="secondary" style={{ fontSize: 10 }}>...and {data.length - MAX_CATS} more</Typography.Text>}
      </div>
    )
  }

  const numGroups = activeGroups.length
  const rowHeight = numGroups * BAR_HEIGHT + (numGroups - 1) * BAR_GAP + 4
  const height = displayData.length * rowHeight + MARGIN.top + MARGIN.bottom
  const width = 280

  const maxCount = Math.max(
    ...displayData.flatMap((d) => activeGroups.map((g) => d.counts[g.id] ?? 0)),
    1,
  )

  const yScale = scaleBand({
    domain: displayData.map((d) => d.label),
    range: [MARGIN.top, height - MARGIN.bottom],
    padding: 0.15,
  })

  const xScale = scaleLinear({
    domain: [0, maxCount],
    range: [0, width - MARGIN.left - MARGIN.right],
  })

  const groupScale = scaleBand({
    domain: activeGroups.map((g) => String(g.id)),
    range: [0, yScale.bandwidth()],
    padding: 0.1,
  })

  return (
    <div style={{ marginBottom: 12 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 4 }}>
        <Typography.Text strong style={{ fontSize: 12 }}>{name}</Typography.Text>
        <Typography.Link style={{ fontSize: 11 }} onClick={() => setShowTable(true)}>Table</Typography.Link>
      </div>
      <svg width={width} height={height}>
        <Group left={MARGIN.left}>
          {displayData.map((d) => {
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
              return s.length > 14 ? s.slice(0, 12) + '\u2026' : s
            }}
          />
        </Group>
      </svg>
      {truncated && <Typography.Text type="secondary" style={{ fontSize: 10 }}>...and {data.length - MAX_CATS} more</Typography.Text>}
    </div>
  )
}
