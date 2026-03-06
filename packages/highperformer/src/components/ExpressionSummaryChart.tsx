import { Group } from '@visx/group'
import { scaleLinear } from '@visx/scale'
import { Bar } from '@visx/shape'
import { AxisBottom, AxisLeft } from '@visx/axis'
import { Typography } from 'antd'
import type { SelectionGroup } from '../store/useAppStore'

interface ExpressionStats {
  mean: number
  median: number
  std: number
  min: number
  max: number
  bins: Uint32Array
  binEdges: Float32Array
}

interface ExpressionSummaryChartProps {
  name: string
  statsByGroup: Map<number, ExpressionStats>
  groups: SelectionGroup[]
  singleGroupId?: number
}

const MARGIN = { top: 8, right: 8, bottom: 24, left: 36 }
const CHART_HEIGHT = 100
const CHART_WIDTH = 260

function formatNum(n: number): string {
  return Math.abs(n) >= 1000 ? n.toExponential(1) : n.toFixed(2)
}

export default function ExpressionSummaryChart({ name, statsByGroup, groups, singleGroupId }: ExpressionSummaryChartProps) {
  const activeGroups = singleGroupId != null
    ? groups.filter((g) => g.id === singleGroupId)
    : groups.filter((g) => statsByGroup.has(g.id))

  if (activeGroups.length === 0) return null

  const innerWidth = CHART_WIDTH - MARGIN.left - MARGIN.right
  const innerHeight = CHART_HEIGHT - MARGIN.top - MARGIN.bottom

  let globalMin = Infinity
  let globalMax = -Infinity
  let maxBinCount = 0
  for (const g of activeGroups) {
    const stats = statsByGroup.get(g.id)
    if (!stats) continue
    if (stats.min < globalMin) globalMin = stats.min
    if (stats.max > globalMax) globalMax = stats.max
    for (let i = 0; i < stats.bins.length; i++) {
      if (stats.bins[i] > maxBinCount) maxBinCount = stats.bins[i]
    }
  }

  const xScale = scaleLinear({
    domain: [globalMin, globalMax],
    range: [0, innerWidth],
  })

  const yScale = scaleLinear({
    domain: [0, maxBinCount],
    range: [innerHeight, 0],
    nice: true,
  })

  const firstStats = statsByGroup.get(activeGroups[0].id)
  const numBins = firstStats?.bins.length ?? 30
  const binWidth = innerWidth / numBins

  return (
    <div style={{ marginBottom: 12 }}>
      <Typography.Text strong style={{ fontSize: 12, display: 'block', marginBottom: 4 }}>{name}</Typography.Text>

      <svg width={CHART_WIDTH} height={CHART_HEIGHT}>
        <Group left={MARGIN.left} top={MARGIN.top}>
          {activeGroups.map((g) => {
            const stats = statsByGroup.get(g.id)
            if (!stats) return null
            const opacity = activeGroups.length > 1 ? 0.5 : 0.8
            return Array.from(stats.bins).map((count, bi) => (
              <Bar
                key={`${g.id}-${bi}`}
                x={bi * binWidth}
                y={yScale(count)}
                width={Math.max(binWidth - 1, 1)}
                height={innerHeight - yScale(count)}
                fill={`rgba(${g.color.join(',')}, ${opacity})`}
              />
            ))
          })}
          <AxisLeft
            scale={yScale}
            numTicks={3}
            stroke="#e8e8e8"
            tickStroke="#e8e8e8"
            tickLabelProps={{ fill: '#999', fontSize: 9 }}
          />
          <AxisBottom
            scale={xScale}
            top={innerHeight}
            numTicks={4}
            stroke="#e8e8e8"
            tickStroke="#e8e8e8"
            tickLabelProps={{ fill: '#999', fontSize: 9 }}
            tickFormat={(v) => formatNum(v as number)}
          />
        </Group>
      </svg>

      <table style={{ width: '100%', fontSize: 10, borderCollapse: 'collapse', marginTop: 4 }}>
        <thead>
          <tr>
            <th style={{ textAlign: 'left', padding: '1px 4px', fontWeight: 600, borderBottom: '1px solid #f0f0f0' }}></th>
            {activeGroups.map((g) => (
              <th key={g.id} style={{ textAlign: 'right', padding: '1px 4px', fontWeight: 600, color: `rgb(${g.color.join(',')})`, borderBottom: '1px solid #f0f0f0' }}>
                G{g.id}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {(['mean', 'median', 'std'] as const).map((stat) => (
            <tr key={stat}>
              <td style={{ padding: '1px 4px', fontWeight: 500 }}>{stat}</td>
              {activeGroups.map((g) => {
                const stats = statsByGroup.get(g.id)
                return (
                  <td key={g.id} style={{ textAlign: 'right', padding: '1px 4px' }}>
                    {stats ? formatNum(stats[stat]) : '—'}
                  </td>
                )
              })}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}
