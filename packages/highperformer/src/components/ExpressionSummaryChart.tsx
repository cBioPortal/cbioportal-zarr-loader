import { useState } from 'react'
import { Group } from '@visx/group'
import { scaleLinear } from '@visx/scale'
import { Bar } from '@visx/shape'
import { AxisBottom, AxisLeft } from '@visx/axis'
import { Typography } from 'antd'
import { ExpandOutlined } from '@ant-design/icons'
import type { SelectionGroup } from '../store/useAppStore'
import ChartModal from './ChartModal'

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
}

function formatNum(n: number): string {
  return Math.abs(n) >= 1000 ? n.toExponential(1) : n.toFixed(2)
}

function ExpressionHistogram({ statsByGroup, activeGroups, width, height }: {
  statsByGroup: Map<number, ExpressionStats>
  activeGroups: SelectionGroup[]
  width: number
  height: number
}) {
  const margin = { top: 8, right: 8, bottom: 24, left: 36 }
  const innerWidth = width - margin.left - margin.right
  const innerHeight = height - margin.top - margin.bottom

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

  const xScale = scaleLinear({ domain: [globalMin, globalMax], range: [0, innerWidth] })
  const yScale = scaleLinear({ domain: [0, maxBinCount], range: [innerHeight, 0], nice: true })

  const firstStats = statsByGroup.get(activeGroups[0].id)
  const numBins = firstStats?.bins.length ?? 30
  const binWidth = innerWidth / numBins

  return (
    <svg width={width} height={height}>
      <Group left={margin.left} top={margin.top}>
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
  )
}

function ExpressionStatsTable({ statsByGroup, activeGroups, fontSize = 10 }: {
  statsByGroup: Map<number, ExpressionStats>
  activeGroups: SelectionGroup[]
  fontSize?: number
}) {
  return (
    <table style={{ width: '100%', fontSize, borderCollapse: 'collapse', marginTop: 4 }}>
      <thead>
        <tr>
          <th style={{ textAlign: 'left', padding: '2px 8px', fontWeight: 600, borderBottom: '1px solid #f0f0f0' }}></th>
          {activeGroups.map((g) => (
            <th key={g.id} style={{ textAlign: 'right', padding: '2px 8px', fontWeight: 600, color: `rgb(${g.color.join(',')})`, borderBottom: '1px solid #f0f0f0' }}>
              G{g.id}
            </th>
          ))}
        </tr>
      </thead>
      <tbody>
        {(['mean', 'median', 'std', 'min', 'max'] as const).map((stat) => (
          <tr key={stat}>
            <td style={{ padding: '2px 8px', fontWeight: 500 }}>{stat}</td>
            {activeGroups.map((g) => {
              const stats = statsByGroup.get(g.id)
              return (
                <td key={g.id} style={{ textAlign: 'right', padding: '2px 8px' }}>
                  {stats ? formatNum(stats[stat]) : '—'}
                </td>
              )
            })}
          </tr>
        ))}
      </tbody>
    </table>
  )
}

export default function ExpressionSummaryChart({ name, statsByGroup, groups }: ExpressionSummaryChartProps) {
  const [modalOpen, setModalOpen] = useState(false)

  const activeGroups = groups.filter((g) => statsByGroup.has(g.id))

  if (activeGroups.length === 0) return null

  return (
    <div style={{ marginBottom: 12 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 4 }}>
        <Typography.Text strong style={{ fontSize: 12 }}>{name}</Typography.Text>
        <ExpandOutlined
          style={{ fontSize: 11, cursor: 'pointer', color: '#1677ff' }}
          onClick={() => setModalOpen(true)}
        />
      </div>

      <ExpressionHistogram statsByGroup={statsByGroup} activeGroups={activeGroups} width={260} height={100} />
      <ExpressionStatsTable statsByGroup={statsByGroup} activeGroups={activeGroups} />

      <ChartModal
        title={name}
        open={modalOpen}
        onClose={() => setModalOpen(false)}
        chart={
          <>
            <ExpressionHistogram statsByGroup={statsByGroup} activeGroups={activeGroups} width={660} height={300} />
            <ExpressionStatsTable statsByGroup={statsByGroup} activeGroups={activeGroups} fontSize={12} />
          </>
        }
        table={<ExpressionStatsTable statsByGroup={statsByGroup} activeGroups={activeGroups} fontSize={12} />}
      />
    </div>
  )
}
