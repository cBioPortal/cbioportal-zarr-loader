import { useMemo } from 'react'
import { Typography } from 'antd'
import type { SelectionGroup } from '../store/useAppStore'
import { ALL_CELLS_GROUP_ID } from '../hooks/useAllCellsSummary'

interface GroupOverviewProps {
  groups: SelectionGroup[]
  totalCells: number
}

interface OverlapStats {
  overlapCount: number
  unionCount: number
  uniqueCounts: Map<number, number>
  pairwiseOverlaps: Map<string, number>
}

function computeOverlap(groups: SelectionGroup[]): OverlapStats {
  const uniqueCounts = new Map<number, number>()
  const pairwiseOverlaps = new Map<string, number>()

  if (groups.length <= 1) {
    const count = groups[0]?.indices.length ?? 0
    if (groups[0]) uniqueCounts.set(groups[0].id, count)
    return { overlapCount: 0, unionCount: count, uniqueCounts, pairwiseOverlaps }
  }

  const membership = new Map<number, number[]>()
  for (const g of groups) {
    for (let i = 0; i < g.indices.length; i++) {
      const idx = g.indices[i]
      const existing = membership.get(idx)
      if (existing) {
        existing.push(g.id)
      } else {
        membership.set(idx, [g.id])
      }
    }
  }

  let overlapCount = 0
  for (const g of groups) uniqueCounts.set(g.id, 0)

  for (const groupIds of membership.values()) {
    if (groupIds.length > 1) {
      overlapCount++
      for (let i = 0; i < groupIds.length; i++) {
        for (let j = i + 1; j < groupIds.length; j++) {
          const key = `${groupIds[i]}-${groupIds[j]}`
          pairwiseOverlaps.set(key, (pairwiseOverlaps.get(key) ?? 0) + 1)
        }
      }
    } else {
      uniqueCounts.set(groupIds[0], (uniqueCounts.get(groupIds[0]) ?? 0) + 1)
    }
  }

  return { overlapCount, unionCount: membership.size, uniqueCounts, pairwiseOverlaps }
}

function circleOverlapArea(r1: number, r2: number, d: number): number {
  if (d >= r1 + r2) return 0
  if (d + r2 <= r1) return Math.PI * r2 * r2
  if (d + r1 <= r2) return Math.PI * r1 * r1

  const a = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
  const h = Math.sqrt(Math.max(0, r1 * r1 - a * a))
  const angle1 = 2 * Math.atan2(h, a)
  const angle2 = 2 * Math.atan2(h, d - a)
  return 0.5 * r1 * r1 * (angle1 - Math.sin(angle1)) + 0.5 * r2 * r2 * (angle2 - Math.sin(angle2))
}

function findDistance(r1: number, r2: number, targetArea: number): number {
  if (targetArea <= 0) return r1 + r2
  const maxArea = Math.PI * Math.min(r1, r2) * Math.min(r1, r2)
  if (targetArea >= maxArea) return Math.abs(r1 - r2)

  let lo = Math.abs(r1 - r2)
  let hi = r1 + r2
  for (let i = 0; i < 50; i++) {
    const mid = (lo + hi) / 2
    if (circleOverlapArea(r1, r2, mid) > targetArea) lo = mid
    else hi = mid
  }
  return (lo + hi) / 2
}

const SVG_WIDTH = 220
const SVG_HEIGHT = 140
const MAX_RADIUS = 50
const MIN_RADIUS = 18

function VennDiagram({ groups, stats, totalCells }: {
  groups: SelectionGroup[]
  stats: OverlapStats
  totalCells: number
}) {
  const cx = SVG_WIDTH / 2
  const cy = SVG_HEIGHT / 2

  const counts = groups.map((g) => g.indices.length)
  const maxCount = Math.max(...counts)
  const radii = counts.map((c) =>
    maxCount > 0 ? MIN_RADIUS + (MAX_RADIUS - MIN_RADIUS) * Math.sqrt(c / maxCount) : MIN_RADIUS
  )

  let circles: { x: number; y: number; r: number }[]

  if (groups.length === 1) {
    circles = [{ x: cx, y: cy, r: radii[0] }]
  } else if (groups.length === 2) {
    const overlapKey = `${groups[0].id}-${groups[1].id}`
    const overlapCount = stats.pairwiseOverlaps.get(overlapKey) ?? 0
    const totalArea1 = Math.PI * radii[0] * radii[0]
    const totalArea2 = Math.PI * radii[1] * radii[1]
    const overlapFraction = counts[0] > 0 && counts[1] > 0
      ? overlapCount / Math.min(counts[0], counts[1])
      : 0
    const targetArea = overlapFraction * Math.min(totalArea1, totalArea2)
    const dist = findDistance(radii[0], radii[1], targetArea)
    circles = [
      { x: cx - dist / 2, y: cy, r: radii[0] },
      { x: cx + dist / 2, y: cy, r: radii[1] },
    ]
  } else {
    // 3 groups: triangulate
    const pairs: [number, number][] = [[0, 1], [0, 2], [1, 2]]
    const distances = pairs.map(([i, j]) => {
      const key1 = `${groups[i].id}-${groups[j].id}`
      const key2 = `${groups[j].id}-${groups[i].id}`
      const oc = stats.pairwiseOverlaps.get(key1) ?? stats.pairwiseOverlaps.get(key2) ?? 0
      const ai = Math.PI * radii[i] * radii[i]
      const aj = Math.PI * radii[j] * radii[j]
      const frac = counts[i] > 0 && counts[j] > 0 ? oc / Math.min(counts[i], counts[j]) : 0
      return findDistance(radii[i], radii[j], frac * Math.min(ai, aj))
    })

    const d01 = distances[0]
    const d02 = distances[1]
    const d12 = distances[2]
    const cosA = d01 > 0 ? (d02 * d02 + d01 * d01 - d12 * d12) / (2 * d02 * d01) : 0
    const clamped = Math.max(-1, Math.min(1, cosA))
    const sinA = Math.sqrt(1 - clamped * clamped)
    const raw: [number, number][] = [
      [0, 0],
      [d01, 0],
      [d02 * clamped, -d02 * sinA],
    ]

    // Scale and center
    const allX = raw.map((p, i) => [p[0] - radii[i], p[0] + radii[i]]).flat()
    const allY = raw.map((p, i) => [p[1] - radii[i], p[1] + radii[i]]).flat()
    const bx0 = Math.min(...allX), bx1 = Math.max(...allX)
    const by0 = Math.min(...allY), by1 = Math.max(...allY)
    const bw = bx1 - bx0, bh = by1 - by0
    const scale = Math.min((SVG_WIDTH - 20) / bw, (SVG_HEIGHT - 20) / bh, 1)
    const ox = cx - ((bx0 + bx1) / 2) * scale
    const oy = cy - ((by0 + by1) / 2) * scale

    circles = raw.map((p, i) => ({ x: p[0] * scale + ox, y: p[1] * scale + oy, r: radii[i] * scale }))
  }

  const colors = groups.map((g) => `rgb(${g.color.join(',')})`)

  return (
    <div>
      <div style={{ display: 'flex', justifyContent: 'center' }}>
        <svg width={SVG_WIDTH} height={SVG_HEIGHT}>
          {circles.map((c, i) => (
            <circle
              key={groups[i].id}
              cx={c.x} cy={c.y} r={c.r}
              fill={colors[i]} fillOpacity={0.3}
              stroke={colors[i]} strokeWidth={1.5}
            />
          ))}
          {/* Unique count labels — pushed outward from centroid */}
          {circles.map((c, i) => {
            const unique = stats.uniqueCounts.get(groups[i].id) ?? 0
            const centroidX = circles.reduce((s, cc) => s + cc.x, 0) / circles.length
            const centroidY = circles.reduce((s, cc) => s + cc.y, 0) / circles.length
            const dx = c.x - centroidX
            const dy = c.y - centroidY
            const len = Math.sqrt(dx * dx + dy * dy) || 1
            const push = groups.length === 1 ? 0 : c.r * 0.4
            const lx = c.x + (dx / len) * push
            const ly = c.y + (dy / len) * push
            return (
              <text key={`u-${groups[i].id}`} x={lx} y={ly} textAnchor="middle" dy="0.35em" fontSize={10} fill="#333">
                {unique.toLocaleString()}
              </text>
            )
          })}
          {/* Overlap count at centroid */}
          {groups.length > 1 && stats.overlapCount > 0 && (
            <text
              x={circles.reduce((s, c) => s + c.x, 0) / circles.length}
              y={circles.reduce((s, c) => s + c.y, 0) / circles.length}
              textAnchor="middle" dy="0.35em" fontSize={10} fontWeight={600} fill="#333"
            >
              {stats.overlapCount.toLocaleString()}
            </text>
          )}
        </svg>
      </div>
      {/* Stats below diagram */}
      <div style={{ display: 'flex', flexDirection: 'column', gap: 2, fontSize: 11, color: '#666', marginTop: 4 }}>
        {groups.map((g) => {
          const count = g.indices.length
          const pct = totalCells > 0 ? ((count / totalCells) * 100).toFixed(1) : '0.0'
          return (
            <div key={g.id} style={{ display: 'flex', justifyContent: 'space-between' }}>
              <span style={{ color: `rgb(${g.color.join(',')})`, fontWeight: 600 }}>
                {g.id === ALL_CELLS_GROUP_ID ? 'All Cells' : `Group ${g.id}`}
              </span>
              <span>
                <span style={{ fontWeight: 500, color: '#333' }}>{count.toLocaleString()}</span>
                {' '}cells ({pct}%)
              </span>
            </div>
          )
        })}
        {groups.length > 1 && (
          <>
            <div style={{ borderTop: '1px solid #f0f0f0', paddingTop: 2, display: 'flex', justifyContent: 'space-between' }}>
              <span>Overlap</span>
              <span>
                <span style={{ fontWeight: 500, color: '#333' }}>{stats.overlapCount.toLocaleString()}</span>
                {' '}cells
              </span>
            </div>
            {stats.unionCount > 0 && (
              <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                <span>Jaccard</span>
                <span style={{ fontWeight: 500, color: '#333' }}>
                  {(stats.overlapCount / stats.unionCount * 100).toFixed(1)}%
                </span>
              </div>
            )}
          </>
        )}
      </div>
    </div>
  )
}

export default function GroupOverview({ groups, totalCells }: GroupOverviewProps) {
  const activeGroups = useMemo(
    () => groups.filter((g) => g.indices.length > 0),
    [groups],
  )

  const stats = useMemo(
    () => computeOverlap(activeGroups),
    [activeGroups],
  )

  if (activeGroups.length === 0) return null

  return (
    <div style={{ marginBottom: 12 }}>
      <Typography.Text strong style={{ fontSize: 12 }}>Groups</Typography.Text>
      <VennDiagram groups={activeGroups} stats={stats} totalCells={totalCells} />
    </div>
  )
}
