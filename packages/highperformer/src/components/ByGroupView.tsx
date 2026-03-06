import { Collapse, Typography } from 'antd'
import useAppStore from '../store/useAppStore'
import CategorySummaryChart from './CategorySummaryChart'
import ExpressionSummaryChart from './ExpressionSummaryChart'
import type { SummaryResult } from '../hooks/useSummaryData'

interface ByGroupViewProps {
  results: SummaryResult[]
}

export default function ByGroupView({ results }: ByGroupViewProps) {
  const selectionGroups = useAppStore((s) => s.selectionGroups)
  const groups = selectionGroups.filter((g) => g.indices.length > 0)

  if (results.length === 0 || groups.length === 0) return null

  const items = groups.map((g) => ({
    key: String(g.id),
    label: (
      <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
        <div style={{
          width: 10,
          height: 10,
          borderRadius: '50%',
          backgroundColor: `rgb(${g.color.join(',')})`,
        }} />
        <Typography.Text strong style={{ fontSize: 12 }}>
          Group {g.id}
        </Typography.Text>
        <Typography.Text type="secondary" style={{ fontSize: 11 }}>
          ({g.indices.length.toLocaleString()} cells)
        </Typography.Text>
      </div>
    ),
    children: (
      <div>
        {results.map((r) => {
          if (r.type === 'category') {
            return (
              <CategorySummaryChart
                key={r.name}
                name={r.name}
                categoryMap={r.categoryMap}
                countsByGroup={r.countsByGroup}
                groups={[g]}
                singleGroupId={g.id}
              />
            )
          }
          return (
            <ExpressionSummaryChart
              key={r.name}
              name={r.name}
              statsByGroup={r.statsByGroup}
              groups={[g]}
              singleGroupId={g.id}
            />
          )
        })}
      </div>
    ),
  }))

  return (
    <Collapse
      ghost
      size="small"
      defaultActiveKey={groups.map((g) => String(g.id))}
      items={items}
    />
  )
}
