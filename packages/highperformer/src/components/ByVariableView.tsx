import useAppStore from '../store/useAppStore'
import CategorySummaryChart from './CategorySummaryChart'
import ExpressionSummaryChart from './ExpressionSummaryChart'
import type { SummaryResult } from '../hooks/useSummaryData'

interface ByVariableViewProps {
  results: SummaryResult[]
}

export default function ByVariableView({ results }: ByVariableViewProps) {
  const selectionGroups = useAppStore((s) => s.selectionGroups)
  const groups = selectionGroups.filter((g) => g.indices.length > 0)

  if (results.length === 0) return null

  return (
    <div>
      {results.map((r) => {
        if (r.type === 'category') {
          return (
            <CategorySummaryChart
              key={r.name}
              name={r.name}
              categoryMap={r.categoryMap}
              countsByGroup={r.countsByGroup}
              groups={groups}
            />
          )
        }
        return (
          <ExpressionSummaryChart
            key={r.name}
            name={r.name}
            statsByGroup={r.statsByGroup}
            groups={groups}
          />
        )
      })}
    </div>
  )
}
