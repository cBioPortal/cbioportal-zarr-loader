import { useMemo, useState } from 'react'
import { Collapse, Segmented, Typography } from 'antd'
import { CloseOutlined } from '@ant-design/icons'
import useAppStore from '../store/useAppStore'
import type { SelectionGroup } from '../store/useAppStore'
import GroupOverview from './GroupOverview'
import VariablePicker from './VariablePicker'
import { useSummaryData } from '../hooks/useSummaryData'
import { useAllCellsSummary, ALL_CELLS_GROUP_ID } from '../hooks/useAllCellsSummary'
import ByVariableView from './ByVariableView'

const ALL_CELLS_COLOR: [number, number, number] = [120, 120, 120]

type SummaryContext = 'all' | 'selections'

export default function SummaryPanel() {
  const summaryPanelOpen = useAppStore((s) => s.summaryPanelOpen)
  const setSummaryPanelOpen = useAppStore((s) => s.setSummaryPanelOpen)
  const selectionGroups = useAppStore((s) => s.selectionGroups)
  const summaryObsColumns = useAppStore((s) => s.summaryObsColumns)
  const summaryGenes = useAppStore((s) => s.summaryGenes)
  const obsColumnNames = useAppStore((s) => s.obsColumnNames)
  const varNames = useAppStore((s) => s.varNames)
  const geneLabelMap = useAppStore((s) => s.geneLabelMap)
  const summaryObsData = useAppStore((s) => s.summaryObsData)
  const summaryGeneData = useAppStore((s) => s.summaryGeneData)
  const summaryObsContinuousData = useAppStore((s) => s.summaryObsContinuousData)
  const embeddingData = useAppStore((s) => s.embeddingData)
  const addSummaryObsColumn = useAppStore((s) => s.addSummaryObsColumn)
  const removeSummaryObsColumn = useAppStore((s) => s.removeSummaryObsColumn)
  const addSummaryGene = useAppStore((s) => s.addSummaryGene)
  const removeSummaryGene = useAppStore((s) => s.removeSummaryGene)

  const [context, setContext] = useState<SummaryContext>('all')

  const groupResults = useSummaryData()
  const allCellsResults = useAllCellsSummary()

  const allCellsGroup = useMemo<SelectionGroup[]>(() => [{
    id: ALL_CELLS_GROUP_ID,
    polygon: [],
    type: 'rectangle',
    indices: new Uint32Array(0),
    color: ALL_CELLS_COLOR,
  }], [])

  if (!summaryPanelOpen) return null

  const hasGroups = selectionGroups.some((g) => g.indices.length > 0)
  const hasVariables = summaryObsColumns.length > 0 || summaryGenes.length > 0

  // Pick results based on context
  const activeResults = context === 'all' ? allCellsResults : groupResults
  const activeGroups = context === 'all' ? allCellsGroup : undefined

  const obsResults = activeResults.filter((r) => r.type === 'category' || summaryObsColumns.includes(r.name))
  const geneResults = activeResults.filter((r) => r.type === 'expression' && !summaryObsColumns.includes(r.name))

  const collapseItems = []
  if (obsResults.length > 0) {
    collapseItems.push({
      key: 'obs',
      label: <Typography.Text strong style={{ fontSize: 12 }}>Obs Summaries</Typography.Text>,
      children: <ByVariableView results={obsResults} groups={activeGroups} />,
    })
  }
  if (geneResults.length > 0) {
    collapseItems.push({
      key: 'genes',
      label: <Typography.Text strong style={{ fontSize: 12 }}>Gene Summaries</Typography.Text>,
      children: <ByVariableView results={geneResults} groups={activeGroups} />,
    })
  }

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
      <div style={{
        padding: '12px 16px',
        borderBottom: '1px solid #f0f0f0',
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center',
      }}>
        <Typography.Text strong style={{ fontSize: 14 }}>Summary</Typography.Text>
        <CloseOutlined
          style={{ fontSize: 12, cursor: 'pointer', color: '#999' }}
          onClick={() => setSummaryPanelOpen(false)}
        />
      </div>

      <div style={{ flex: 1, overflow: 'auto', padding: '12px 16px' }}>
        <VariablePicker
          label="Obs Columns"
          options={obsColumnNames}
          selected={summaryObsColumns}
          onAdd={addSummaryObsColumn}
          onRemove={removeSummaryObsColumn}
          loading={new Set(summaryObsColumns.filter((c) => !summaryObsData.has(c) && !summaryObsContinuousData.has(c)))}
        />

        <VariablePicker
          label="Genes"
          variant="search"
          options={varNames}
          selected={summaryGenes}
          onAdd={addSummaryGene}
          onRemove={removeSummaryGene}
          labelMap={geneLabelMap}
          loading={new Set(summaryGenes.filter((g) => !summaryGeneData.has(g)))}
        />

        <div style={{ display: 'flex', justifyContent: 'center', margin: '8px 0' }}>
          <Segmented
            size="small"
            value={context}
            onChange={(v) => setContext(v as SummaryContext)}
            options={[
              { label: 'All Cells', value: 'all' },
              { label: `Selections${hasGroups ? '' : ' (none)'}`, value: 'selections', disabled: !hasGroups },
            ]}
            style={{ fontSize: 11 }}
          />
        </div>

        {context === 'selections' && hasGroups && (
          <GroupOverview groups={selectionGroups} totalCells={embeddingData?.numPoints ?? 0} />
        )}

        {collapseItems.length > 0 && (
          <Collapse
            defaultActiveKey={['obs', 'genes']}
            ghost
            size="small"
            style={{ marginTop: 4 }}
            items={collapseItems}
          />
        )}

        {!hasVariables && (
          <Typography.Text type="secondary" style={{ fontSize: 12 }}>
            Add obs columns or genes above to see summaries.
          </Typography.Text>
        )}
      </div>
    </div>
  )
}
