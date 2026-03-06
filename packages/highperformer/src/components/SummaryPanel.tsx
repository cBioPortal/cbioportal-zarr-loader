import { Collapse, Typography } from 'antd'
import { CloseOutlined } from '@ant-design/icons'
import useAppStore from '../store/useAppStore'
import VariablePicker from './VariablePicker'
import { useSummaryData } from '../hooks/useSummaryData'
import ByVariableView from './ByVariableView'

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
  const addSummaryObsColumn = useAppStore((s) => s.addSummaryObsColumn)
  const removeSummaryObsColumn = useAppStore((s) => s.removeSummaryObsColumn)
  const addSummaryGene = useAppStore((s) => s.addSummaryGene)
  const removeSummaryGene = useAppStore((s) => s.removeSummaryGene)

  const results = useSummaryData()

  if (!summaryPanelOpen) return null

  const hasGroups = selectionGroups.some((g) => g.indices.length > 0)
  const obsResults = results.filter((r) => r.type === 'category' || summaryObsColumns.includes(r.name))
  const geneResults = results.filter((r) => r.type === 'expression' && !summaryObsColumns.includes(r.name))

  const collapseItems = []
  if (obsResults.length > 0) {
    collapseItems.push({
      key: 'obs',
      label: <Typography.Text strong style={{ fontSize: 12 }}>Obs Summaries</Typography.Text>,
      children: <ByVariableView results={obsResults} />,
    })
  }
  if (geneResults.length > 0) {
    collapseItems.push({
      key: 'genes',
      label: <Typography.Text strong style={{ fontSize: 12 }}>Gene Summaries</Typography.Text>,
      children: <ByVariableView results={geneResults} />,
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
        {!hasGroups ? (
          <Typography.Text type="secondary" style={{ fontSize: 12 }}>
            Draw a selection to see summaries.
          </Typography.Text>
        ) : (
          <>
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

            {collapseItems.length > 0 && (
              <Collapse
                defaultActiveKey={['obs', 'genes']}
                ghost
                size="small"
                style={{ marginTop: 12 }}
                items={collapseItems}
              />
            )}
          </>
        )}
      </div>
    </div>
  )
}
