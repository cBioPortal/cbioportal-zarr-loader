import { startTransition, useMemo, useState } from 'react'
import { Collapse, Segmented, Tooltip, Typography } from 'antd'
import { BarChartOutlined, CloseOutlined, InfoCircleOutlined, PieChartOutlined, SearchOutlined } from '@ant-design/icons'
import useAppStore from '../store/useAppStore'
import type { SelectionGroup } from '../store/useAppStore'
import GroupOverview from './GroupOverview'
import VariablePicker from './VariablePicker'
import ByVariableView from './ByVariableView'
import CategoryDotPlot from './CategoryDotPlot'
import ExpressionDotPlot from './ExpressionDotPlot'
import { ALL_CELLS_GROUP_ID } from '../constants'

const ALL_CELLS_COLOR: [number, number, number] = [120, 120, 120]

type SummaryContext = 'all' | 'selections' | 'compare'

const collapsedIconStyle: React.CSSProperties = {
  fontSize: 16,
  color: '#999',
  padding: '12px 0',
  display: 'flex',
  justifyContent: 'center',
  cursor: 'pointer',
}

interface SummaryPanelProps {
  collapsed?: boolean
  onExpand?: () => void
}

export default function SummaryPanel({ collapsed, onExpand }: SummaryPanelProps) {
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

  // Build groups with real indices for the active context
  const numPoints = embeddingData?.numPoints ?? 0

  const allCellsGroups = useMemo<SelectionGroup[]>(() => {
    if (numPoints === 0) return []
    const indices = new Uint32Array(numPoints)
    for (let i = 0; i < numPoints; i++) indices[i] = i
    return [{
      id: ALL_CELLS_GROUP_ID,
      polygon: [],
      type: 'rectangle' as const,
      indices,
      color: ALL_CELLS_COLOR,
    }]
  }, [numPoints])

  const activeSelectionGroups = useMemo(
    () => selectionGroups.filter((g) => g.indices.length > 0),
    [selectionGroups],
  )

  if (!summaryPanelOpen) return null

  // Collapsed view — icon strip
  if (collapsed) {
    const hasObs = summaryObsColumns.length > 0
    const hasGenes = summaryGenes.length > 0

    return (
      <div
        onClick={onExpand}
        style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          paddingTop: 12,
          height: '100%',
        }}
      >
        <Tooltip title="Add variables" placement="left">
          <div style={collapsedIconStyle}><SearchOutlined /></div>
        </Tooltip>
        {hasObs && (
          <Tooltip title="Category Distributions" placement="left">
            <div style={collapsedIconStyle}><PieChartOutlined /></div>
          </Tooltip>
        )}
        {hasGenes && (
          <Tooltip title="Expression Distributions" placement="left">
            <div style={collapsedIconStyle}><BarChartOutlined /></div>
          </Tooltip>
        )}
      </div>
    )
  }

  const hasGroups = selectionGroups.some((g) => g.indices.length > 0)
  const hasMultipleGroups = activeSelectionGroups.length >= 2
  const hasVariables = summaryObsColumns.length > 0 || summaryGenes.length > 0

  // Both contexts are always mounted — components are pure store readers,
  // no hooks dispatching workers. CSS display toggle for instant switching.
  const renderCharts = (groups: SelectionGroup[], visible: boolean) => {
    const style = visible ? undefined : { display: 'none' } as React.CSSProperties
    const collapseItems = []
    if (summaryObsColumns.length > 0) {
      collapseItems.push({
        key: 'obs',
        label: <Typography.Text strong style={{ fontSize: 12 }}>Category Distributions</Typography.Text>,
        children: <ByVariableView type="obs" groups={groups} onRemove={removeSummaryObsColumn} />,
      })
    }
    if (summaryGenes.length > 0) {
      collapseItems.push({
        key: 'genes',
        label: <Typography.Text strong style={{ fontSize: 12 }}>Expression Distributions</Typography.Text>,
        children: <ByVariableView type="genes" groups={groups} onRemove={removeSummaryGene} />,
      })
    }
    // Show expression dot plot when both categorical obs and genes are selected
    const hasCategoricalObs = summaryObsColumns.some((c) => summaryObsData.has(c))
    if (hasCategoricalObs && summaryGenes.length > 0) {
      collapseItems.push({
        key: 'exprcat',
        label: <Typography.Text strong style={{ fontSize: 12 }}>Gene Expression by Category</Typography.Text>,
        children: groups.length === 1
          ? <ExpressionDotPlot groups={groups} />
          : (
            <div style={{ display: 'flex', flexDirection: 'column', gap: 12 }}>
              {groups.map((g) => (
                <div key={g.id}>
                  <Typography.Text style={{ fontSize: 11, color: `rgb(${g.color.join(',')})`, fontWeight: 600, display: 'block', marginBottom: 4 }}>
                    Group {g.id}
                  </Typography.Text>
                  <ExpressionDotPlot groups={[g]} />
                </div>
              ))}
            </div>
          ),
      })
    }

    if (collapseItems.length === 0) return null
    return (
      <div style={style}>
        <Collapse
          defaultActiveKey={['obs', 'genes', 'exprcat']}
          ghost
          size="small"
          style={{ marginTop: 4 }}
          items={collapseItems}
        />
      </div>
    )
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
            onChange={(v) => startTransition(() => setContext(v as SummaryContext))}
            options={[
              { label: 'All Cells', value: 'all' },
              { label: `Selections${hasGroups ? '' : ' (none)'}`, value: 'selections', disabled: !hasGroups },
              { label: 'Compare', value: 'compare', disabled: !hasMultipleGroups },
            ]}
            style={{ fontSize: 11 }}
          />
        </div>

        <div style={{ display: context === 'all' ? 'block' : 'none' }}>
          <Typography.Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 4 }}>
            Showing summaries for all {numPoints.toLocaleString()} cells
          </Typography.Text>
        </div>

        {hasGroups && (
          <div style={{ display: context === 'selections' ? 'block' : 'none' }}>
            <GroupOverview groups={selectionGroups} totalCells={numPoints} />
          </div>
        )}

        {renderCharts(allCellsGroups, context === 'all')}
        {hasGroups && renderCharts(activeSelectionGroups, context === 'selections')}

        {hasMultipleGroups && (
          <div style={{ display: context === 'compare' ? 'block' : 'none' }}>
            <GroupOverview groups={selectionGroups} totalCells={numPoints} />
            <div style={{ display: 'flex', alignItems: 'center', gap: 6, marginBottom: 8 }}>
              <Typography.Text type="secondary" style={{ fontSize: 12 }}>
                Comparing category distributions across groups
              </Typography.Text>
              <Tooltip
                title="Each row is a category. Dots show the percentage of that category within each group. Lines connect the dots — longer lines mean bigger differences between groups."
                placement="bottom"
              >
                <InfoCircleOutlined style={{ fontSize: 11, color: '#999', cursor: 'help' }} />
              </Tooltip>
            </div>
            {summaryObsColumns.map((col) => {
              const catData = summaryObsData.get(col)
              if (!catData) return null
              return (
                <CategoryDotPlot
                  key={col}
                  name={col}
                  categoryMap={catData.categoryMap}
                  groups={activeSelectionGroups}
                  onRemove={() => removeSummaryObsColumn(col)}
                />
              )
            })}
          </div>
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
