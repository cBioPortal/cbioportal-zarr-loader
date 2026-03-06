import { useEffect, useRef, useState } from 'react'
import useAppStore, { getPool } from '../store/useAppStore'
import type { SelectionGroup } from '../store/useAppStore'
import type { CategorySummaryResponse, ExpressionSummaryResponse } from '../workers/summary.schemas'
import type { RGB } from '../utils/colors'

export interface CategorySummaryResult {
  type: 'category'
  name: string
  categoryMap: { label: string; color: RGB }[]
  countsByGroup: Map<number, Uint32Array>
}

export interface ExpressionSummaryResult {
  type: 'expression'
  name: string
  statsByGroup: Map<number, {
    mean: number
    median: number
    std: number
    min: number
    max: number
    bins: Uint32Array
    binEdges: Float32Array
  }>
}

export type SummaryResult = CategorySummaryResult | ExpressionSummaryResult

const NUM_BINS = 30

export function useSummaryData(): SummaryResult[] {
  const selectionGroups = useAppStore((s) => s.selectionGroups)
  const pinnedObsColumns = useAppStore((s) => s.pinnedObsColumns)
  const pinnedGenes = useAppStore((s) => s.pinnedGenes)
  const pinnedObsData = useAppStore((s) => s.pinnedObsData)
  const pinnedGeneData = useAppStore((s) => s.pinnedGeneData)

  const colorMode = useAppStore((s) => s.colorMode)
  const selectedObsColumn = useAppStore((s) => s.selectedObsColumn)
  const selectedGene = useAppStore((s) => s.selectedGene)
  const _categoryCodes = useAppStore((s) => s._categoryCodes)
  const categoryMap = useAppStore((s) => s.categoryMap)
  const _expressionData = useAppStore((s) => s._expressionData)
  const geneLabelMap = useAppStore((s) => s.geneLabelMap)

  const [results, setResults] = useState<SummaryResult[]>([])
  const versionRef = useRef(0)

  useEffect(() => {
    const groupsWithIndices = selectionGroups.filter((g) => g.indices.length > 0)
    if (groupsWithIndices.length === 0) {
      setResults([])
      return
    }

    versionRef.current++
    const version = versionRef.current

    const tasks: Promise<SummaryResult>[] = []

    // Active coloring summary
    if (colorMode === 'category' && selectedObsColumn && _categoryCodes && categoryMap.length > 0) {
      tasks.push(computeCategorySummary(selectedObsColumn, _categoryCodes, categoryMap, groupsWithIndices, version, versionRef))
    } else if (colorMode === 'gene' && selectedGene && _expressionData) {
      const displayName = geneLabelMap?.get(selectedGene) ?? selectedGene
      tasks.push(computeExpressionSummary(displayName, _expressionData, groupsWithIndices, version, versionRef))
    }

    // Pinned obs columns
    for (const name of pinnedObsColumns) {
      const data = pinnedObsData.get(name)
      if (!data) continue
      if (colorMode === 'category' && name === selectedObsColumn) continue
      tasks.push(computeCategorySummary(name, data.codes, data.categoryMap, groupsWithIndices, version, versionRef))
    }

    // Pinned genes
    for (const name of pinnedGenes) {
      const data = pinnedGeneData.get(name)
      if (!data) continue
      if (colorMode === 'gene' && name === selectedGene) continue
      const displayName = geneLabelMap?.get(name) ?? name
      tasks.push(computeExpressionSummary(displayName, data, groupsWithIndices, version, versionRef))
    }

    Promise.all(tasks).then((resolved) => {
      if (versionRef.current !== version) return
      setResults(resolved)
    })
  }, [selectionGroups, pinnedObsColumns, pinnedGenes, pinnedObsData, pinnedGeneData,
      colorMode, selectedObsColumn, selectedGene, _categoryCodes, categoryMap, _expressionData, geneLabelMap])

  return results
}

async function computeCategorySummary(
  name: string,
  codes: Uint8Array,
  categoryMap: { label: string; color: RGB }[],
  groups: SelectionGroup[],
  version: number,
  versionRef: { current: number },
): Promise<CategorySummaryResult> {
  const countsByGroup = new Map<number, Uint32Array>()

  await Promise.all(groups.map(async (group) => {
    const response = await getPool().dispatch<CategorySummaryResponse>({
      type: 'summarizeCategory',
      codes,
      indices: group.indices,
      numCategories: categoryMap.length,
      version,
    })
    if (versionRef.current === version) {
      countsByGroup.set(group.id, response.counts)
    }
  }))

  return { type: 'category', name, categoryMap, countsByGroup }
}

async function computeExpressionSummary(
  name: string,
  expression: Float32Array,
  groups: SelectionGroup[],
  version: number,
  versionRef: { current: number },
): Promise<ExpressionSummaryResult> {
  const statsByGroup = new Map<number, {
    mean: number; median: number; std: number; min: number; max: number
    bins: Uint32Array; binEdges: Float32Array
  }>()

  await Promise.all(groups.map(async (group) => {
    const response = await getPool().dispatch<ExpressionSummaryResponse>({
      type: 'summarizeExpression',
      expression,
      indices: group.indices,
      numBins: NUM_BINS,
      version,
    })
    if (versionRef.current === version) {
      statsByGroup.set(group.id, {
        mean: response.mean,
        median: response.median,
        std: response.std,
        min: response.min,
        max: response.max,
        bins: response.bins,
        binEdges: response.binEdges,
      })
    }
  }))

  return { type: 'expression', name, statsByGroup }
}
