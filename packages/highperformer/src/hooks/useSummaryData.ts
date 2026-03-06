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
  const summaryObsColumns = useAppStore((s) => s.summaryObsColumns)
  const summaryGenes = useAppStore((s) => s.summaryGenes)
  const summaryObsData = useAppStore((s) => s.summaryObsData)
  const summaryObsContinuousData = useAppStore((s) => s.summaryObsContinuousData)
  const summaryGeneData = useAppStore((s) => s.summaryGeneData)
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

    // Pinned obs columns
    for (const name of summaryObsColumns) {
      const catData = summaryObsData.get(name)
      if (catData) {
        tasks.push(computeCategorySummary(name, catData.codes, catData.categoryMap, groupsWithIndices, version, versionRef))
        continue
      }
      const contData = summaryObsContinuousData.get(name)
      if (contData) {
        tasks.push(computeExpressionSummary(name, contData, groupsWithIndices, version, versionRef))
      }
    }

    // Pinned genes
    for (const name of summaryGenes) {
      const data = summaryGeneData.get(name)
      if (!data) continue
      const displayName = geneLabelMap?.get(name) ?? name
      tasks.push(computeExpressionSummary(displayName, data, groupsWithIndices, version, versionRef))
    }

    Promise.all(tasks).then((resolved) => {
      if (versionRef.current !== version) return
      setResults(resolved)
    })
  }, [selectionGroups, summaryObsColumns, summaryGenes, summaryObsData, summaryObsContinuousData, summaryGeneData, geneLabelMap])

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
