import { useEffect, useMemo, useRef, useState } from 'react'
import useAppStore, { getPool } from '../store/useAppStore'
import type { CategorySummaryResponse, ExpressionSummaryResponse } from '../workers/summary.schemas'
import type { SummaryResult, CategorySummaryResult, ExpressionSummaryResult } from './useSummaryData'
import type { RGB } from '../utils/colors'

const ALL_CELLS_GROUP_ID = -1
const NUM_BINS = 30

export function useAllCellsSummary(): SummaryResult[] {
  const embeddingData = useAppStore((s) => s.embeddingData)
  const summaryObsColumns = useAppStore((s) => s.summaryObsColumns)
  const summaryGenes = useAppStore((s) => s.summaryGenes)
  const summaryObsData = useAppStore((s) => s.summaryObsData)
  const summaryObsContinuousData = useAppStore((s) => s.summaryObsContinuousData)
  const summaryGeneData = useAppStore((s) => s.summaryGeneData)
  const geneLabelMap = useAppStore((s) => s.geneLabelMap)

  const numPoints = embeddingData?.numPoints ?? 0

  // Build sequential indices array once, reuse across computations
  const allIndices = useMemo(() => {
    if (numPoints === 0) return new Uint32Array(0)
    const arr = new Uint32Array(numPoints)
    for (let i = 0; i < numPoints; i++) arr[i] = i
    return arr
  }, [numPoints])

  const [results, setResults] = useState<SummaryResult[]>([])
  const versionRef = useRef(0)

  useEffect(() => {
    if (numPoints === 0) {
      setResults([])
      return
    }

    const hasVariables = summaryObsColumns.length > 0 || summaryGenes.length > 0
    if (!hasVariables) {
      setResults([])
      return
    }

    versionRef.current++
    const version = versionRef.current

    const tasks: Promise<SummaryResult>[] = []

    for (const name of summaryObsColumns) {
      const catData = summaryObsData.get(name)
      if (catData) {
        tasks.push(computeAllCellsCategory(name, catData.codes, catData.categoryMap, allIndices, version, versionRef))
        continue
      }
      const contData = summaryObsContinuousData.get(name)
      if (contData) {
        tasks.push(computeAllCellsExpression(name, name, contData, allIndices, version, versionRef))
      }
    }

    for (const name of summaryGenes) {
      const data = summaryGeneData.get(name)
      if (!data) continue
      const displayName = geneLabelMap?.get(name) ?? name
      tasks.push(computeAllCellsExpression(displayName, name, data, allIndices, version, versionRef))
    }

    Promise.all(tasks).then((resolved) => {
      if (versionRef.current !== version) return
      setResults(resolved)
    })
  }, [numPoints, allIndices, summaryObsColumns, summaryGenes, summaryObsData, summaryObsContinuousData, summaryGeneData, geneLabelMap])

  return results
}

async function computeAllCellsCategory(
  name: string,
  codes: Uint8Array,
  categoryMap: { label: string; color: RGB }[],
  indices: Uint32Array,
  version: number,
  versionRef: { current: number },
): Promise<CategorySummaryResult> {
  const countsByGroup = new Map<number, Uint32Array>()

  const response = await getPool().dispatch<CategorySummaryResponse>({
    type: 'summarizeCategory',
    codes,
    indices,
    numCategories: categoryMap.length,
    version,
  })
  if (versionRef.current === version) {
    countsByGroup.set(ALL_CELLS_GROUP_ID, response.counts)
  }

  return { type: 'category', name, categoryMap, countsByGroup }
}

async function computeAllCellsExpression(
  name: string,
  dataKey: string,
  expression: Float32Array,
  indices: Uint32Array,
  version: number,
  versionRef: { current: number },
): Promise<ExpressionSummaryResult> {
  const statsByGroup = new Map<number, {
    mean: number; median: number; std: number; min: number; max: number
    q1: number; q3: number; whiskerLow: number; whiskerHigh: number
    bins: Uint32Array; binEdges: Float32Array
    kdeX: Float32Array; kdeDensity: Float32Array; clippedCount: number
  }>()

  const response = await getPool().dispatch<ExpressionSummaryResponse>({
    type: 'summarizeExpression',
    expression,
    indices,
    numBins: NUM_BINS,
    version,
  })
  if (versionRef.current === version) {
    statsByGroup.set(ALL_CELLS_GROUP_ID, {
      mean: response.mean,
      median: response.median,
      std: response.std,
      min: response.min,
      max: response.max,
      q1: response.q1,
      q3: response.q3,
      whiskerLow: response.whiskerLow,
      whiskerHigh: response.whiskerHigh,
      bins: response.bins,
      binEdges: response.binEdges,
      kdeX: response.kdeX,
      kdeDensity: response.kdeDensity,
      clippedCount: response.clippedCount ?? 0,
    })
  }

  return { type: 'expression', name, dataKey, statsByGroup }
}

export { ALL_CELLS_GROUP_ID }
