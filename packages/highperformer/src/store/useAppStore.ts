import { create } from 'zustand'
import { AnnDataStore, GENE_SYMBOL_COLUMNS } from '@cbioportal-zarr-loader/zarrstore'
import type { ArrayResult } from '@cbioportal-zarr-loader/zarrstore'
import { WorkerPool } from '../pool/WorkerPool'
import UniversalWorker from '../workers/universal.worker.ts?worker'
import type { ColorBufferResponse } from '../workers/colorBuffer.schemas'
import type { RGB } from '../utils/colors'
import { encodeCategories, MAX_CATEGORIES } from '../utils/categoryEncoding'

export interface EmbeddingBounds {
  minX: number
  maxX: number
  minY: number
  maxY: number
}

export interface EmbeddingData {
  positions: Float32Array
  numPoints: number
  bounds: EmbeddingBounds
}

export type ColorMode = 'category' | 'gene'

export interface AppState {
  // Dataset
  datasetUrl: string | null
  adata: AnnDataStore | null
  loading: boolean

  // Metadata derived from adata (cached on load)
  nObs: number | null
  nVar: number | null
  obsmKeys: string[]

  // Embedding selection
  selectedEmbedding: string | null

  // Rendering controls
  pointRadius: number
  opacity: number
  antialiasing: boolean
  collisionEnabled: boolean
  collisionRadiusScale: number
  setPointRadius: (v: number) => void
  setOpacity: (v: number) => void
  setAntialiasing: (v: boolean) => void
  setCollisionEnabled: (v: boolean) => void
  setCollisionRadiusScale: (v: number) => void

  // Embedding data — typed arrays for direct GPU upload
  embeddingData: EmbeddingData | null
  embeddingLoading: boolean
  _embeddingAbort: AbortController | null

  // Color buffer — Uint8Array(numPoints * 4), RGBA per point
  colorBuffer: Uint8Array | null
  colorBufferLoading: boolean

  // Color By state
  colorMode: ColorMode
  selectedObsColumn: string | null
  selectedGene: string | null
  colorScaleName: string
  obsColumnNames: string[]
  varNames: string[]
  categoryMap: { label: string; color: RGB }[]
  expressionRange: { min: number; max: number } | null
  categoryWarning: string | null

  // Gene label resolution
  varColumns: string[]
  geneLabelColumn: string | null
  geneLabelMap: Map<string, string> | null

  // Internal — cached data for rebuilds (not for UI consumption)
  _categoryCodes: Uint8Array | null
  _expressionData: Float32Array | null
  _colorAbort: AbortController | null

  // Actions
  openDataset: (url: string) => Promise<void>
  setSelectedEmbedding: (key: string) => void
  fetchEmbedding: (key: string) => Promise<void>
  rebuildColorBuffer: () => void
  setColorMode: (mode: ColorMode) => void
  selectObsColumn: (name: string) => void
  clearObsColumn: () => void
  selectGene: (name: string) => void
  clearGene: () => void
  setColorScaleName: (name: string) => void
  setGeneLabelColumn: (col: string | null) => void
  _resolveGeneLabels: () => Promise<void>
}

const DEFAULT_RGB: RGB = [100, 150, 255]

function computeBounds(positions: Float32Array, numPoints: number): EmbeddingBounds {
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity
  for (let i = 0; i < numPoints; i++) {
    const x = positions[i * 2]
    const y = positions[i * 2 + 1]
    if (x < minX) minX = x
    if (x > maxX) maxX = x
    if (y < minY) minY = y
    if (y > maxY) maxY = y
  }
  return { minX, maxX, minY, maxY }
}

function computeRange(data: Float32Array): { min: number; max: number } {
  let min = Infinity, max = -Infinity
  for (let i = 0; i < data.length; i++) {
    if (data[i] < min) min = data[i]
    if (data[i] > max) max = data[i]
  }
  return { min, max }
}

// Version counter for stale response detection
let colorBuildVersion = 0
export function getColorBuildVersion(): number { return colorBuildVersion }
export function resetColorBuildVersion(): void { colorBuildVersion = 0 }

// Singleton pool — created lazily
let pool: WorkerPool | null = null
function getPool(): WorkerPool {
  if (!pool) pool = new WorkerPool(() => new UniversalWorker())
  return pool
}

// Debounce timer for opacity-driven color buffer rebuilds
let debounceTimer: ReturnType<typeof setTimeout> | null = null
const DEBOUNCE_MS = 150

const useAppStore = create<AppState>((set, get) => ({
  // Dataset
  datasetUrl: null,
  adata: null,
  loading: false,

  // Metadata derived from adata (cached on load)
  nObs: null,
  nVar: null,
  obsmKeys: [],

  // Embedding selection
  selectedEmbedding: null,

  // Rendering controls
  pointRadius: 1,
  opacity: 0.3,
  antialiasing: true,
  collisionEnabled: false,
  collisionRadiusScale: 0,
  setPointRadius: (v) => set({ pointRadius: v }),
  setOpacity: (v) => {
    set({ opacity: v, colorBufferLoading: true })
    if (debounceTimer) clearTimeout(debounceTimer)
    debounceTimer = setTimeout(() => get().rebuildColorBuffer(), DEBOUNCE_MS)
  },
  setAntialiasing: (v) => set({ antialiasing: v }),
  setCollisionEnabled: (v) => set({ collisionEnabled: v }),
  setCollisionRadiusScale: (v) => set({ collisionRadiusScale: v }),

  // Embedding data — typed arrays for direct GPU upload
  embeddingData: null,
  embeddingLoading: false,
  _embeddingAbort: null,

  // Color buffer — Uint8Array(numPoints * 4), RGBA per point
  colorBuffer: null,
  colorBufferLoading: false,

  // Color By state
  colorMode: 'category',
  selectedObsColumn: null,
  selectedGene: null,
  colorScaleName: 'viridis',
  obsColumnNames: [],
  varNames: [],
  categoryMap: [],
  expressionRange: null,
  categoryWarning: null,

  // Gene label resolution
  varColumns: [],
  geneLabelColumn: null,
  geneLabelMap: null,

  // Internal
  _categoryCodes: null,
  _expressionData: null,
  _colorAbort: null,

  // Actions
  openDataset: async (url) => {
    if (url === get().datasetUrl && get().adata) return
    set({
      datasetUrl: url, loading: true, adata: null, nObs: null, nVar: null, obsmKeys: [],
      selectedEmbedding: null, embeddingData: null, colorBuffer: null,
      colorMode: 'category', selectedObsColumn: null, selectedGene: null,
      obsColumnNames: [], varNames: [], categoryMap: [], expressionRange: null,
      categoryWarning: null, _categoryCodes: null, _expressionData: null,
      varColumns: [], geneLabelColumn: null, geneLabelMap: null,
    })
    try {
      const adata = await AnnDataStore.open(url)
      const obsmKeys = adata.obsmKeys()
      const umap = obsmKeys.find((k) => /umap/i.test(k))
      const defaultKey = umap ?? obsmKeys[0] ?? null
      set({
        adata,
        nObs: adata.nObs,
        nVar: adata.nVar,
        obsmKeys,
        selectedEmbedding: defaultKey,
        loading: false,
      })
      if (defaultKey) get().fetchEmbedding(defaultKey)
    } catch {
      set({ loading: false })
    }
  },

  setSelectedEmbedding: (key) => {
    set({ selectedEmbedding: key })
    get().fetchEmbedding(key)
  },

  fetchEmbedding: async (key) => {
    const { adata, _embeddingAbort } = get()
    if (!adata || !key) return

    // Abort any in-flight fetch
    if (_embeddingAbort) _embeddingAbort.abort()
    const abortController = new AbortController()
    set({ embeddingLoading: true, _embeddingAbort: abortController })

    try {
      const result = await adata.obsm(key, abortController.signal, 2) as ArrayResult
      // Ensure Float32Array — deck.gl doesn't support Float16Array as a vertex attribute
      const positions = result.data instanceof Float32Array
        ? result.data
        : new Float32Array(result.data as ArrayLike<number>)
      const numPoints = result.shape[0]
      const bounds = computeBounds(positions, numPoints)
      set({
        embeddingData: { positions, numPoints, bounds },
        embeddingLoading: false,
        _embeddingAbort: null,
      })
      // Build default color buffer for the new embedding
      get().rebuildColorBuffer()

      // Fetch obs column names, var names, and var columns in the background
      Promise.all([adata.obsColumns(), adata.varNames(), adata.varColumns()]).then(([obsColumnNames, varNamesRaw, varCols]) => {
        const varNames = varNamesRaw.map((v) => String(v ?? ''))

        // Auto-detect gene label column
        const colsLower = varCols.map((c) => c.toLowerCase())
        let detectedCol: string | null = null
        for (const candidate of GENE_SYMBOL_COLUMNS) {
          // Exact match first
          const exactIdx = varCols.indexOf(candidate)
          if (exactIdx !== -1) { detectedCol = varCols[exactIdx]; break }
          // Case-insensitive fallback
          const lowerIdx = colsLower.indexOf(candidate.toLowerCase())
          if (lowerIdx !== -1) { detectedCol = varCols[lowerIdx]; break }
        }

        set({ obsColumnNames, varNames, varColumns: varCols, geneLabelColumn: detectedCol })
        if (detectedCol) get()._resolveGeneLabels()
      })
    } catch (err) {
      if (err instanceof Error && err.name !== 'AbortError') {
        set({ embeddingLoading: false, _embeddingAbort: null })
      }
    }
  },

  rebuildColorBuffer: () => {
    const { embeddingData, opacity, colorMode, _categoryCodes, _expressionData, expressionRange, colorScaleName } = get()
    if (!embeddingData) return

    colorBuildVersion++
    const version = colorBuildVersion

    let message: Record<string, unknown>

    if (colorMode === 'category' && _categoryCodes) {
      message = {
        type: 'buildFromCategories',
        numPoints: embeddingData.numPoints,
        categories: _categoryCodes,
        alpha: opacity,
        version,
      }
    } else if (colorMode === 'gene' && _expressionData && expressionRange) {
      message = {
        type: 'buildFromExpression',
        numPoints: embeddingData.numPoints,
        expression: _expressionData,
        min: expressionRange.min,
        max: expressionRange.max,
        alpha: opacity,
        scaleName: colorScaleName,
        version,
      }
    } else {
      message = {
        type: 'buildDefault',
        numPoints: embeddingData.numPoints,
        rgb: DEFAULT_RGB,
        alpha: opacity,
        version,
      }
    }

    getPool()
      .dispatch<ColorBufferResponse>(message)
      .then((response) => {
        if (version !== colorBuildVersion) return // stale
        set({ colorBuffer: response.buffer, colorBufferLoading: false })
      })
  },

  setColorMode: (mode) => {
    set({
      colorMode: mode,
      categoryWarning: null,
    })
    // When switching modes, rebuild if cached data exists for the new mode
    if (mode === 'category' && get()._categoryCodes) {
      get().rebuildColorBuffer()
    } else if (mode === 'gene' && get()._expressionData) {
      get().rebuildColorBuffer()
    }
  },

  selectObsColumn: (name) => {
    const { adata, _colorAbort } = get()
    if (!adata) return

    if (_colorAbort) _colorAbort.abort()
    const abortController = new AbortController()
    set({
      selectedObsColumn: name,
      colorBufferLoading: true,
      categoryWarning: null,
      _colorAbort: abortController,
    })

    adata.obsColumn(name, abortController.signal).then((values) => {
      const valuesArray = Array.isArray(values) ? values : Array.from(values as Iterable<number>)
      const { codes, categoryMap, uniqueCount } = encodeCategories(valuesArray as (string | number | null)[])

      if (uniqueCount > MAX_CATEGORIES) {
        set({
          categoryWarning: `This column has ${uniqueCount} unique values (likely continuous). Please choose a categorical column.`,
          colorBufferLoading: false,
          _categoryCodes: null,
          categoryMap: [],
          _colorAbort: null,
        })
        return
      }

      set({
        _categoryCodes: codes,
        categoryMap,
        categoryWarning: null,
        _colorAbort: null,
      })
      get().rebuildColorBuffer()
    })
  },

  clearObsColumn: () => {
    const { _colorAbort } = get()
    if (_colorAbort) _colorAbort.abort()
    set({
      selectedObsColumn: null,
      _categoryCodes: null,
      categoryMap: [],
      categoryWarning: null,
      _colorAbort: null,
    })
    get().rebuildColorBuffer()
  },

  selectGene: (name) => {
    const { adata, _colorAbort } = get()
    if (!adata) return

    if (_colorAbort) _colorAbort.abort()
    const abortController = new AbortController()
    set({
      selectedGene: name,
      colorBufferLoading: true,
      _colorAbort: abortController,
    })

    adata.geneExpression(name, abortController.signal).then((expression) => {
      const data = expression instanceof Float32Array
        ? expression
        : new Float32Array(expression as ArrayLike<number>)
      const range = computeRange(data)
      set({
        _expressionData: data,
        expressionRange: range,
        _colorAbort: null,
      })
      get().rebuildColorBuffer()
    })
  },

  clearGene: () => {
    const { _colorAbort } = get()
    if (_colorAbort) _colorAbort.abort()
    set({
      selectedGene: null,
      _expressionData: null,
      expressionRange: null,
      _colorAbort: null,
    })
    get().rebuildColorBuffer()
  },

  setColorScaleName: (name) => {
    set({ colorScaleName: name })
    if (get().colorMode === 'gene' && get()._expressionData) {
      get().rebuildColorBuffer()
    }
  },

  setGeneLabelColumn: (col) => {
    set({ geneLabelColumn: col, geneLabelMap: null })
    if (col) {
      get()._resolveGeneLabels()
    }
  },

  _resolveGeneLabels: async () => {
    const { adata, geneLabelColumn, varNames } = get()
    if (!adata || !geneLabelColumn) return

    try {
      const [symbols, varNamesRaw] = await Promise.all([
        adata.varColumn(geneLabelColumn),
        varNames.length > 0 ? Promise.resolve(varNames) : adata.varNames().then((raw) => raw.map((v) => String(v ?? ''))),
      ])
      const symbolsArr = Array.isArray(symbols) ? symbols : Array.from(symbols as Iterable<unknown>)
      const map = new Map<string, string>()
      for (let i = 0; i < varNamesRaw.length; i++) {
        const varIndex = String(varNamesRaw[i])
        const label = String(symbolsArr[i] ?? '')
        if (label && label !== varIndex) {
          map.set(varIndex, label)
        }
      }
      // Only update if the column hasn't changed while we were fetching
      if (get().geneLabelColumn === geneLabelColumn) {
        set({ geneLabelMap: map })
      }
    } catch {
      // Silently fail — labels are a nicety, not critical
    }
  },
}))

export default useAppStore
