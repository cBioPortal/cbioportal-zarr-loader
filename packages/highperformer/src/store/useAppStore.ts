import { create } from 'zustand'
import { AnnDataStore } from '@cbioportal-zarr-loader/zarrstore'
import type { ArrayResult } from '@cbioportal-zarr-loader/zarrstore'
import ColorBufferWorker from '../workers/colorBuffer.worker.ts?worker'
import type { RGB } from '../utils/colors'

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

export interface AppState {
  // Dataset
  datasetUrl: string | null
  adata: AnnDataStore | null
  loading: boolean

  // Metadata derived from adata (cached on load)
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

  // Actions
  openDataset: (url: string) => Promise<void>
  setSelectedEmbedding: (key: string) => void
  fetchEmbedding: (key: string) => Promise<void>
  rebuildColorBuffer: () => void
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

// Singleton worker — created lazily, shared across all store actions
let colorWorker: Worker | null = null
function getColorWorker(): Worker {
  if (!colorWorker) colorWorker = new ColorBufferWorker()
  return colorWorker
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

  // Actions
  openDataset: async (url) => {
    if (url === get().datasetUrl && get().adata) return
    set({ datasetUrl: url, loading: true, adata: null, obsmKeys: [], selectedEmbedding: null, embeddingData: null, colorBuffer: null })
    try {
      const adata = await AnnDataStore.open(url)
      const obsmKeys = adata.obsmKeys()
      const umap = obsmKeys.find((k) => /umap/i.test(k))
      const defaultKey = umap ?? obsmKeys[0] ?? null
      set({
        adata,
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
    } catch (err) {
      if (err instanceof Error && err.name !== 'AbortError') {
        set({ embeddingLoading: false, _embeddingAbort: null })
      }
    }
  },

  rebuildColorBuffer: () => {
    const { embeddingData, opacity } = get()
    if (!embeddingData) return

    const worker = getColorWorker()

    // One-shot listener for this build — replaces any previous listener
    worker.onmessage = (e) => {
      if (e.data.type === 'colorBuffer') {
        set({ colorBuffer: e.data.buffer, colorBufferLoading: false })
      }
    }

    worker.postMessage({
      type: 'buildDefault',
      numPoints: embeddingData.numPoints,
      rgb: DEFAULT_RGB,
      alpha: opacity,
    })
  },
}))

export default useAppStore
