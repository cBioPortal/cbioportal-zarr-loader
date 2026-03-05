import { describe, it, expect, beforeEach, vi, afterEach } from 'vitest'

// Mock the WorkerPool before importing the store
const mockDispatch = vi.fn()
vi.mock('../pool/WorkerPool', () => {
  return {
    WorkerPool: class MockPool {
      dispatch = mockDispatch
      dispose() {}
    },
  }
})

// Mock the worker import
vi.mock('../workers/universal.worker.ts?worker', () => {
  return { default: class MockWorker {} }
})

const { default: useAppStore, getColorBuildVersion, resetColorBuildVersion } = await import('./useAppStore')

describe('useAppStore', () => {
  beforeEach(() => {
    vi.useFakeTimers()
    useAppStore.setState(useAppStore.getInitialState())
    mockDispatch.mockClear()
    resetColorBuildVersion()
    mockDispatch.mockResolvedValue({ type: 'colorBuffer', buffer: new Uint8Array(8), version: 1 })
  })

  afterEach(() => {
    vi.useRealTimers()
  })

  describe('initial state', () => {
    it('has correct default values', () => {
      const state = useAppStore.getState()
      expect(state.pointRadius).toBe(1)
      expect(state.opacity).toBe(0.3)
      expect(state.antialiasing).toBe(true)
      expect(state.collisionEnabled).toBe(false)
      expect(state.collisionRadiusScale).toBe(0)
      expect(state.datasetUrl).toBeNull()
      expect(state.adata).toBeNull()
      expect(state.loading).toBe(false)
      expect(state.embeddingData).toBeNull()
      expect(state.colorBuffer).toBeNull()
      expect(state.colorBufferLoading).toBe(false)
      expect(state.obsmKeys).toEqual([])
      expect(state.selectedEmbedding).toBeNull()
    })

    it('has correct color state defaults', () => {
      const state = useAppStore.getState()
      expect(state.colorMode).toBe('category')
      expect(state.selectedObsColumn).toBeNull()
      expect(state.selectedGene).toBeNull()
      expect(state.colorScaleName).toBe('viridis')
      expect(state.obsColumnNames).toEqual([])
      expect(state.varNames).toEqual([])
      expect(state.categoryMap).toEqual([])
      expect(state.expressionRange).toBeNull()
      expect(state.categoryWarning).toBeNull()
    })
  })

  describe('rendering control setters', () => {
    it('setPointRadius updates pointRadius', () => {
      useAppStore.getState().setPointRadius(3)
      expect(useAppStore.getState().pointRadius).toBe(3)
    })

    it('setAntialiasing updates antialiasing', () => {
      useAppStore.getState().setAntialiasing(false)
      expect(useAppStore.getState().antialiasing).toBe(false)
    })

    it('setCollisionEnabled updates collisionEnabled', () => {
      useAppStore.getState().setCollisionEnabled(true)
      expect(useAppStore.getState().collisionEnabled).toBe(true)
    })

    it('setCollisionRadiusScale updates collisionRadiusScale', () => {
      useAppStore.getState().setCollisionRadiusScale(5)
      expect(useAppStore.getState().collisionRadiusScale).toBe(5)
    })
  })

  describe('setOpacity', () => {
    it('updates opacity immediately', () => {
      useAppStore.getState().setOpacity(0.8)
      expect(useAppStore.getState().opacity).toBe(0.8)
    })

    it('sets colorBufferLoading to true immediately', () => {
      useAppStore.getState().setOpacity(0.5)
      expect(useAppStore.getState().colorBufferLoading).toBe(true)
    })

    it('debounces rebuildColorBuffer — does not fire immediately', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
      })
      mockDispatch.mockClear()

      useAppStore.getState().setOpacity(0.5)
      expect(mockDispatch).not.toHaveBeenCalled()
    })

    it('fires rebuildColorBuffer after debounce delay', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
      })
      mockDispatch.mockClear()
      mockDispatch.mockResolvedValue({ type: 'colorBuffer', buffer: new Uint8Array(8), version: 1 })

      useAppStore.getState().setOpacity(0.5)
      vi.advanceTimersByTime(150)

      expect(mockDispatch).toHaveBeenCalledWith({
        type: 'buildDefault',
        numPoints: 2,
        rgb: [100, 150, 255],
        alpha: 0.5,
        version: expect.any(Number),
      })
    })

    it('coalesces rapid calls — only last opacity value fires', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
      })
      mockDispatch.mockClear()
      mockDispatch.mockResolvedValue({ type: 'colorBuffer', buffer: new Uint8Array(8), version: 1 })

      useAppStore.getState().setOpacity(0.3)
      useAppStore.getState().setOpacity(0.5)
      useAppStore.getState().setOpacity(0.8)
      vi.advanceTimersByTime(150)

      expect(mockDispatch).toHaveBeenCalledTimes(1)
      expect(mockDispatch).toHaveBeenCalledWith({
        type: 'buildDefault',
        numPoints: 2,
        rgb: [100, 150, 255],
        alpha: 0.8,
        version: expect.any(Number),
      })
    })
  })

  describe('rebuildColorBuffer', () => {
    it('does nothing when embeddingData is null', () => {
      useAppStore.getState().rebuildColorBuffer()
      expect(mockDispatch).not.toHaveBeenCalled()
    })

    it('dispatches buildDefault when no cached data for active mode', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1, 2, 2]),
          numPoints: 3,
          bounds: { minX: 0, maxX: 2, minY: 0, maxY: 2 },
        },
        opacity: 0.7,
        colorMode: 'category',
        _categoryCodes: null,
      })

      useAppStore.getState().rebuildColorBuffer()

      expect(mockDispatch).toHaveBeenCalledWith({
        type: 'buildDefault',
        numPoints: 3,
        rgb: [100, 150, 255],
        alpha: 0.7,
        version: 1,
      })
    })

    it('dispatches buildFromCategories in category mode', () => {
      const codes = new Uint8Array([0, 1, 0])
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1, 2, 2]),
          numPoints: 3,
          bounds: { minX: 0, maxX: 2, minY: 0, maxY: 2 },
        },
        opacity: 0.7,
        colorMode: 'category',
        _categoryCodes: codes,
      })

      useAppStore.getState().rebuildColorBuffer()

      expect(mockDispatch).toHaveBeenCalledWith({
        type: 'buildFromCategories',
        numPoints: 3,
        categories: codes,
        alpha: 0.7,
        version: 1,
      })
    })

    it('dispatches buildFromExpression in gene mode', () => {
      const expression = new Float32Array([0.5, 1.0, 0.0])
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1, 2, 2]),
          numPoints: 3,
          bounds: { minX: 0, maxX: 2, minY: 0, maxY: 2 },
        },
        opacity: 0.7,
        colorMode: 'gene',
        colorScaleName: 'viridis',
        _expressionData: expression,
        expressionRange: { min: 0.0, max: 1.0 },
      })

      useAppStore.getState().rebuildColorBuffer()

      expect(mockDispatch).toHaveBeenCalledWith({
        type: 'buildFromExpression',
        numPoints: 3,
        expression,
        min: 0.0,
        max: 1.0,
        alpha: 0.7,
        scaleName: 'viridis',
        version: 1,
      })
    })

    it('matching version response updates store', async () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorBufferLoading: true,
      })

      const fakeBuffer = new Uint8Array(8)
      mockDispatch.mockResolvedValue({ type: 'colorBuffer', buffer: fakeBuffer, version: 1 })

      useAppStore.getState().rebuildColorBuffer()

      await vi.waitFor(() => {
        expect(useAppStore.getState().colorBufferLoading).toBe(false)
      })
      expect(useAppStore.getState().colorBuffer).toBe(fakeBuffer)
    })

    it('discards stale responses with outdated version', async () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorBufferLoading: true,
      })

      let resolveFirst!: (v: unknown) => void
      mockDispatch.mockImplementationOnce(() => new Promise((r) => { resolveFirst = r }))

      useAppStore.getState().rebuildColorBuffer()

      const currentBuffer = new Uint8Array([10, 20, 30, 40, 50, 60, 70, 80])
      mockDispatch.mockResolvedValueOnce({ type: 'colorBuffer', buffer: currentBuffer, version: 2 })

      useAppStore.getState().rebuildColorBuffer()

      await vi.waitFor(() => {
        expect(useAppStore.getState().colorBuffer).toBe(currentBuffer)
      })

      const staleBuffer = new Uint8Array([1, 2, 3, 4, 5, 6, 7, 8])
      resolveFirst({ type: 'colorBuffer', buffer: staleBuffer, version: 1 })

      await vi.waitFor(() => {
        expect(useAppStore.getState().colorBuffer).toBe(currentBuffer)
      })
    })
  })

  describe('setColorMode', () => {
    it('switches mode and rebuilds if cached data exists', () => {
      const codes = new Uint8Array([0, 1])
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorMode: 'gene',
        _categoryCodes: codes,
      })
      mockDispatch.mockClear()

      useAppStore.getState().setColorMode('category')

      expect(useAppStore.getState().colorMode).toBe('category')
      expect(useAppStore.getState().categoryWarning).toBeNull()
      expect(mockDispatch).toHaveBeenCalledWith(
        expect.objectContaining({ type: 'buildFromCategories' }),
      )
    })

    it('does not rebuild when switching to mode without cached data', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorMode: 'category',
        _expressionData: null,
      })
      mockDispatch.mockClear()

      useAppStore.getState().setColorMode('gene')

      expect(useAppStore.getState().colorMode).toBe('gene')
      expect(mockDispatch).not.toHaveBeenCalled()
    })
  })

  describe('selectObsColumn', () => {
    it('fetches column, encodes categories, and rebuilds buffer', async () => {
      const mockAdata = {
        obsColumn: vi.fn().mockResolvedValue(['A', 'B', 'A']),
      }
      useAppStore.setState({
        adata: mockAdata as any,
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1, 2, 2]),
          numPoints: 3,
          bounds: { minX: 0, maxX: 2, minY: 0, maxY: 2 },
        },
        colorMode: 'category',
      })
      mockDispatch.mockClear()

      useAppStore.getState().selectObsColumn('cluster')

      await vi.waitFor(() => {
        expect(useAppStore.getState().selectedObsColumn).toBe('cluster')
      })
      expect(mockAdata.obsColumn).toHaveBeenCalledWith('cluster', expect.any(AbortSignal))
      expect(useAppStore.getState().categoryMap.length).toBe(2)
      expect(useAppStore.getState().categoryWarning).toBeNull()
      expect(mockDispatch).toHaveBeenCalledWith(
        expect.objectContaining({ type: 'buildFromCategories' }),
      )
    })

    it('shows warning when unique values exceed threshold', async () => {
      const manyValues = Array.from({ length: 1001 }, (_, i) => `val_${i}`)
      const mockAdata = {
        obsColumn: vi.fn().mockResolvedValue(manyValues),
      }
      useAppStore.setState({
        adata: mockAdata as any,
        embeddingData: {
          positions: new Float32Array(1001 * 2),
          numPoints: 1001,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorMode: 'category',
      })
      mockDispatch.mockClear()

      useAppStore.getState().selectObsColumn('continuous_col')

      await vi.waitFor(() => {
        expect(useAppStore.getState().categoryWarning).toContain('1001')
      })
      // Should NOT dispatch a color buffer build
      expect(mockDispatch).not.toHaveBeenCalled()
    })
  })

  describe('clearObsColumn', () => {
    it('clears selection, cached data, and rebuilds with default', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorMode: 'category',
        selectedObsColumn: 'cluster',
        _categoryCodes: new Uint8Array([0, 1]),
        categoryMap: [{ label: 'A', color: [0, 0, 0] }],
        categoryWarning: 'some warning',
      })
      mockDispatch.mockClear()

      useAppStore.getState().clearObsColumn()

      expect(useAppStore.getState().selectedObsColumn).toBeNull()
      expect(useAppStore.getState()._categoryCodes).toBeNull()
      expect(useAppStore.getState().categoryMap).toEqual([])
      expect(useAppStore.getState().categoryWarning).toBeNull()
      expect(mockDispatch).toHaveBeenCalledWith(
        expect.objectContaining({ type: 'buildDefault' }),
      )
    })
  })

  describe('selectGene', () => {
    it('fetches expression data, computes range, and rebuilds buffer', async () => {
      const expression = new Float32Array([0.0, 2.5, 5.0])
      const mockAdata = {
        geneExpression: vi.fn().mockResolvedValue(expression),
      }
      useAppStore.setState({
        adata: mockAdata as any,
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1, 2, 2]),
          numPoints: 3,
          bounds: { minX: 0, maxX: 2, minY: 0, maxY: 2 },
        },
        colorMode: 'gene',
        colorScaleName: 'viridis',
      })
      mockDispatch.mockClear()

      useAppStore.getState().selectGene('TP53')

      await vi.waitFor(() => {
        expect(useAppStore.getState().selectedGene).toBe('TP53')
      })
      expect(mockAdata.geneExpression).toHaveBeenCalledWith('TP53', expect.any(AbortSignal))
      expect(useAppStore.getState().expressionRange).toEqual({ min: 0.0, max: 5.0 })
      expect(mockDispatch).toHaveBeenCalledWith(
        expect.objectContaining({
          type: 'buildFromExpression',
          expression,
          min: 0.0,
          max: 5.0,
          scaleName: 'viridis',
        }),
      )
    })
  })

  describe('clearGene', () => {
    it('clears selection, cached data, and rebuilds with default', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorMode: 'gene',
        selectedGene: 'TP53',
        _expressionData: new Float32Array([0, 1]),
        expressionRange: { min: 0, max: 1 },
      })
      mockDispatch.mockClear()

      useAppStore.getState().clearGene()

      expect(useAppStore.getState().selectedGene).toBeNull()
      expect(useAppStore.getState()._expressionData).toBeNull()
      expect(useAppStore.getState().expressionRange).toBeNull()
      expect(mockDispatch).toHaveBeenCalledWith(
        expect.objectContaining({ type: 'buildDefault' }),
      )
    })
  })

  describe('setColorScaleName', () => {
    it('updates scale and rebuilds when in gene mode with expression data', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorMode: 'gene',
        _expressionData: new Float32Array([0, 1]),
        expressionRange: { min: 0, max: 1 },
      })
      mockDispatch.mockClear()

      useAppStore.getState().setColorScaleName('magma')

      expect(useAppStore.getState().colorScaleName).toBe('magma')
      expect(mockDispatch).toHaveBeenCalledWith(
        expect.objectContaining({ type: 'buildFromExpression', scaleName: 'magma' }),
      )
    })
  })

  describe('gene label resolution', () => {
    it('auto-detects gene_symbol column and populates geneLabelMap', async () => {
      const mockAdata = {
        obsmKeys: () => ['X_umap'],
        obsm: vi.fn().mockResolvedValue({ data: new Float32Array([0, 0, 1, 1]), shape: [2, 2] }),
        obsColumns: vi.fn().mockResolvedValue(['n_cells', 'cluster']),
        varNames: vi.fn().mockResolvedValue(['ENSG001', 'ENSG002']),
        varColumns: vi.fn().mockResolvedValue(['n_cells', 'gene_symbol']),
        varColumn: vi.fn().mockResolvedValue(['TP53', 'BRCA1']),
      }
      useAppStore.setState({ adata: mockAdata as any })

      await useAppStore.getState().fetchEmbedding('X_umap')

      await vi.waitFor(() => {
        expect(useAppStore.getState().geneLabelColumn).toBe('gene_symbol')
      })
      await vi.waitFor(() => {
        expect(useAppStore.getState().geneLabelMap).not.toBeNull()
      })
      const map = useAppStore.getState().geneLabelMap!
      expect(map.get('ENSG001')).toBe('TP53')
      expect(map.get('ENSG002')).toBe('BRCA1')
    })

    it('sets geneLabelColumn to null when no symbol column exists', async () => {
      const mockAdata = {
        obsmKeys: () => ['X_umap'],
        obsm: vi.fn().mockResolvedValue({ data: new Float32Array([0, 0, 1, 1]), shape: [2, 2] }),
        obsColumns: vi.fn().mockResolvedValue(['n_cells']),
        varNames: vi.fn().mockResolvedValue(['ENSG001', 'ENSG002']),
        varColumns: vi.fn().mockResolvedValue(['n_cells']),
      }
      useAppStore.setState({ adata: mockAdata as any })

      await useAppStore.getState().fetchEmbedding('X_umap')

      await vi.waitFor(() => {
        expect(useAppStore.getState().varColumns).toEqual(['n_cells'])
      })
      expect(useAppStore.getState().geneLabelColumn).toBeNull()
      expect(useAppStore.getState().geneLabelMap).toBeNull()
    })

    it('user override via setGeneLabelColumn fetches new labels', async () => {
      const mockAdata = {
        varColumn: vi.fn().mockResolvedValue(['GeneA', 'GeneB']),
        varNames: vi.fn().mockResolvedValue(['IDX1', 'IDX2']),
      }
      useAppStore.setState({
        adata: mockAdata as any,
        varNames: ['IDX1', 'IDX2'],
        varColumns: ['feature_name', 'alt_name'],
        geneLabelColumn: null,
        geneLabelMap: null,
      })

      useAppStore.getState().setGeneLabelColumn('alt_name')

      await vi.waitFor(() => {
        expect(useAppStore.getState().geneLabelMap).not.toBeNull()
      })
      expect(useAppStore.getState().geneLabelColumn).toBe('alt_name')
      expect(mockAdata.varColumn).toHaveBeenCalledWith('alt_name')
      expect(useAppStore.getState().geneLabelMap!.get('IDX1')).toBe('GeneA')
    })

    it('setGeneLabelColumn(null) clears the map', () => {
      useAppStore.setState({
        geneLabelColumn: 'gene_symbol',
        geneLabelMap: new Map([['ENSG001', 'TP53']]),
      })

      useAppStore.getState().setGeneLabelColumn(null)

      expect(useAppStore.getState().geneLabelColumn).toBeNull()
      expect(useAppStore.getState().geneLabelMap).toBeNull()
    })

    it('selectGene still passes var index to geneExpression', async () => {
      const expression = new Float32Array([0.5, 1.0])
      const mockAdata = {
        geneExpression: vi.fn().mockResolvedValue(expression),
      }
      useAppStore.setState({
        adata: mockAdata as any,
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorMode: 'gene' as const,
        geneLabelMap: new Map([['ENSG001', 'TP53']]),
      })
      mockDispatch.mockClear()

      useAppStore.getState().selectGene('ENSG001')

      await vi.waitFor(() => {
        expect(mockAdata.geneExpression).toHaveBeenCalledWith('ENSG001', expect.any(AbortSignal))
      })
    })
  })
})
