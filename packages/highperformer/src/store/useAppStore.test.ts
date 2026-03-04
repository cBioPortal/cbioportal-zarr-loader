import { describe, it, expect, beforeEach, vi, afterEach } from 'vitest'

// Mock the worker import before importing the store
const mockPostMessage = vi.fn()
let mockWorkerInstance: { onmessage: ((e: MessageEvent) => void) | null } | null = null
vi.mock('../workers/colorBuffer.worker.ts?worker', () => {
  return {
    default: class MockWorker {
      onmessage: ((e: MessageEvent) => void) | null = null
      constructor() {
        mockWorkerInstance = this
      }
      postMessage(...args: unknown[]) {
        mockPostMessage(...args)
      }
      terminate() {}
    },
  }
})

const { default: useAppStore } = await import('./useAppStore')

describe('useAppStore', () => {
  beforeEach(() => {
    vi.useFakeTimers()
    useAppStore.setState(useAppStore.getInitialState())
    mockPostMessage.mockClear()
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
      mockPostMessage.mockClear()

      useAppStore.getState().setOpacity(0.5)

      // Worker should NOT have been called yet (debounce pending)
      expect(mockPostMessage).not.toHaveBeenCalled()
    })

    it('fires rebuildColorBuffer after debounce delay', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
      })
      mockPostMessage.mockClear()

      useAppStore.getState().setOpacity(0.5)
      vi.advanceTimersByTime(150)

      expect(mockPostMessage).toHaveBeenCalledWith({
        type: 'buildDefault',
        numPoints: 2,
        rgb: [100, 150, 255],
        alpha: 0.5,
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
      mockPostMessage.mockClear()

      useAppStore.getState().setOpacity(0.3)
      useAppStore.getState().setOpacity(0.5)
      useAppStore.getState().setOpacity(0.8)
      vi.advanceTimersByTime(150)

      expect(mockPostMessage).toHaveBeenCalledTimes(1)
      expect(mockPostMessage).toHaveBeenCalledWith({
        type: 'buildDefault',
        numPoints: 2,
        rgb: [100, 150, 255],
        alpha: 0.8,
      })
    })
  })

  describe('rebuildColorBuffer', () => {
    it('does nothing when embeddingData is null', () => {
      useAppStore.getState().rebuildColorBuffer()
      expect(mockPostMessage).not.toHaveBeenCalled()
    })

    it('posts buildDefault message to worker', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1, 2, 2]),
          numPoints: 3,
          bounds: { minX: 0, maxX: 2, minY: 0, maxY: 2 },
        },
        opacity: 0.7,
      })

      useAppStore.getState().rebuildColorBuffer()

      expect(mockPostMessage).toHaveBeenCalledWith({
        type: 'buildDefault',
        numPoints: 3,
        rgb: [100, 150, 255],
        alpha: 0.7,
      })
    })

    it('worker response sets colorBufferLoading to false', () => {
      useAppStore.setState({
        embeddingData: {
          positions: new Float32Array([0, 0, 1, 1]),
          numPoints: 2,
          bounds: { minX: 0, maxX: 1, minY: 0, maxY: 1 },
        },
        colorBufferLoading: true,
      })

      // Trigger rebuildColorBuffer to set up worker.onmessage
      useAppStore.getState().rebuildColorBuffer()

      // Simulate the worker posting back a colorBuffer result
      const fakeBuffer = new Uint8Array(8)
      mockWorkerInstance!.onmessage!({ data: { type: 'colorBuffer', buffer: fakeBuffer } } as MessageEvent)

      expect(useAppStore.getState().colorBufferLoading).toBe(false)
      expect(useAppStore.getState().colorBuffer).toBe(fakeBuffer)
    })
  })
})
