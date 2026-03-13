import { describe, it, expect, beforeEach, vi } from 'vitest'
import { interpolateColorScale, COLOR_SCALES, CATEGORICAL_COLORS } from '../utils/colors'

interface PostedMessage {
  _poolTaskId: number
  type: string
  buffer: Uint8Array
  version: number
}

const postedMessages: PostedMessage[] = []
globalThis.postMessage = (data: PostedMessage) => {
  const cloned = { ...data }
  if (data.buffer instanceof Uint8Array) {
    cloned.buffer = new Uint8Array(data.buffer)
  }
  postedMessages.push(cloned)
}

await import('./universal.worker')

const handler = globalThis.onmessage as (e: { data: unknown }) => void

describe('universal worker', () => {
  beforeEach(() => {
    postedMessages.length = 0
  })

  describe('colorBuffer routing', () => {
    it('handles buildDefault and echoes _poolTaskId', () => {
      handler({ data: { _poolTaskId: 5, type: 'buildDefault', numPoints: 2, rgb: [100, 150, 255], alpha: 0.5, version: 1 } })

      expect(postedMessages).toHaveLength(1)
      expect(postedMessages[0]._poolTaskId).toBe(5)
      expect(postedMessages[0].type).toBe('colorBuffer')
      expect(postedMessages[0].buffer).toBeInstanceOf(Uint8Array)
      expect(postedMessages[0].buffer.length).toBe(2 * 4)
    })

    it('handles buildFromExpression', () => {
      handler({
        data: {
          _poolTaskId: 6,
          type: 'buildFromExpression',
          numPoints: 1,
          expression: new Float32Array([0.5]),
          min: 0,
          max: 1,
          alpha: 1.0,
          scaleName: 'viridis',
          version: 1,
        },
      })

      expect(postedMessages).toHaveLength(1)
      expect(postedMessages[0]._poolTaskId).toBe(6)
      expect(postedMessages[0].type).toBe('colorBuffer')
    })

    it('handles buildFromCategories', () => {
      handler({
        data: {
          _poolTaskId: 7,
          type: 'buildFromCategories',
          numPoints: 1,
          categories: new Uint8Array([0]),
          alpha: 1.0,
          highlightedCodes: null,
          version: 1,
        },
      })

      expect(postedMessages).toHaveLength(1)
      expect(postedMessages[0]._poolTaskId).toBe(7)
      expect(postedMessages[0].type).toBe('colorBuffer')
    })
  })

  describe('validation', () => {
    it('ignores messages with unknown type', () => {
      const spy = vi.spyOn(console, 'warn').mockImplementation(() => {})
      handler({ data: { _poolTaskId: 1, type: 'unknown' } })
      expect(postedMessages).toHaveLength(0)
      expect(spy).toHaveBeenCalled()
      spy.mockRestore()
    })

    it('ignores messages missing required fields', () => {
      const spy = vi.spyOn(console, 'warn').mockImplementation(() => {})
      handler({ data: { _poolTaskId: 1, type: 'buildDefault' } })
      expect(postedMessages).toHaveLength(0)
      expect(spy).toHaveBeenCalled()
      spy.mockRestore()
    })
  })
})
