import { describe, it, expect, beforeEach } from 'vitest'
import { interpolateColorScale, COLOR_SCALES, CATEGORICAL_COLORS } from '../utils/colors'

interface PostedMessage {
  type: string
  buffer: Uint8Array
}

// Capture posted messages before the buffer gets transferred/detached.
// The worker calls self.postMessage(data, [buf.buffer]) which would detach the typed array.
const postedMessages: PostedMessage[] = []
globalThis.postMessage = (data: PostedMessage) => {
  // Clone the buffer before it could be detached
  const cloned = { ...data }
  if (data.buffer instanceof Uint8Array) {
    cloned.buffer = new Uint8Array(data.buffer)
  }
  postedMessages.push(cloned)
}

// Dynamic import so our postMessage stub is in place first
await import('./colorBuffer.worker')

// The worker sets self.onmessage — grab a reference to it
const handler = globalThis.onmessage as (e: { data: unknown }) => void

describe('colorBuffer worker', () => {
  beforeEach(() => {
    postedMessages.length = 0
  })

  describe('buildDefault', () => {
    it('produces uniform RGBA buffer', () => {
      handler({ data: { type: 'buildDefault', numPoints: 3, rgb: [100, 150, 255], alpha: 0.5 } })

      expect(postedMessages).toHaveLength(1)
      const { type, buffer } = postedMessages[0]
      expect(type).toBe('colorBuffer')
      expect(buffer).toBeInstanceOf(Uint8Array)
      expect(buffer.length).toBe(3 * 4)

      const expectedAlpha = Math.round(0.5 * 255)
      for (let i = 0; i < 3; i++) {
        const off = i * 4
        expect(buffer[off]).toBe(100)
        expect(buffer[off + 1]).toBe(150)
        expect(buffer[off + 2]).toBe(255)
        expect(buffer[off + 3]).toBe(expectedAlpha)
      }
    })

    it('computes alpha correctly from float', () => {
      handler({ data: { type: 'buildDefault', numPoints: 1, rgb: [0, 0, 0], alpha: 1.0 } })
      expect(postedMessages[0].buffer[3]).toBe(255)

      postedMessages.length = 0
      handler({ data: { type: 'buildDefault', numPoints: 1, rgb: [0, 0, 0], alpha: 0.0 } })
      expect(postedMessages[0].buffer[3]).toBe(0)
    })
  })

  describe('buildFromExpression', () => {
    it('maps expression values through color scale', () => {
      const expression = new Float32Array([0, 5, 10])
      handler({
        data: {
          type: 'buildFromExpression',
          numPoints: 3,
          expression,
          min: 0,
          max: 10,
          alpha: 1.0,
          scaleName: 'viridis',
        },
      })

      const { buffer } = postedMessages[0]
      expect(buffer).toBeInstanceOf(Uint8Array)
      expect(buffer.length).toBe(3 * 4)

      // First point (t=0) should match viridis start
      const [r0, g0, b0] = interpolateColorScale(0, COLOR_SCALES.viridis)
      expect(buffer[0]).toBe(r0)
      expect(buffer[1]).toBe(g0)
      expect(buffer[2]).toBe(b0)

      // Last point (t=1) should match viridis end
      const [r2, g2, b2] = interpolateColorScale(1, COLOR_SCALES.viridis)
      expect(buffer[8]).toBe(r2)
      expect(buffer[9]).toBe(g2)
      expect(buffer[10]).toBe(b2)
    })

    it('defaults to viridis for unknown scale name', () => {
      const expression = new Float32Array([0])
      handler({
        data: {
          type: 'buildFromExpression',
          numPoints: 1,
          expression,
          min: 0,
          max: 1,
          alpha: 1.0,
          scaleName: 'nonexistent',
        },
      })

      const { buffer } = postedMessages[0]
      const [r, g, b] = interpolateColorScale(0, COLOR_SCALES.viridis)
      expect(buffer[0]).toBe(r)
      expect(buffer[1]).toBe(g)
      expect(buffer[2]).toBe(b)
    })
  })

  describe('buildFromCategories', () => {
    it('maps category indices to categorical colors', () => {
      const categories = new Uint8Array([0, 1, 2])
      handler({
        data: { type: 'buildFromCategories', numPoints: 3, categories, alpha: 0.8 },
      })

      const { buffer } = postedMessages[0]
      const expectedAlpha = Math.round(0.8 * 255)

      for (let i = 0; i < 3; i++) {
        const off = i * 4
        const color = CATEGORICAL_COLORS[i]
        expect(buffer[off]).toBe(color[0])
        expect(buffer[off + 1]).toBe(color[1])
        expect(buffer[off + 2]).toBe(color[2])
        expect(buffer[off + 3]).toBe(expectedAlpha)
      }
    })

    it('wraps category indices beyond palette length', () => {
      const categories = new Uint8Array([15]) // 15 % 15 = 0
      handler({
        data: { type: 'buildFromCategories', numPoints: 1, categories, alpha: 1.0 },
      })

      const { buffer } = postedMessages[0]
      const color = CATEGORICAL_COLORS[0]
      expect(buffer[0]).toBe(color[0])
      expect(buffer[1]).toBe(color[1])
      expect(buffer[2]).toBe(color[2])
    })
  })
})
