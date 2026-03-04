import { describe, it, expect } from 'vitest'
import { interpolateColorScale, COLOR_SCALES, CATEGORICAL_COLORS } from '../utils/colors'
import { handleColorBufferMessage } from './colorBuffer.handler'

describe('colorBuffer handler', () => {
  describe('buildDefault', () => {
    it('produces uniform RGBA buffer', () => {
      const result = handleColorBufferMessage({
        type: 'buildDefault',
        numPoints: 3,
        rgb: [100, 150, 255],
        alpha: 0.5,
        version: 1,
      })

      expect(result.type).toBe('colorBuffer')
      expect(result.buffer).toBeInstanceOf(Uint8Array)
      expect(result.buffer.length).toBe(3 * 4)
      expect(result.version).toBe(1)

      const expectedAlpha = Math.round(0.5 * 255)
      for (let i = 0; i < 3; i++) {
        const off = i * 4
        expect(result.buffer[off]).toBe(100)
        expect(result.buffer[off + 1]).toBe(150)
        expect(result.buffer[off + 2]).toBe(255)
        expect(result.buffer[off + 3]).toBe(expectedAlpha)
      }
    })

    it('computes alpha correctly from float', () => {
      const r1 = handleColorBufferMessage({ type: 'buildDefault', numPoints: 1, rgb: [0, 0, 0], alpha: 1.0, version: 1 })
      expect(r1.buffer[3]).toBe(255)

      const r2 = handleColorBufferMessage({ type: 'buildDefault', numPoints: 1, rgb: [0, 0, 0], alpha: 0.0, version: 1 })
      expect(r2.buffer[3]).toBe(0)
    })

    it('echoes version in response', () => {
      const result = handleColorBufferMessage({ type: 'buildDefault', numPoints: 1, rgb: [0, 0, 0], alpha: 1.0, version: 42 })
      expect(result.version).toBe(42)
    })
  })

  describe('buildFromExpression', () => {
    it('maps expression values through color scale', () => {
      const expression = new Float32Array([0, 5, 10])
      const result = handleColorBufferMessage({
        type: 'buildFromExpression',
        numPoints: 3,
        expression,
        min: 0,
        max: 10,
        alpha: 1.0,
        scaleName: 'viridis',
        version: 1,
      })

      const { buffer } = result
      expect(buffer).toBeInstanceOf(Uint8Array)
      expect(buffer.length).toBe(3 * 4)

      const [r0, g0, b0] = interpolateColorScale(0, COLOR_SCALES.viridis)
      expect(buffer[0]).toBe(r0)
      expect(buffer[1]).toBe(g0)
      expect(buffer[2]).toBe(b0)

      const [r2, g2, b2] = interpolateColorScale(1, COLOR_SCALES.viridis)
      expect(buffer[8]).toBe(r2)
      expect(buffer[9]).toBe(g2)
      expect(buffer[10]).toBe(b2)
    })

    it('defaults to viridis for unknown scale name', () => {
      const result = handleColorBufferMessage({
        type: 'buildFromExpression',
        numPoints: 1,
        expression: new Float32Array([0]),
        min: 0,
        max: 1,
        alpha: 1.0,
        scaleName: 'nonexistent',
        version: 1,
      })

      const [r, g, b] = interpolateColorScale(0, COLOR_SCALES.viridis)
      expect(result.buffer[0]).toBe(r)
      expect(result.buffer[1]).toBe(g)
      expect(result.buffer[2]).toBe(b)
    })

    it('echoes version in response', () => {
      const result = handleColorBufferMessage({
        type: 'buildFromExpression',
        numPoints: 1,
        expression: new Float32Array([0]),
        min: 0,
        max: 1,
        alpha: 1.0,
        scaleName: 'viridis',
        version: 7,
      })
      expect(result.version).toBe(7)
    })
  })

  describe('buildFromCategories', () => {
    it('maps category indices to categorical colors', () => {
      const categories = new Uint8Array([0, 1, 2])
      const result = handleColorBufferMessage({
        type: 'buildFromCategories',
        numPoints: 3,
        categories,
        alpha: 0.8,
        version: 1,
      })

      const { buffer } = result
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
      const result = handleColorBufferMessage({
        type: 'buildFromCategories',
        numPoints: 1,
        categories: new Uint8Array([15]),
        alpha: 1.0,
        version: 1,
      })

      const color = CATEGORICAL_COLORS[0]
      expect(result.buffer[0]).toBe(color[0])
      expect(result.buffer[1]).toBe(color[1])
      expect(result.buffer[2]).toBe(color[2])
    })

    it('echoes version in response', () => {
      const result = handleColorBufferMessage({
        type: 'buildFromCategories',
        numPoints: 1,
        categories: new Uint8Array([0]),
        alpha: 1.0,
        version: 99,
      })
      expect(result.version).toBe(99)
    })
  })
})
