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
        highlightedCodes: null,
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
        highlightedCodes: null,
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
        highlightedCodes: null,
        version: 99,
      })
      expect(result.version).toBe(99)
    })

    it('returns null radiusBuffer when no highlights', () => {
      const result = handleColorBufferMessage({
        type: 'buildFromCategories',
        numPoints: 2,
        categories: new Uint8Array([0, 1]),
        alpha: 0.5,
        highlightedCodes: null,
        version: 1,
      })
      expect(result.radiusBuffer).toBeNull()
    })

    it('grays out non-highlighted points and produces radiusBuffer', () => {
      const result = handleColorBufferMessage({
        type: 'buildFromCategories',
        numPoints: 3,
        categories: new Uint8Array([0, 1, 2]),
        alpha: 0.5,
        highlightedCodes: [1],
        version: 1,
      })

      const { buffer, radiusBuffer } = result
      expect(radiusBuffer).toBeInstanceOf(Float32Array)
      expect(radiusBuffer!.length).toBe(3)

      // Point 0: code 0, NOT highlighted → gray, radius 0.5
      expect(buffer[0]).toBe(200)
      expect(buffer[1]).toBe(200)
      expect(buffer[2]).toBe(200)
      expect(buffer[3]).toBe(Math.round(0.5 * 255))
      expect(radiusBuffer![0]).toBe(0.5)

      // Point 1: code 1, highlighted → category color, full alpha, radius 1.0
      const color1 = CATEGORICAL_COLORS[1]
      expect(buffer[4]).toBe(color1[0])
      expect(buffer[5]).toBe(color1[1])
      expect(buffer[6]).toBe(color1[2])
      expect(buffer[7]).toBe(255)
      expect(radiusBuffer![1]).toBe(1.0)

      // Point 2: code 2, NOT highlighted → gray, radius 0.5
      expect(buffer[8]).toBe(200)
      expect(buffer[9]).toBe(200)
      expect(buffer[10]).toBe(200)
      expect(radiusBuffer![2]).toBe(0.5)
    })

    it('highlights multiple codes simultaneously', () => {
      const result = handleColorBufferMessage({
        type: 'buildFromCategories',
        numPoints: 3,
        categories: new Uint8Array([0, 1, 2]),
        alpha: 0.5,
        highlightedCodes: [0, 2],
        version: 1,
      })

      const { buffer, radiusBuffer } = result

      // Point 0: highlighted
      const color0 = CATEGORICAL_COLORS[0]
      expect(buffer[0]).toBe(color0[0])
      expect(radiusBuffer![0]).toBe(1.0)

      // Point 1: NOT highlighted → gray
      expect(buffer[4]).toBe(200)
      expect(radiusBuffer![1]).toBe(0.5)

      // Point 2: highlighted
      const color2 = CATEGORICAL_COLORS[2]
      expect(buffer[8]).toBe(color2[0])
      expect(radiusBuffer![2]).toBe(1.0)
    })

    it('treats empty highlightedCodes array same as null', () => {
      const result = handleColorBufferMessage({
        type: 'buildFromCategories',
        numPoints: 2,
        categories: new Uint8Array([0, 1]),
        alpha: 0.8,
        highlightedCodes: [],
        version: 1,
      })

      // No highlights active — normal category colors, no radius buffer
      expect(result.radiusBuffer).toBeNull()
      const color0 = CATEGORICAL_COLORS[0]
      expect(result.buffer[0]).toBe(color0[0])
      expect(result.buffer[3]).toBe(Math.round(0.8 * 255))
    })
  })
})
