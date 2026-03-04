import { describe, it, expect } from 'vitest'
import { interpolateColorScale, COLOR_SCALES, CATEGORICAL_COLORS } from './colors'

describe('interpolateColorScale', () => {
  const scale = COLOR_SCALES.viridis

  it('returns the first color at t=0', () => {
    expect(interpolateColorScale(0, scale)).toEqual([68, 1, 84])
  })

  it('returns the last color at t=1', () => {
    expect(interpolateColorScale(1, scale)).toEqual([253, 231, 37])
  })

  it('interpolates between colors at t=0.5', () => {
    const [r, g, b] = interpolateColorScale(0.5, scale)
    // t=0.5 maps to index 4.5 in a 10-stop scale, so it's between stops 4 and 5
    expect(r).toBeGreaterThan(0)
    expect(g).toBeGreaterThan(0)
    expect(b).toBeGreaterThan(0)
    // Should be between the two middle colors
    expect(r).toBeLessThan(255)
    expect(g).toBeLessThan(255)
    expect(b).toBeLessThan(255)
  })

  it('clamps t values below 0', () => {
    expect(interpolateColorScale(-0.5, scale)).toEqual(interpolateColorScale(0, scale))
  })

  it('clamps t values above 1', () => {
    expect(interpolateColorScale(1.5, scale)).toEqual(interpolateColorScale(1, scale))
  })

  it('returns integer RGB values', () => {
    const [r, g, b] = interpolateColorScale(0.33, scale)
    expect(Number.isInteger(r)).toBe(true)
    expect(Number.isInteger(g)).toBe(true)
    expect(Number.isInteger(b)).toBe(true)
  })
})

describe('COLOR_SCALES', () => {
  it('has viridis and magma entries', () => {
    expect(COLOR_SCALES).toHaveProperty('viridis')
    expect(COLOR_SCALES).toHaveProperty('magma')
  })

  it('each scale has 10 color stops', () => {
    expect(COLOR_SCALES.viridis).toHaveLength(10)
    expect(COLOR_SCALES.magma).toHaveLength(10)
  })

  it('each stop is an RGB triplet', () => {
    for (const stop of COLOR_SCALES.viridis) {
      expect(stop).toHaveLength(3)
      stop.forEach((v) => {
        expect(v).toBeGreaterThanOrEqual(0)
        expect(v).toBeLessThanOrEqual(255)
      })
    }
  })
})

describe('CATEGORICAL_COLORS', () => {
  it('has 15 colors', () => {
    expect(CATEGORICAL_COLORS).toHaveLength(15)
  })

  it('each entry is an RGB triplet', () => {
    for (const color of CATEGORICAL_COLORS) {
      expect(color).toHaveLength(3)
      color.forEach((v) => {
        expect(v).toBeGreaterThanOrEqual(0)
        expect(v).toBeLessThanOrEqual(255)
      })
    }
  })
})
