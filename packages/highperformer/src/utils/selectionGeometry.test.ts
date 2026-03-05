import { describe, it, expect } from 'vitest'
import { pointInPolygon, simplifyPolygon } from './selectionGeometry'

describe('pointInPolygon', () => {
  const square: [number, number][] = [[0, 0], [10, 0], [10, 10], [0, 10]]

  it('returns true for point inside polygon', () => {
    expect(pointInPolygon(5, 5, square)).toBe(true)
  })

  it('returns false for point outside polygon', () => {
    expect(pointInPolygon(15, 5, square)).toBe(false)
  })

  it('returns false for point above polygon', () => {
    expect(pointInPolygon(5, 15, square)).toBe(false)
  })

  it('works with triangle', () => {
    const triangle: [number, number][] = [[0, 0], [10, 0], [5, 10]]
    expect(pointInPolygon(5, 5, triangle)).toBe(true)
    expect(pointInPolygon(1, 9, triangle)).toBe(false)
  })

  it('works with concave polygon', () => {
    // L-shape
    const lShape: [number, number][] = [[0, 0], [10, 0], [10, 5], [5, 5], [5, 10], [0, 10]]
    expect(pointInPolygon(2, 2, lShape)).toBe(true)  // bottom-left
    expect(pointInPolygon(2, 8, lShape)).toBe(true)  // top-left
    expect(pointInPolygon(8, 2, lShape)).toBe(true)  // bottom-right
    expect(pointInPolygon(8, 8, lShape)).toBe(false) // top-right (outside L)
  })
})

describe('simplifyPolygon', () => {
  it('returns polygon unchanged if 3 or fewer points', () => {
    const tri: [number, number][] = [[0, 0], [1, 0], [0, 1]]
    expect(simplifyPolygon(tri)).toEqual(tri)
  })

  it('removes collinear points', () => {
    const line: [number, number][] = [[0, 0], [1, 0], [2, 0], [3, 0], [4, 0]]
    const result = simplifyPolygon(line)
    expect(result.length).toBeLessThan(line.length)
    expect(result[0]).toEqual([0, 0])
    expect(result[result.length - 1]).toEqual([4, 0])
  })

  it('preserves corners', () => {
    // Square with many intermediate points on each edge
    const dense: [number, number][] = [
      [0, 0], [2, 0], [4, 0], [6, 0], [8, 0], [10, 0],
      [10, 2], [10, 4], [10, 6], [10, 8], [10, 10],
      [8, 10], [6, 10], [4, 10], [2, 10], [0, 10],
      [0, 8], [0, 6], [0, 4], [0, 2],
    ]
    const result = simplifyPolygon(dense)
    expect(result.length).toBeLessThanOrEqual(5) // 4 corners + maybe start/end
  })
})
