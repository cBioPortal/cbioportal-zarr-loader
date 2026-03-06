import { describe, it, expect } from 'vitest'
import { handleSelectionMessage } from './selection.handler'

describe('selection handler', () => {
  // Square polygon: (0,0) to (10,10)
  const square: [number, number][] = [[0, 0], [10, 0], [10, 10], [0, 10]]

  describe('rectangle selection', () => {
    it('finds points inside rectangle bounds', () => {
      // 3 points: (5,5) inside, (15,5) outside, (5,15) outside
      const positions = new Float32Array([5, 5, 15, 5, 5, 15])
      const result = handleSelectionMessage({
        type: 'pointsInPolygon',
        positions,
        numPoints: 3,
        polygon: square,
        selectionType: 'rectangle',
        version: 1,
      })

      expect(result.type).toBe('selectionResult')
      expect(result.indices).toBeInstanceOf(Uint32Array)
      expect(Array.from(result.indices)).toEqual([0])
      expect(result.version).toBe(1)
    })

    it('returns empty indices when no points match', () => {
      const positions = new Float32Array([20, 20, 30, 30])
      const result = handleSelectionMessage({
        type: 'pointsInPolygon',
        positions,
        numPoints: 2,
        polygon: square,
        selectionType: 'rectangle',
        version: 1,
      })

      expect(result.indices.length).toBe(0)
    })

    it('finds all points when all are inside', () => {
      const positions = new Float32Array([1, 1, 5, 5, 9, 9])
      const result = handleSelectionMessage({
        type: 'pointsInPolygon',
        positions,
        numPoints: 3,
        polygon: square,
        selectionType: 'rectangle',
        version: 1,
      })

      expect(Array.from(result.indices)).toEqual([0, 1, 2])
    })
  })

  describe('lasso selection', () => {
    it('finds points inside lasso polygon', () => {
      const triangle: [number, number][] = [[0, 0], [10, 0], [5, 10]]
      // (5,5) inside triangle, (8,8) outside triangle
      const positions = new Float32Array([5, 5, 8, 8])
      const result = handleSelectionMessage({
        type: 'pointsInPolygon',
        positions,
        numPoints: 2,
        polygon: triangle,
        selectionType: 'lasso',
        version: 2,
      })

      expect(Array.from(result.indices)).toEqual([0])
      expect(result.version).toBe(2)
    })

    it('works with concave polygon', () => {
      // L-shape
      const lShape: [number, number][] = [[0, 0], [10, 0], [10, 5], [5, 5], [5, 10], [0, 10]]
      // (2,2) inside, (8,8) outside (top-right of L)
      const positions = new Float32Array([2, 2, 8, 8])
      const result = handleSelectionMessage({
        type: 'pointsInPolygon',
        positions,
        numPoints: 2,
        polygon: lShape,
        selectionType: 'lasso',
        version: 1,
      })

      expect(Array.from(result.indices)).toEqual([0])
    })
  })

  it('echoes version in response', () => {
    const positions = new Float32Array([5, 5])
    const result = handleSelectionMessage({
      type: 'pointsInPolygon',
      positions,
      numPoints: 1,
      polygon: square,
      selectionType: 'rectangle',
      version: 99,
    })
    expect(result.version).toBe(99)
  })
})
