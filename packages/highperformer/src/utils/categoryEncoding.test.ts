import { describe, it, expect } from 'vitest'
import { encodeCategories, MAX_CATEGORIES } from './categoryEncoding'
import { CATEGORICAL_COLORS } from './colors'

describe('encodeCategories', () => {
  it('encodes string values to sequential integer codes', () => {
    const values = ['A', 'B', 'A', 'C', 'B']
    const result = encodeCategories(values)

    expect(result.codes).toBeInstanceOf(Uint8Array)
    expect(result.codes.length).toBe(5)
    // Same strings get same codes
    expect(result.codes[0]).toBe(result.codes[2]) // A === A
    expect(result.codes[1]).toBe(result.codes[4]) // B === B
    // Different strings get different codes
    expect(result.codes[0]).not.toBe(result.codes[1])
  })

  it('builds a categoryMap with labels and colors', () => {
    const values = ['cat', 'dog', 'cat']
    const result = encodeCategories(values)

    expect(result.categoryMap).toHaveLength(2)
    expect(result.categoryMap[0].label).toBe('cat')
    expect(result.categoryMap[0].color).toEqual(CATEGORICAL_COLORS[0])
    expect(result.categoryMap[1].label).toBe('dog')
    expect(result.categoryMap[1].color).toEqual(CATEGORICAL_COLORS[1])
  })

  it('returns uniqueCount', () => {
    const values = ['A', 'B', 'C', 'A']
    const result = encodeCategories(values)
    expect(result.uniqueCount).toBe(3)
  })

  it('handles numeric values by converting to string', () => {
    const values = [1, 2, 1, 3] as unknown as (string | number | null)[]
    const result = encodeCategories(values)
    expect(result.uniqueCount).toBe(3)
    expect(result.codes.length).toBe(4)
  })

  it('handles null values with a "null" category', () => {
    const values = ['A', null, 'A'] as (string | number | null)[]
    const result = encodeCategories(values)
    expect(result.uniqueCount).toBe(2)
    expect(result.codes.length).toBe(3)
  })

  it('wraps colors via modulo for > 15 categories', () => {
    const values = Array.from({ length: 20 }, (_, i) => `cat_${i}`)
    const result = encodeCategories(values)
    expect(result.uniqueCount).toBe(20)
    // Color at index 15 wraps to CATEGORICAL_COLORS[0]
    expect(result.categoryMap[15].color).toEqual(CATEGORICAL_COLORS[0])
  })

  it('exports MAX_CATEGORIES as 1000', () => {
    expect(MAX_CATEGORIES).toBe(1000)
  })
})
