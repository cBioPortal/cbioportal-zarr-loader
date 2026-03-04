import { interpolateColorScale, COLOR_SCALES, CATEGORICAL_COLORS } from '../utils/colors'
import type { RGB } from '../utils/colors'

// In a Web Worker context, workerSelf.postMessage takes (message, transfer[])
const workerSelf = self as unknown as {
  onmessage: ((e: MessageEvent) => void) | null
  postMessage(message: unknown, transfer: Transferable[]): void
}

interface BuildDefaultMsg {
  type: 'buildDefault'
  numPoints: number
  rgb: RGB
  alpha: number
}

interface BuildFromExpressionMsg {
  type: 'buildFromExpression'
  numPoints: number
  expression: Float32Array
  min: number
  max: number
  alpha: number
  scaleName: string
}

interface BuildFromCategoriesMsg {
  type: 'buildFromCategories'
  numPoints: number
  categories: Uint8Array
  alpha: number
}

type WorkerMessage = BuildDefaultMsg | BuildFromExpressionMsg | BuildFromCategoriesMsg

/**
 * Web Worker for building color buffers off the main thread.
 *
 * Message types:
 *   buildDefault  — uniform RGBA for all points
 *   buildFromExpression — continuous color scale from Float32Array values
 *   buildFromCategories — categorical colors from integer category indices
 */
workerSelf.onmessage = (e: MessageEvent<WorkerMessage>) => {
  const { type } = e.data

  if (type === 'buildDefault') {
    const { numPoints, rgb, alpha } = e.data
    const buf = new Uint8Array(numPoints * 4)
    const a = Math.round(alpha * 255)
    for (let i = 0; i < numPoints; i++) {
      const off = i * 4
      buf[off] = rgb[0]
      buf[off + 1] = rgb[1]
      buf[off + 2] = rgb[2]
      buf[off + 3] = a
    }
    workerSelf.postMessage({ type: 'colorBuffer', buffer: buf }, [buf.buffer] as Transferable[])
    return
  }

  if (type === 'buildFromExpression') {
    const { numPoints, expression, min, max, alpha, scaleName } = e.data
    const scale = COLOR_SCALES[scaleName] || COLOR_SCALES.viridis
    const buf = new Uint8Array(numPoints * 4)
    const a = Math.round(alpha * 255)
    const range = max - min || 1
    for (let i = 0; i < numPoints; i++) {
      const t = (expression[i] - min) / range
      const [r, g, b] = interpolateColorScale(t, scale)
      const off = i * 4
      buf[off] = r
      buf[off + 1] = g
      buf[off + 2] = b
      buf[off + 3] = a
    }
    workerSelf.postMessage({ type: 'colorBuffer', buffer: buf }, [buf.buffer] as Transferable[])
    return
  }

  if (type === 'buildFromCategories') {
    const { numPoints, categories, alpha } = e.data
    const buf = new Uint8Array(numPoints * 4)
    const a = Math.round(alpha * 255)
    const numColors = CATEGORICAL_COLORS.length
    for (let i = 0; i < numPoints; i++) {
      const color = CATEGORICAL_COLORS[categories[i] % numColors]
      const off = i * 4
      buf[off] = color[0]
      buf[off + 1] = color[1]
      buf[off + 2] = color[2]
      buf[off + 3] = a
    }
    workerSelf.postMessage({ type: 'colorBuffer', buffer: buf }, [buf.buffer] as Transferable[])
    return
  }
}
