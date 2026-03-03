/**
 * Ensure a typed array is Float32Array.
 * Returns the input unchanged if already Float32; otherwise creates a new Float32Array copy.
 *
 * @param {TypedArray} data
 * @returns {Float32Array}
 */
export function ensureFloat32(data) {
  if (data instanceof Float32Array) return data;
  return new Float32Array(data);
}
