import { describe, it, expect } from "vitest";
import { ensureFloat32 } from "./ensureFloat32";

describe("ensureFloat32", () => {
  it("returns Float32Array unchanged", () => {
    const f32 = new Float32Array([1, 2, 3]);
    expect(ensureFloat32(f32)).toBe(f32); // same reference
  });

  it("converts Float64Array to Float32Array", () => {
    const f64 = new Float64Array([1.5, 2.5, 3.5]);
    const result = ensureFloat32(f64);
    expect(result).toBeInstanceOf(Float32Array);
    expect(Array.from(result)).toEqual([
      Math.fround(1.5),
      Math.fround(2.5),
      Math.fround(3.5),
    ]);
  });

  it("converts Int32Array to Float32Array", () => {
    const i32 = new Int32Array([1, 2, 3]);
    const result = ensureFloat32(i32);
    expect(result).toBeInstanceOf(Float32Array);
    expect(Array.from(result)).toEqual([1, 2, 3]);
  });
});
