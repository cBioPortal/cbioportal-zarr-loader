import { describe, it, expect, vi } from "vitest";
import { renderHook } from "@testing-library/react";
import useSelectionInteraction from "./useSelectionInteraction";

describe("useSelectionInteraction", () => {
  const baseMocks = {
    deckRef: { current: null },
    setSelectedPoints: vi.fn(),
    setSelectionGeometry: vi.fn(),
    clearSelectedPoints: vi.fn(),
  };

  it("accepts points array (backward compat)", () => {
    const points = [
      { position: [1, 2], index: 0 },
      { position: [3, 4], index: 1 },
    ];
    const { result } = renderHook(() =>
      useSelectionInteraction({ ...baseMocks, points }),
    );
    expect(result.current.selectMode).toBe("pan");
  });

  it("accepts positionBuffer and stride", () => {
    const positionBuffer = new Float32Array([1, 2, 3, 4]);
    const { result } = renderHook(() =>
      useSelectionInteraction({
        ...baseMocks,
        positionBuffer,
        stride: 2,
        numPoints: 2,
      }),
    );
    expect(result.current.selectMode).toBe("pan");
  });
});
