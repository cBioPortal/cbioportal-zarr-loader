import { describe, it, expect } from "vitest";
import {
  pointInPolygon,
  computeRange,
  computeViewState,
  buildScatterplotPoints,
  buildSelectionSummary,
} from "./scatterplotUtils";
import { CATEGORICAL_COLORS } from "./colors";

// ---------------------------------------------------------------------------
// pointInPolygon
// ---------------------------------------------------------------------------

// Unit square: (0,0), (4,0), (4,4), (0,4)
const square = [[0, 0], [4, 0], [4, 4], [0, 4]];

// Triangle: (0,0), (6,0), (3,6)
const triangle = [[0, 0], [6, 0], [3, 6]];

// L-shaped concave polygon
const lShape = [[0, 0], [4, 0], [4, 2], [2, 2], [2, 4], [0, 4]];

describe("pointInPolygon", () => {
  it("detects point inside a square", () => {
    expect(pointInPolygon(2, 2, square)).toBe(true);
  });

  it("detects point outside a square", () => {
    expect(pointInPolygon(5, 5, square)).toBe(false);
  });

  it("detects point inside a triangle", () => {
    expect(pointInPolygon(3, 2, triangle)).toBe(true);
  });

  it("detects point outside a triangle", () => {
    expect(pointInPolygon(0, 4, triangle)).toBe(false);
  });

  it("detects point inside concave polygon", () => {
    expect(pointInPolygon(1, 3, lShape)).toBe(true);
  });

  it("detects point in concavity as outside", () => {
    // (3, 3) is in the notch of the L — should be outside
    expect(pointInPolygon(3, 3, lShape)).toBe(false);
  });

  it("handles point far outside", () => {
    expect(pointInPolygon(100, 100, square)).toBe(false);
    expect(pointInPolygon(-1, -1, square)).toBe(false);
  });

  it("handles negative coordinates", () => {
    const centered = [[-2, -2], [2, -2], [2, 2], [-2, 2]];
    expect(pointInPolygon(0, 0, centered)).toBe(true);
    expect(pointInPolygon(-1, -1, centered)).toBe(true);
    expect(pointInPolygon(3, 3, centered)).toBe(false);
  });

  it("handles a degenerate polygon (fewer than 3 vertices)", () => {
    // A line segment isn't a polygon — point should be outside
    expect(pointInPolygon(1, 0, [[0, 0], [2, 0]])).toBe(false);
  });

  it("handles an empty polygon", () => {
    expect(pointInPolygon(0, 0, [])).toBe(false);
  });
});

// ---------------------------------------------------------------------------
// computeRange
// ---------------------------------------------------------------------------

describe("computeRange", () => {
  it("returns null for null input", () => {
    expect(computeRange(null)).toBeNull();
  });

  it("returns null for undefined input", () => {
    expect(computeRange(undefined)).toBeNull();
  });

  it("returns null for empty array", () => {
    expect(computeRange([])).toBeNull();
  });

  it("returns correct min and max for a regular array", () => {
    expect(computeRange([3, 1, 4, 1, 5, 9, 2, 6])).toEqual({ min: 1, max: 9 });
  });

  it("handles a single element", () => {
    expect(computeRange([42])).toEqual({ min: 42, max: 42 });
  });

  it("handles negative values", () => {
    expect(computeRange([-5, -1, -10, -3])).toEqual({ min: -10, max: -1 });
  });

  it("handles mixed positive and negative", () => {
    expect(computeRange([-2, 0, 3])).toEqual({ min: -2, max: 3 });
  });

  it("works with Float32Array", () => {
    const arr = new Float32Array([0.1, 0.9, 0.5]);
    const result = computeRange(arr);
    expect(result.min).toBeCloseTo(0.1);
    expect(result.max).toBeCloseTo(0.9);
  });

  it("handles all identical values", () => {
    expect(computeRange([7, 7, 7])).toEqual({ min: 7, max: 7 });
  });
});

// ---------------------------------------------------------------------------
// computeViewState
// ---------------------------------------------------------------------------

const squareBounds = { minX: 0, maxX: 100, minY: 0, maxY: 100 };
const wideBounds = { minX: -50, maxX: 50, minY: -10, maxY: 10 };
const squareContainer = { width: 600, height: 600 };

describe("computeViewState", () => {
  it("returns default state when bounds is null", () => {
    const result = computeViewState(null, squareContainer);
    expect(result).toEqual({ target: [0, 0], zoom: 0 });
  });

  it("centers target on the midpoint of bounds", () => {
    const { target } = computeViewState(squareBounds, squareContainer);
    expect(target).toEqual([50, 50]);
  });

  it("centers target with negative bounds", () => {
    const { target } = computeViewState(wideBounds, squareContainer);
    expect(target).toEqual([0, 0]);
  });

  it("uses the smaller container dimension for zoom", () => {
    const wide = { width: 1000, height: 400 };
    const tall = { width: 400, height: 1000 };
    const zoomWide = computeViewState(squareBounds, wide).zoom;
    const zoomTall = computeViewState(squareBounds, tall).zoom;
    // Both should use 400 as viewSize, so zoom should be the same
    expect(zoomWide).toBeCloseTo(zoomTall);
  });

  it("computes higher zoom for smaller data range", () => {
    const smallBounds = { minX: 0, maxX: 10, minY: 0, maxY: 10 };
    const largeBounds = { minX: 0, maxX: 1000, minY: 0, maxY: 1000 };
    const zoomSmall = computeViewState(smallBounds, squareContainer).zoom;
    const zoomLarge = computeViewState(largeBounds, squareContainer).zoom;
    expect(zoomSmall).toBeGreaterThan(zoomLarge);
  });

  it("clamps zoom to minimum of -5", () => {
    const hugeBounds = { minX: 0, maxX: 1e10, minY: 0, maxY: 1e10 };
    const { zoom } = computeViewState(hugeBounds, squareContainer);
    expect(zoom).toBe(-5);
  });

  it("clamps zoom to maximum of 10", () => {
    const tinyBounds = { minX: 0, maxX: 0.001, minY: 0, maxY: 0.001 };
    const { zoom } = computeViewState(tinyBounds, squareContainer);
    expect(zoom).toBe(10);
  });

  it("uses the larger data axis for zoom calculation", () => {
    const { zoom: zoomWide } = computeViewState(wideBounds, squareContainer);
    const { zoom: zoomSquare } = computeViewState(squareBounds, squareContainer);
    expect(zoomWide).toBeCloseTo(zoomSquare);
  });
});

// ---------------------------------------------------------------------------
// buildScatterplotPoints
// ---------------------------------------------------------------------------

// 4 points, 2 columns (x, y): row-major flat array
// Point 0: (1, 2), Point 1: (3, 4), Point 2: (5, 6), Point 3: (7, 8)
const data4x2 = new Float32Array([1, 2, 3, 4, 5, 6, 7, 8]);
const shape4x2 = [4, 2];

describe("buildScatterplotPoints", () => {
  it("returns empty result for null data", () => {
    const result = buildScatterplotPoints({ data: null, shape: [4, 2] });
    expect(result.points).toEqual([]);
    expect(result.bounds).toBeNull();
  });

  it("returns empty result for null shape", () => {
    const result = buildScatterplotPoints({ data: data4x2, shape: null });
    expect(result.points).toEqual([]);
    expect(result.bounds).toBeNull();
  });

  it("builds correct number of points", () => {
    const { points } = buildScatterplotPoints({ data: data4x2, shape: shape4x2 });
    expect(points).toHaveLength(4);
  });

  it("reads x from column 0 and flips y from column 1", () => {
    const { points } = buildScatterplotPoints({ data: data4x2, shape: shape4x2 });
    expect(points[0].position).toEqual([1, -2]);
    expect(points[1].position).toEqual([3, -4]);
    expect(points[2].position).toEqual([5, -6]);
    expect(points[3].position).toEqual([7, -8]);
  });

  it("preserves original row index on each point", () => {
    const { points } = buildScatterplotPoints({ data: data4x2, shape: shape4x2 });
    expect(points.map((p) => p.index)).toEqual([0, 1, 2, 3]);
  });

  it("computes correct bounds with flipped y", () => {
    const { bounds } = buildScatterplotPoints({ data: data4x2, shape: shape4x2 });
    expect(bounds.minX).toBe(1);
    expect(bounds.maxX).toBe(7);
    // y is flipped: -2, -4, -6, -8
    expect(bounds.minY).toBe(-8);
    expect(bounds.maxY).toBe(-2);
  });

  it("assigns category 'All' when no colorData", () => {
    const { points } = buildScatterplotPoints({ data: data4x2, shape: shape4x2 });
    for (const pt of points) {
      expect(pt.category).toBe("All");
    }
  });

  it("assigns categories from colorData", () => {
    const colorData = ["A", "B", "A", "C"];
    const { points, categoryColorMap } = buildScatterplotPoints({
      data: data4x2,
      shape: shape4x2,
      colorData,
    });
    expect(points[0].category).toBe("A");
    expect(points[1].category).toBe("B");
    expect(points[3].category).toBe("C");
    expect(Object.keys(categoryColorMap)).toEqual(["A", "B", "C"]);
  });

  it("maps category colors from CATEGORICAL_COLORS palette", () => {
    const colorData = ["X", "Y", "X", "Y"];
    const { categoryColorMap } = buildScatterplotPoints({
      data: data4x2,
      shape: shape4x2,
      colorData,
    });
    expect(categoryColorMap["X"]).toEqual(CATEGORICAL_COLORS[0]);
    expect(categoryColorMap["Y"]).toEqual(CATEGORICAL_COLORS[1]);
  });

  it("wraps color index when categories exceed palette size", () => {
    const n = CATEGORICAL_COLORS.length + 1;
    const flat = new Float32Array(n * 2);
    for (let i = 0; i < n; i++) {
      flat[i * 2] = i;
      flat[i * 2 + 1] = i;
    }
    const colorData = Array.from({ length: n }, (_, i) => `cat${i}`);
    const { categoryColorMap } = buildScatterplotPoints({
      data: flat,
      shape: [n, 2],
      colorData,
    });
    expect(categoryColorMap[`cat${CATEGORICAL_COLORS.length}`]).toEqual(CATEGORICAL_COLORS[0]);
  });

  it("attaches expression values from geneExpression", () => {
    const geneExpression = new Float32Array([0.1, 0.5, 0.9, 0.0]);
    const { points } = buildScatterplotPoints({
      data: data4x2,
      shape: shape4x2,
      geneExpression,
    });
    expect(points[0].expression).toBeCloseTo(0.1);
    expect(points[2].expression).toBeCloseTo(0.9);
  });

  it("sets expression to null when no geneExpression", () => {
    const { points } = buildScatterplotPoints({ data: data4x2, shape: shape4x2 });
    for (const pt of points) {
      expect(pt.expression).toBeNull();
    }
  });

  describe("downsampling via maxPoints", () => {
    const data10 = new Float32Array(20);
    for (let i = 0; i < 10; i++) {
      data10[i * 2] = i;
      data10[i * 2 + 1] = i * 10;
    }
    const shape10 = [10, 2];

    it("returns all points when maxPoints >= total", () => {
      const { points } = buildScatterplotPoints({ data: data10, shape: shape10, maxPoints: 100 });
      expect(points).toHaveLength(10);
    });

    it("downsamples when maxPoints < total", () => {
      const { points } = buildScatterplotPoints({ data: data10, shape: shape10, maxPoints: 5 });
      expect(points).toHaveLength(5);
      expect(points.map((p) => p.index)).toEqual([0, 2, 4, 6, 8]);
    });

    it("preserves correct positions after downsampling", () => {
      const { points } = buildScatterplotPoints({ data: data10, shape: shape10, maxPoints: 5 });
      expect(points[0].position[0]).toBe(0);
      expect(points[1].position[0]).toBe(2);
    });
  });
});

// ---------------------------------------------------------------------------
// buildSelectionSummary
// ---------------------------------------------------------------------------

const selPoints = [
  { index: 0, category: "A", expression: 1.0 },
  { index: 1, category: "B", expression: 2.0 },
  { index: 2, category: "A", expression: 3.0 },
  { index: 3, category: "C", expression: 4.0 },
  { index: 4, category: "B", expression: 5.0 },
];

const defaults = {
  points: selPoints,
  hasColorData: false,
  hasGeneExpression: false,
  tooltipData: {},
};

describe("buildSelectionSummary", () => {
  it("returns null when nothing is selected", () => {
    const result = buildSelectionSummary({ ...defaults, selectedSet: new Set() });
    expect(result).toBeNull();
  });

  describe("categoryBreakdown", () => {
    it("is null when hasColorData is false", () => {
      const result = buildSelectionSummary({
        ...defaults,
        selectedSet: new Set([0, 1]),
        hasColorData: false,
      });
      expect(result.categoryBreakdown).toBeNull();
    });

    it("counts categories for selected points only", () => {
      const result = buildSelectionSummary({
        ...defaults,
        selectedSet: new Set([0, 2, 3]),
        hasColorData: true,
      });
      expect(result.categoryBreakdown).toEqual([["A", 2], ["C", 1]]);
    });

    it("sorts categories by count descending", () => {
      const result = buildSelectionSummary({
        ...defaults,
        selectedSet: new Set([0, 1, 2, 3, 4]),
        hasColorData: true,
      });
      const counts = result.categoryBreakdown.map(([, c]) => c);
      for (let i = 1; i < counts.length; i++) {
        expect(counts[i]).toBeLessThanOrEqual(counts[i - 1]);
      }
    });
  });

  describe("expressionStats", () => {
    it("is null when hasGeneExpression is false", () => {
      const result = buildSelectionSummary({
        ...defaults,
        selectedSet: new Set([0, 1]),
        hasGeneExpression: false,
      });
      expect(result.expressionStats).toBeNull();
    });

    it("computes min, max, mean, count for selected points", () => {
      const result = buildSelectionSummary({
        ...defaults,
        selectedSet: new Set([1, 2, 3]),
        hasGeneExpression: true,
      });
      expect(result.expressionStats).toEqual({
        min: 2.0,
        max: 4.0,
        mean: 3.0,
        count: 3,
      });
    });

    it("skips points with null expression", () => {
      const mixedPoints = [
        { index: 0, category: "A", expression: null },
        { index: 1, category: "A", expression: 10 },
        { index: 2, category: "A", expression: 20 },
      ];
      const result = buildSelectionSummary({
        ...defaults,
        points: mixedPoints,
        selectedSet: new Set([0, 1, 2]),
        hasGeneExpression: true,
      });
      expect(result.expressionStats).toEqual({
        min: 10,
        max: 20,
        mean: 15,
        count: 2,
      });
    });

    it("is null when all selected expressions are null", () => {
      const nullPoints = [
        { index: 0, category: "A", expression: null },
        { index: 1, category: "A", expression: null },
      ];
      const result = buildSelectionSummary({
        ...defaults,
        points: nullPoints,
        selectedSet: new Set([0, 1]),
        hasGeneExpression: true,
      });
      expect(result.expressionStats).toBeNull();
    });
  });

  describe("tooltipBreakdowns", () => {
    it("is empty object when no tooltipData", () => {
      const result = buildSelectionSummary({
        ...defaults,
        selectedSet: new Set([0]),
        tooltipData: {},
      });
      expect(result.tooltipBreakdowns).toEqual({});
    });

    it("counts values per tooltip column for selected points", () => {
      const tooltipData = {
        tissue: ["brain", "liver", "brain", "lung", "liver"],
      };
      const result = buildSelectionSummary({
        ...defaults,
        selectedSet: new Set([0, 1, 2]),
        tooltipData,
      });
      expect(result.tooltipBreakdowns.tissue).toEqual([["brain", 2], ["liver", 1]]);
    });

    it("handles multiple tooltip columns", () => {
      const tooltipData = {
        tissue: ["brain", "liver", "brain"],
        stage: ["I", "II", "I"],
      };
      const result = buildSelectionSummary({
        ...defaults,
        points: selPoints.slice(0, 3),
        selectedSet: new Set([0, 1, 2]),
        tooltipData,
      });
      expect(Object.keys(result.tooltipBreakdowns)).toEqual(["tissue", "stage"]);
    });

    it("sorts breakdown values by count descending", () => {
      const tooltipData = {
        status: ["x", "y", "x", "x", "y"],
      };
      const result = buildSelectionSummary({
        ...defaults,
        selectedSet: new Set([0, 1, 2, 3, 4]),
        tooltipData,
      });
      expect(result.tooltipBreakdowns.status[0]).toEqual(["x", 3]);
      expect(result.tooltipBreakdowns.status[1]).toEqual(["y", 2]);
    });
  });
});
