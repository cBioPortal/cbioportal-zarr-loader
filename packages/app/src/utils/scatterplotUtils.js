import { CATEGORICAL_COLORS } from "./colors";

/**
 * Ray-casting point-in-polygon test.
 * Returns true if the point (x, y) lies inside the polygon.
 *
 * @param {number} x
 * @param {number} y
 * @param {Array<[number, number]>} polygon - Array of [x, y] vertices
 * @returns {boolean}
 */
export function pointInPolygon(x, y, polygon) {
  let inside = false;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i][0], yi = polygon[i][1];
    const xj = polygon[j][0], yj = polygon[j][1];
    if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
      inside = !inside;
    }
  }
  return inside;
}

/**
 * Compute the min and max of a numeric array.
 *
 * @param {Float32Array|number[]|null} values
 * @returns {{ min: number, max: number } | null} null if values is falsy or empty
 */
export function computeRange(values) {
  if (!values || values.length === 0) return null;
  let min = Infinity, max = -Infinity;
  for (let i = 0; i < values.length; i++) {
    const val = values[i];
    if (val < min) min = val;
    if (val > max) max = val;
  }
  return { min, max };
}

/**
 * Compute the initial view state (center target and zoom level) for an
 * orthographic scatterplot given data bounds and container dimensions.
 *
 * @param {Object|null} bounds - { minX, maxX, minY, maxY }
 * @param {{ width: number, height: number }} containerSize
 * @returns {{ target: [number, number], zoom: number }}
 */
export function computeViewState(bounds, containerSize) {
  if (!bounds) return { target: [0, 0], zoom: 0 };

  const centerX = (bounds.minX + bounds.maxX) / 2;
  const centerY = (bounds.minY + bounds.maxY) / 2;
  const rangeX = bounds.maxX - bounds.minX;
  const rangeY = bounds.maxY - bounds.minY;
  const maxRange = Math.max(rangeX, rangeY);

  const viewSize = Math.min(containerSize.width, containerSize.height);
  const zoom = Math.log2(viewSize / maxRange) - 0.1;

  return {
    target: [centerX, centerY],
    zoom: Math.max(-5, Math.min(zoom, 10)),
  };
}

/**
 * Build scatter plot points, category color map, and bounds from raw embedding data.
 *
 * @param {Object} options
 * @param {Float32Array|number[]} options.data - Flat array of embedding values (row-major)
 * @param {[number, number]} options.shape - [nRows, nCols]
 * @param {number} [options.maxPoints=Infinity] - Downsample to at most this many points
 * @param {ArrayLike<*>|null} [options.colorData=null] - Per-row category values for coloring
 * @param {Float32Array|number[]|null} [options.geneExpression=null] - Per-row expression values
 * @returns {{ points: Array, categoryColorMap: Object, bounds: Object|null }}
 */
export function buildScatterplotPoints({
  data,
  shape,
  maxPoints = Infinity,
  colorData = null,
  geneExpression = null,
}) {
  if (!data || !shape) return { points: [], categoryColorMap: {}, bounds: null };

  const cols = shape[1];
  const step = Math.max(1, Math.floor(shape[0] / maxPoints));

  // Build category color map (for obs columns)
  const categories = new Map();
  if (colorData) {
    for (let i = 0; i < shape[0]; i += step) {
      const cat = String(colorData[i]);
      if (!categories.has(cat)) {
        categories.set(cat, (categories.size % CATEGORICAL_COLORS.length));
      }
    }
  }

  // Build points array and calculate bounds
  const pts = [];
  let minX = Infinity, maxX = -Infinity;
  let minY = Infinity, maxY = -Infinity;

  for (let i = 0; i < shape[0]; i += step) {
    const x = data[i * cols];
    const y = -data[i * cols + 1]; // Flip Y axis

    minX = Math.min(minX, x);
    maxX = Math.max(maxX, x);
    minY = Math.min(minY, y);
    maxY = Math.max(maxY, y);

    const cat = colorData ? String(colorData[i]) : "All";
    const exprValue = geneExpression ? geneExpression[i] : null;

    pts.push({
      position: [x, y],
      category: cat,
      colorIndex: categories.get(cat) ?? 0,
      index: i,
      expression: exprValue,
    });
  }

  // Convert categories map to object for legend
  const colorMap = {};
  categories.forEach((colorIdx, cat) => {
    colorMap[cat] = CATEGORICAL_COLORS[colorIdx];
  });

  return {
    points: pts,
    categoryColorMap: colorMap,
    bounds: { minX, maxX, minY, maxY },
  };
}

/**
 * Compute a summary of the currently selected points, including
 * category breakdowns, expression statistics, and tooltip column breakdowns.
 *
 * @param {Object} options
 * @param {Set<number>} options.selectedSet - Set of selected point indices
 * @param {Array<{index: number, category: string, expression: number|null}>} options.points
 * @param {boolean} options.hasColorData - Whether color column data is active
 * @param {boolean} options.hasGeneExpression - Whether gene expression data is active
 * @param {Object<string, ArrayLike>} options.tooltipData - Map of column name to per-row values
 * @returns {null|{categoryBreakdown: Array|null, expressionStats: Object|null, tooltipBreakdowns: Object}}
 */
export function buildSelectionSummary({
  selectedSet,
  points,
  hasColorData,
  hasGeneExpression,
  tooltipData,
}) {
  if (selectedSet.size === 0) return null;

  // Category breakdown
  let categoryBreakdown = null;
  if (hasColorData) {
    const counts = {};
    for (const pt of points) {
      if (!selectedSet.has(pt.index)) continue;
      counts[pt.category] = (counts[pt.category] || 0) + 1;
    }
    categoryBreakdown = Object.entries(counts).sort((a, b) => b[1] - a[1]);
  }

  // Expression stats
  let expressionStats = null;
  if (hasGeneExpression) {
    let min = Infinity, max = -Infinity, sum = 0, count = 0;
    for (const pt of points) {
      if (!selectedSet.has(pt.index) || pt.expression == null) continue;
      const v = pt.expression;
      if (v < min) min = v;
      if (v > max) max = v;
      sum += v;
      count++;
    }
    if (count > 0) {
      expressionStats = { min, max, mean: sum / count, count };
    }
  }

  // Tooltip column breakdowns
  const tooltipBreakdowns = {};
  for (const [col, values] of Object.entries(tooltipData)) {
    const counts = {};
    for (const pt of points) {
      if (!selectedSet.has(pt.index)) continue;
      const val = String(values[pt.index]);
      counts[val] = (counts[val] || 0) + 1;
    }
    tooltipBreakdowns[col] = Object.entries(counts).sort((a, b) => b[1] - a[1]);
  }

  return { categoryBreakdown, expressionStats, tooltipBreakdowns };
}
