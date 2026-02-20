import { ascending, quantileSorted } from "d3-array";

/**
 * Compute boxplot statistics for each category group.
 *
 * @param {Array<Object>} data - Raw per-cell data objects (e.g. [{ cell_type: "A", EGFR: 1.2 }, ...])
 * @param {string} categoryField - Key used for grouping (e.g. "cell_type")
 * @param {string} valueField - Key for the numeric value (e.g. "EGFR")
 * @returns {{ groups: string[], stats: Array<{ group: string, min: number, q1: number, median: number, q3: number, max: number, whiskerLow: number, whiskerHigh: number, outliers: number[], count: number }> }}
 */
export function computeBoxplotStats(data, categoryField, valueField) {
  if (!data || data.length === 0) {
    return { groups: [], stats: [] };
  }

  // Group values by category
  const grouped = new Map();
  for (const d of data) {
    const key = d[categoryField];
    if (key == null) continue;
    const val = d[valueField];
    if (val == null || Number.isNaN(val)) continue;
    if (!grouped.has(key)) grouped.set(key, []);
    grouped.get(key).push(val);
  }

  const groups = [...grouped.keys()].sort();
  const stats = [];

  for (const group of groups) {
    const values = grouped.get(group);
    values.sort(ascending);

    const count = values.length;
    const min = values[0];
    const max = values[count - 1];
    const q1 = quantileSorted(values, 0.25);
    const median = quantileSorted(values, 0.5);
    const q3 = quantileSorted(values, 0.75);
    const iqr = q3 - q1;
    const whiskerLow = Math.max(min, q1 - 1.5 * iqr);
    const whiskerHigh = Math.min(max, q3 + 1.5 * iqr);

    // Outliers are values outside the whisker bounds
    const outliers = values.filter((v) => v < whiskerLow || v > whiskerHigh);

    stats.push({
      group,
      min,
      q1,
      median,
      q3,
      max,
      whiskerLow,
      whiskerHigh,
      outliers,
      count,
    });
  }

  return { groups, stats };
}
