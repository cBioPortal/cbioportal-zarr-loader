// Categorical color palette (similar to D3 category10)
export const CATEGORICAL_COLORS = [
  [31, 119, 180],   // #1f77b4
  [255, 127, 14],   // #ff7f0e
  [44, 160, 44],    // #2ca02c
  [214, 39, 40],    // #d62728
  [148, 103, 189],  // #9467bd
  [140, 86, 75],    // #8c564b
  [227, 119, 194],  // #e377c2
  [127, 127, 127],  // #7f7f7f
  [188, 189, 34],   // #bcbd22
  [23, 190, 207],   // #17becf
  [174, 199, 232],  // #aec7e8
  [255, 187, 120],  // #ffbb78
  [152, 223, 138],  // #98df8a
  [255, 152, 150],  // #ff9896
  [197, 176, 213],  // #c5b0d5
];

// Viridis color scale for continuous data
export const VIRIDIS = [
  [68, 1, 84],
  [72, 40, 120],
  [62, 74, 137],
  [49, 104, 142],
  [38, 130, 142],
  [31, 158, 137],
  [53, 183, 121],
  [109, 205, 89],
  [180, 222, 44],
  [253, 231, 37],
];

/**
 * Interpolate through the Viridis color scale.
 * @param {number} t - Value between 0 and 1
 * @returns {[number, number, number]} RGB color array
 */
export function interpolateViridis(t) {
  const clampedT = Math.max(0, Math.min(1, t));
  const idx = clampedT * (VIRIDIS.length - 1);
  const lower = Math.floor(idx);
  const upper = Math.min(lower + 1, VIRIDIS.length - 1);
  const frac = idx - lower;
  return [
    Math.round(VIRIDIS[lower][0] + (VIRIDIS[upper][0] - VIRIDIS[lower][0]) * frac),
    Math.round(VIRIDIS[lower][1] + (VIRIDIS[upper][1] - VIRIDIS[lower][1]) * frac),
    Math.round(VIRIDIS[lower][2] + (VIRIDIS[upper][2] - VIRIDIS[lower][2]) * frac),
  ];
}

/**
 * Convert RGB array to CSS rgb string.
 * @param {[number, number, number]} rgb - RGB color array
 * @returns {string} CSS rgb string
 */
export function rgbToString(rgb) {
  return `rgb(${rgb.join(",")})`;
}

/**
 * Generate a CSS gradient string from the Viridis scale.
 * @param {string} direction - CSS gradient direction (e.g., "to bottom", "to right")
 * @returns {string} CSS linear-gradient string
 */
export function viridisGradient(direction = "to bottom") {
  const colors = VIRIDIS.slice().reverse().map(c => rgbToString(c)).join(", ");
  return `linear-gradient(${direction}, ${colors})`;
}
