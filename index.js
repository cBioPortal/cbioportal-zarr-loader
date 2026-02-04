import * as zarr from "zarrita";
import { ZarrStore } from "./src/ZarrStore.js";

const z = await ZarrStore.open(
  "http://localhost:3000/spectrum_all_cells.zarr",
);
console.log("Root attributes:", z.attrs);

// Open the X array (expression matrix)
const X = await zarr.open(z.root.resolve("X"), { kind: "array" });
console.log("\nX array:");
console.log("  Shape:", X.shape);
console.log("  Chunks:", X.chunks);
console.log("  Dtype:", X.dtype);

// Read a small slice: first 5 cells x first 10 genes
const slice = await zarr.get(X, [zarr.slice(5), zarr.slice(10)]);
console.log("\nX[0:5, 0:10]:");
console.log("  Shape:", slice.shape);
console.log("  Data:", slice.data);
