import { describe, it, expect } from "vitest";
import { ZarrStore } from "./ZarrStore.js";

const URL = "http://localhost:3000/spectrum_all_cells.zarr";

describe("ZarrStore", () => {
  it("opens a zarr store and exposes store, root, and attrs", async () => {
    const z = await ZarrStore.open(URL);

    expect(z.store).toBeDefined();
    expect(z.root).toBeDefined();
    expect(z.attrs).toBeDefined();
  });

  it("attrs contains AnnData encoding metadata", async () => {
    const z = await ZarrStore.open(URL);

    expect(z.attrs).toHaveProperty("encoding-type", "anndata");
    expect(z.attrs).toHaveProperty("encoding-version");
  });
});
