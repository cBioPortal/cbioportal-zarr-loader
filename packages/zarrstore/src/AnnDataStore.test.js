import { describe, it, expect } from "vitest";
import { AnnDataStore } from "./AnnDataStore.js";
import { ZarrStore } from "./ZarrStore.js";

const URL = `${globalThis.__TEST_BASE_URL__}/pbmc3k.zarr`;

describe("AnnDataStore", () => {
  describe("open", () => {
    it("opens a store and exposes shape metadata", async () => {
      const adata = await AnnDataStore.open(URL);

      expect(adata.shape).toBeDefined();
      expect(adata.shape).toHaveLength(2);
      expect(adata.nObs).toBe(adata.shape[0]);
      expect(adata.nVar).toBe(adata.shape[1]);
    });

    it("exposes root attrs with anndata encoding-type", async () => {
      const adata = await AnnDataStore.open(URL);

      expect(adata.attrs).toHaveProperty("encoding-type", "anndata");
    });

    it("exposes the underlying ZarrStore", async () => {
      const adata = await AnnDataStore.open(URL);

      expect(adata.zarrStore).toBeInstanceOf(ZarrStore);
    });
  });

  describe("fromZarrStore", () => {
    it("creates AnnDataStore from an existing ZarrStore", async () => {
      const zs = await ZarrStore.open(URL);
      const adata = await AnnDataStore.fromZarrStore(zs);

      expect(adata.shape).toBeDefined();
      expect(adata.zarrStore).toBe(zs);
    });
  });

  describe("X matrix", () => {
    it("reads the full X matrix", async () => {
      const adata = await AnnDataStore.open(URL);
      const x = await adata.X();

      expect(x).toBeDefined();
      // sparse will have format/data/indices/indptr/shape
      // dense will have data/shape
      expect(x.data).toBeDefined();
      expect(x.shape).toBeDefined();
    });
  });

  describe("obs", () => {
    it("reads the full obs dataframe", async () => {
      const adata = await AnnDataStore.open(URL);
      const obs = await adata.obs();

      expect(obs.index).toBeDefined();
      expect(obs.columns).toBeDefined();
      expect(obs.columnOrder).toBeDefined();
      expect(Array.isArray(obs.index)).toBe(true);
      expect(obs.index.length).toBe(adata.nObs);
    });

    it("reads obs column names", async () => {
      const adata = await AnnDataStore.open(URL);
      const cols = await adata.obsColumns();

      expect(Array.isArray(cols)).toBe(true);
      expect(cols.length).toBeGreaterThan(0);
    });

    it("reads a single obs column", async () => {
      const adata = await AnnDataStore.open(URL);
      const cols = await adata.obsColumns();
      const col = await adata.obsColumn(cols[0]);

      expect(col).toBeDefined();
      expect(col.length).toBe(adata.nObs);
    });

    it("reads obs names (index)", async () => {
      const adata = await AnnDataStore.open(URL);
      const names = await adata.obsNames();

      expect(Array.isArray(names)).toBe(true);
      expect(names.length).toBe(adata.nObs);
    });
  });

  describe("var", () => {
    it("reads the full var dataframe", async () => {
      const adata = await AnnDataStore.open(URL);
      const v = await adata.var();

      expect(v.index).toBeDefined();
      expect(v.columns).toBeDefined();
      expect(v.columnOrder).toBeDefined();
      expect(Array.isArray(v.index)).toBe(true);
      expect(v.index.length).toBe(adata.nVar);
    });

    it("reads var names (gene names)", async () => {
      const adata = await AnnDataStore.open(URL);
      const names = await adata.varNames();

      expect(Array.isArray(names)).toBe(true);
      expect(names.length).toBe(adata.nVar);
    });
  });

  describe("obsm", () => {
    it("lists obsm keys", async () => {
      const adata = await AnnDataStore.open(URL);
      const keys = adata.obsmKeys();

      expect(Array.isArray(keys)).toBe(true);
    });

    it("reads an obsm entry if keys exist", async () => {
      const adata = await AnnDataStore.open(URL);
      const keys = adata.obsmKeys();

      if (keys.length > 0) {
        const entry = await adata.obsm(keys[0]);
        expect(entry).toBeDefined();
        expect(entry.data).toBeDefined();
      }
    });
  });

  describe("layers", () => {
    it("lists layer keys", async () => {
      const adata = await AnnDataStore.open(URL);
      const keys = adata.layerKeys();

      expect(Array.isArray(keys)).toBe(true);
    });
  });

  describe("uns", () => {
    it("lists uns keys", async () => {
      const adata = await AnnDataStore.open(URL);
      const keys = adata.unsKeys();

      expect(Array.isArray(keys)).toBe(true);
    });
  });
});
