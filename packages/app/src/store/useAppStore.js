import { create } from "zustand";
import { AnnDataStore } from "@cbioportal-zarr-loader/zarrstore";

const useAppStore = create((set, get) => ({
  // Core data
  adata: null,
  metadata: null,
  loading: true,
  error: null,

  // Obs column state
  selectedObsColumn: null,
  obsColumnData: null,
  obsColumnLoading: false,
  obsColumnTime: null,

  // Var column state
  selectedVarColumn: null,
  varColumnData: null,
  varColumnLoading: false,
  varColumnTime: null,

  // Obsm state
  selectedObsm: null,
  obsmData: null,
  obsmLoading: false,
  obsmTime: null,

  // Gene expression state (for scatterplot coloring)
  selectedGene: null,
  geneExpression: null,
  geneLoading: false,

  // Obs column for scatterplot coloring
  colorColumn: null,
  colorData: null,
  colorLoading: false,

  // Cached indices
  obsIndex: null,
  varIndex: null,

  // Actions
  initialize: async (url) => {
    const timings = {};

    try {
      let start = performance.now();
      const store = await AnnDataStore.open(url);
      timings.open = performance.now() - start;

      start = performance.now();
      const obsColumns = await store.obsColumns();
      timings.obsColumns = performance.now() - start;

      start = performance.now();
      const varColumns = await store.varColumns();
      timings.varColumns = performance.now() - start;

      start = performance.now();
      const obsmKeys = store.obsmKeys();
      timings.obsmKeys = performance.now() - start;

      start = performance.now();
      const layerKeys = store.layerKeys();
      timings.layerKeys = performance.now() - start;

      start = performance.now();
      const varNames = await store.varNames();
      timings.varNames = performance.now() - start;

      // Try to get feature_name column for display
      let geneNames = varNames;
      try {
        const featureNames = await store.varColumn("feature_name");
        if (featureNames && featureNames.length === varNames.length) {
          geneNames = featureNames;
        }
      } catch {
        // feature_name column doesn't exist
      }

      // Get chunk size from X array
      let chunks = null;
      try {
        const xArray = await store.zarrStore.openArray("X");
        chunks = xArray.chunks;
      } catch {
        try {
          const dataArray = await store.zarrStore.openArray("X/data");
          chunks = dataArray.chunks;
        } catch {
          // Couldn't get chunks
        }
      }

      set({
        adata: store,
        metadata: { obsColumns, varColumns, obsmKeys, layerKeys, varNames, geneNames, timings, chunks },
        loading: false,
        error: null,
      });
    } catch (err) {
      set({ error: err.message, loading: false });
      console.error(err);
    }
  },

  fetchObsColumn: async (colName) => {
    const { adata, obsIndex } = get();
    if (!adata) return;

    set({ selectedObsColumn: colName, obsColumnLoading: true, obsColumnData: null });

    try {
      const start = performance.now();
      let index = obsIndex;
      if (!index) {
        index = await adata.obsNames();
        set({ obsIndex: index });
      }
      const values = await adata.obsColumn(colName);
      set({
        obsColumnTime: performance.now() - start,
        obsColumnData: { values, index },
        obsColumnLoading: false,
      });
    } catch (err) {
      console.error(err);
      set({ obsColumnData: { error: err.message }, obsColumnLoading: false });
    }
  },

  fetchVarColumn: async (colName) => {
    const { adata, varIndex } = get();
    if (!adata) return;

    set({ selectedVarColumn: colName, varColumnLoading: true, varColumnData: null });

    try {
      const start = performance.now();
      let index = varIndex;
      if (!index) {
        index = await adata.varNames();
        set({ varIndex: index });
      }
      const values = await adata.varColumn(colName);
      set({
        varColumnTime: performance.now() - start,
        varColumnData: { values, index },
        varColumnLoading: false,
      });
    } catch (err) {
      console.error(err);
      set({ varColumnData: { error: err.message }, varColumnLoading: false });
    }
  },

  fetchObsm: async (key) => {
    const { adata, obsIndex } = get();
    if (!adata) return;

    set({ selectedObsm: key, obsmLoading: true, obsmData: null });

    try {
      const start = performance.now();
      let index = obsIndex;
      if (!index) {
        index = await adata.obsNames();
        set({ obsIndex: index });
      }
      const result = await adata.obsm(key);
      set({
        obsmTime: performance.now() - start,
        obsmData: { ...result, index },
        obsmLoading: false,
      });
    } catch (err) {
      console.error(err);
      set({ obsmData: { error: err.message }, obsmLoading: false });
    }
  },

  // Scatterplot coloring actions
  setColorColumn: async (colName) => {
    const { adata } = get();

    if (!colName) {
      set({ colorColumn: null, colorData: null });
      return;
    }

    // Clear gene selection when obs column is selected
    set({
      colorColumn: colName,
      colorLoading: true,
      colorData: null,
      selectedGene: null,
      geneExpression: null,
    });

    try {
      const values = await adata.obsColumn(colName);
      set({ colorData: values, colorLoading: false });
    } catch (err) {
      console.error(err);
      set({ colorData: null, colorLoading: false });
    }
  },

  setSelectedGene: async (geneName) => {
    const { adata, metadata } = get();

    if (!geneName) {
      set({ selectedGene: null, geneExpression: null });
      return;
    }

    // Clear obs column when gene is selected
    set({
      selectedGene: geneName,
      geneLoading: true,
      geneExpression: null,
      colorColumn: null,
      colorData: null,
    });

    try {
      const { varNames, geneNames } = metadata;

      // Map from feature_name to var index if needed
      let queryName = geneName;
      if (geneNames !== varNames) {
        const idx = geneNames.indexOf(geneName);
        if (idx !== -1) {
          queryName = varNames[idx];
        }
      }

      const values = await adata.geneExpression(queryName);
      set({ geneExpression: values, geneLoading: false });
    } catch (err) {
      console.error(err);
      set({ geneExpression: null, geneLoading: false });
    }
  },

  clearGeneSelection: () => {
    set({ selectedGene: null, geneExpression: null });
  },
}));

export default useAppStore;
