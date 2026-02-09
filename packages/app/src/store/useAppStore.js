import { create } from "zustand";
import { AnnDataStore } from "@cbioportal-zarr-loader/zarrstore";

const useAppStore = create((set, get) => ({
  // Core data
  adata: null,
  metadata: null,
  loading: true,
  error: null,

  // Obs column state (multi-select)
  obsColumnsSelected: [],
  obsColumnsData: {},
  obsColumnLoading: null,
  obsColumnTime: null,

  // Var column state (multi-select)
  varColumnsSelected: [],
  varColumnsData: {},
  varColumnLoading: null,
  varColumnTime: null,

  // Obsm state
  selectedObsm: null,
  obsmData: null,
  obsmLoading: false,
  obsmTime: null,

  // Obsm streaming state
  obsmStreamingData: null,
  obsmStreamingLoading: false,
  obsmStreamingTime: null,
  obsmStreamingProgress: null,

  // Gene expression state (for scatterplot coloring)
  selectedGene: null,
  geneExpression: null,
  geneLoading: false,

  // Obs column for scatterplot coloring
  colorColumn: null,
  colorData: null,
  colorLoading: false,

  // Tooltip obs columns for scatterplot
  tooltipColumns: [],
  tooltipData: {},
  tooltipColumnLoading: null,

  // Plots tab state (independent of scatterplot coloring)
  plotGene: null,
  plotGeneExpression: null,
  plotGeneLoading: false,
  plotObsColumn: null,
  plotObsData: null,
  plotObsLoading: false,

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
      const obsNames = await store.obsNames();
      timings.obsNames = performance.now() - start;

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
        obsIndex: obsNames,
        varIndex: varNames,
        loading: false,
        error: null,
      });
    } catch (err) {
      set({ error: err.message, loading: false });
      console.error(err);
    }
  },

  toggleObsColumn: async (colName) => {
    const { adata, obsColumnsSelected, obsColumnsData } = get();
    if (!adata) return;

    // Toggle off — remove column
    if (obsColumnsSelected.includes(colName)) {
      const { [colName]: _, ...rest } = obsColumnsData;
      set({
        obsColumnsSelected: obsColumnsSelected.filter((c) => c !== colName),
        obsColumnsData: rest,
      });
      return;
    }

    // Toggle on — fetch and add column
    set({ obsColumnLoading: colName });

    try {
      const start = performance.now();
      const values = await adata.obsColumn(colName);
      const { obsColumnsSelected: current, obsColumnsData: currentData } = get();
      set({
        obsColumnTime: performance.now() - start,
        obsColumnsSelected: [...current, colName],
        obsColumnsData: { ...currentData, [colName]: values },
        obsColumnLoading: null,
      });
    } catch (err) {
      console.error(err);
      set({ obsColumnLoading: null });
    }
  },

  toggleVarColumn: async (colName) => {
    const { adata, varColumnsSelected, varColumnsData } = get();
    if (!adata) return;

    // Toggle off — remove column
    if (varColumnsSelected.includes(colName)) {
      const { [colName]: _, ...rest } = varColumnsData;
      set({
        varColumnsSelected: varColumnsSelected.filter((c) => c !== colName),
        varColumnsData: rest,
      });
      return;
    }

    // Toggle on — fetch and add column
    set({ varColumnLoading: colName });

    try {
      const start = performance.now();
      const values = await adata.varColumn(colName);
      const { varColumnsSelected: current, varColumnsData: currentData } = get();
      set({
        varColumnTime: performance.now() - start,
        varColumnsSelected: [...current, colName],
        varColumnsData: { ...currentData, [colName]: values },
        varColumnLoading: null,
      });
    } catch (err) {
      console.error(err);
      set({ varColumnLoading: null });
    }
  },

  clearObsColumns: () => {
    set({ obsColumnsSelected: [], obsColumnsData: {}, obsColumnTime: null });
  },

  clearVarColumns: () => {
    set({ varColumnsSelected: [], varColumnsData: {}, varColumnTime: null });
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

  fetchObsmStreaming: async (key) => {
    const { adata } = get();
    if (!adata) return;

    set({
      obsmStreamingData: null,
      obsmStreamingLoading: true,
      obsmStreamingTime: null,
      obsmStreamingProgress: 0,
    });

    try {
      const start = performance.now();
      let buffer = null;
      let nDims = 0;
      let loadedRows = 0;

      for await (const batch of adata.obsmStreaming(key)) {
        const { data, shape, offset, total } = batch;
        nDims = shape[1];

        // Pre-allocate on first batch
        if (!buffer) {
          buffer = new data.constructor(total * nDims);
        }

        // Copy batch data into buffer
        buffer.set(data, offset * nDims);
        loadedRows = offset + shape[0];

        set({
          obsmStreamingData: { data: buffer, shape: [loadedRows, nDims] },
          obsmStreamingLoading: false,
          obsmStreamingProgress: loadedRows / total,
        });
      }

      set({
        obsmStreamingTime: performance.now() - start,
        obsmStreamingProgress: null,
      });
    } catch (err) {
      console.error(err);
      set({
        obsmStreamingData: { error: err.message },
        obsmStreamingLoading: false,
        obsmStreamingProgress: null,
      });
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

  toggleTooltipColumn: async (colName) => {
    const { adata, tooltipColumns, tooltipData } = get();
    if (!adata) return;

    if (tooltipColumns.includes(colName)) {
      const { [colName]: _, ...rest } = tooltipData;
      set({
        tooltipColumns: tooltipColumns.filter((c) => c !== colName),
        tooltipData: rest,
      });
      return;
    }

    set({ tooltipColumnLoading: colName });
    try {
      const values = await adata.obsColumn(colName);
      const { tooltipColumns: current, tooltipData: currentData } = get();
      set({
        tooltipColumns: [...current, colName],
        tooltipData: { ...currentData, [colName]: values },
        tooltipColumnLoading: null,
      });
    } catch (err) {
      console.error(err);
      set({ tooltipColumnLoading: null });
    }
  },

  clearTooltipColumns: () => {
    set({ tooltipColumns: [], tooltipData: {} });
  },

  // Plots tab actions (independent of scatterplot coloring)
  setPlotGene: async (geneName) => {
    const { adata, metadata } = get();
    if (!adata || !geneName) return;

    console.log("[PlotsTab] setPlotGene called:", geneName);
    set({ plotGene: geneName, plotGeneLoading: true, plotGeneExpression: null });

    try {
      const { varNames, geneNames } = metadata;

      let queryName = geneName;
      if (geneNames !== varNames) {
        const idx = geneNames.indexOf(geneName);
        if (idx !== -1) {
          queryName = varNames[idx];
        }
      }

      console.log("[PlotsTab] Fetching gene expression for queryName:", queryName);
      const values = await adata.geneExpression(queryName);
      console.log("[PlotsTab] Gene expression fetched, length:", values?.length, "sample:", values?.slice(0, 5));
      set({ plotGeneExpression: values, plotGeneLoading: false });
    } catch (err) {
      console.error("[PlotsTab] Gene expression fetch error:", err);
      set({ plotGeneExpression: null, plotGeneLoading: false });
    }
  },

  clearPlotGene: () => {
    set({ plotGene: null, plotGeneExpression: null });
  },

  setPlotObsColumn: async (colName) => {
    const { adata } = get();
    if (!adata || !colName) return;

    console.log("[PlotsTab] setPlotObsColumn called:", colName);
    set({ plotObsColumn: colName, plotObsLoading: true, plotObsData: null });

    try {
      console.log("[PlotsTab] Fetching obs column:", colName);
      const values = await adata.obsColumn(colName);
      console.log("[PlotsTab] Obs column fetched, length:", values?.length, "sample:", values?.slice(0, 5));
      set({ plotObsData: values, plotObsLoading: false });
    } catch (err) {
      console.error("[PlotsTab] Obs column fetch error:", err);
      set({ plotObsData: null, plotObsLoading: false });
    }
  },

  clearPlotObsColumn: () => {
    set({ plotObsColumn: null, plotObsData: null });
  },
}));

export default useAppStore;
