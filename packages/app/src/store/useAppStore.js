import { create } from "zustand";
import { AnnDataStore } from "@cbioportal-zarr-loader/zarrstore";
import {
  FilterSchema,
  findMatchingIndices,
  resolveInitialView,
  resolveViewWithDefaults,
} from "../utils/filterUtils";
import {
  pointInPolygon,
  buildScatterplotPoints,
} from "../utils/scatterplotUtils";
import { fetchFeatureFlags } from "../utils/featureFlags";

const useAppStore = create((set, get) => ({
  // Feature flags (fetched from GitHub, local fallback)
  featureFlags: {},

  // Core data
  url: null,
  adata: null,
  metadata: null,
  loading: true,
  error: null,
  pendingFilterConfig: null,

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
  colorScaleName: "viridis",

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

  // Dotplot tab state
  dotplotGenes: [],
  dotplotGeneExpressions: {},
  dotplotGeneLoading: null,

  // Selection state
  selectedPointIndices: [],
  selectionGeometry: null, // { type: "rectangle", bounds: [x1,y1,x2,y2] } or { type: "lasso", polygon: [[x,y],...] }

  // View config state (from JSON config or postMessage)
  filterJson: "",
  viewConfigDefaults: {},
  appliedSelections: [],
  activeSelectionIndex: undefined,

  // Cached indices
  obsIndex: null,
  varIndex: null,

  // Actions
  initialize: async (url) => {
    const timings = {};

    try {
      let start = performance.now();
      const [store, featureFlags] = await Promise.all([
        AnnDataStore.open(url),
        fetchFeatureFlags(),
      ]);
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
        url,
        adata: store,
        metadata: { obsColumns, varColumns, obsmKeys, layerKeys, varNames, geneNames, timings, chunks },
        obsIndex: obsNames,
        varIndex: varNames,
        featureFlags,
        loading: false,
        error: null,
      });

      // Drain queued config that arrived before initialization completed
      const { pendingFilterConfig } = get();
      if (pendingFilterConfig) {
        console.debug("[CZL:postMessage] Applying queued config after initialization");
        set({ pendingFilterConfig: null });
        await get().applyFilterConfig(pendingFilterConfig);
      }
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

  setColorScaleName: (name) => set({ colorScaleName: name }),

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

    console.debug("[PlotsTab] setPlotGene called:", geneName);
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

      console.debug("[PlotsTab] Fetching gene expression for queryName:", queryName);
      const values = await adata.geneExpression(queryName);
      console.debug("[PlotsTab] Gene expression fetched, length:", values?.length, "sample:", values?.slice(0, 5));
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

    console.debug("[PlotsTab] setPlotObsColumn called:", colName);
    set({ plotObsColumn: colName, plotObsLoading: true, plotObsData: null });

    try {
      console.debug("[PlotsTab] Fetching obs column:", colName);
      const values = await adata.obsColumn(colName);
      console.debug("[PlotsTab] Obs column fetched, length:", values?.length, "sample:", values?.slice(0, 5));
      set({ plotObsData: values, plotObsLoading: false });
    } catch (err) {
      console.error("[PlotsTab] Obs column fetch error:", err);
      set({ plotObsData: null, plotObsLoading: false });
    }
  },

  clearPlotObsColumn: () => {
    set({ plotObsColumn: null, plotObsData: null });
  },

  // Dotplot tab actions
  toggleDotplotGene: async (geneName) => {
    const { adata, metadata, dotplotGenes, dotplotGeneExpressions } = get();
    if (!adata || !geneName) return;

    // Toggle off — remove gene
    if (dotplotGenes.includes(geneName)) {
      const { [geneName]: _, ...rest } = dotplotGeneExpressions;
      set({
        dotplotGenes: dotplotGenes.filter((g) => g !== geneName),
        dotplotGeneExpressions: rest,
      });
      return;
    }

    // Toggle on — fetch and add gene
    set({ dotplotGeneLoading: geneName });

    try {
      const { varNames, geneNames } = metadata;
      let queryName = geneName;
      if (geneNames !== varNames) {
        const idx = geneNames.indexOf(geneName);
        if (idx !== -1) queryName = varNames[idx];
      }
      const values = await adata.geneExpression(queryName);
      const { dotplotGenes: current, dotplotGeneExpressions: currentData } = get();
      set({
        dotplotGenes: [...current, geneName],
        dotplotGeneExpressions: { ...currentData, [geneName]: values },
        dotplotGeneLoading: null,
      });
    } catch (err) {
      console.error("[DotplotTab] Gene expression fetch error:", err);
      set({ dotplotGeneLoading: null });
    }
  },

  clearDotplotGenes: () => {
    set({ dotplotGenes: [], dotplotGeneExpressions: {} });
  },

  setSelectedPoints: (indices) => {
    set({ selectedPointIndices: indices });
  },

  setSelectionGeometry: (geometry) => set({ selectionGeometry: geometry }),

  clearSelectedPoints: () => {
    set({ selectedPointIndices: [], selectionGeometry: null });
  },

  // Apply a single resolved view (embedding, tooltips, color, selection)
  applyView: async (view) => {
    const defaults = get().viewConfigDefaults;
    const resolved = resolveViewWithDefaults(view, defaults);
    const { adata } = get();

    // Load embedding (skip if already on the same key)
    if (resolved.embeddingKey && resolved.embeddingKey !== get().selectedObsm) {
      await get().fetchObsm(resolved.embeddingKey);
    }

    const selection = resolved.selection;
    const selType = selection.type || "category";

    // Build the set of obs columns to fetch in parallel
    const tooltipCols = selType === "category"
      ? [...new Set([selection.target, ...resolved.activeTooltips])]
      : [...resolved.activeTooltips];

    // Determine color_by fetch
    const colorByCategoryCol = resolved.colorBy?.type === "category" ? resolved.colorBy.value : null;
    let geneQueryName = null;
    if (resolved.colorBy?.type === "gene") {
      const { geneNames, varNames } = get().metadata;
      const match = geneNames.find(g => g.toLowerCase() === resolved.colorBy.value.toLowerCase());
      if (match) {
        geneQueryName = match;
        if (geneNames !== varNames) {
          const idx = geneNames.indexOf(match);
          if (idx !== -1) geneQueryName = varNames[idx];
        }
      }
    }

    // Fetch all obs columns + color/gene in parallel
    const allObsCols = [...new Set([...tooltipCols, ...(colorByCategoryCol ? [colorByCategoryCol] : [])])];
    const fetches = allObsCols.map(col => adata.obsColumn(col).catch(() => null));
    if (geneQueryName) {
      fetches.push(adata.geneExpression(geneQueryName).catch(() => null));
    }

    const results = await Promise.all(fetches);

    // Unpack obs column results
    const obsResults = {};
    allObsCols.forEach((col, i) => { obsResults[col] = results[i]; });

    // Apply tooltips — single store write
    const newTooltipData = {};
    for (const col of tooltipCols) {
      if (obsResults[col]) newTooltipData[col] = obsResults[col];
    }
    set({ tooltipColumns: tooltipCols, tooltipData: newTooltipData, tooltipColumnLoading: null });

    // Apply color_by
    if (colorByCategoryCol && obsResults[colorByCategoryCol]) {
      set({
        colorColumn: colorByCategoryCol,
        colorData: obsResults[colorByCategoryCol],
        colorLoading: false,
        selectedGene: null,
        geneExpression: null,
      });
    } else if (geneQueryName) {
      const geneValues = results[results.length - 1];
      if (geneValues) {
        set({
          selectedGene: resolved.colorBy.value,
          geneExpression: geneValues,
          geneLoading: false,
          colorColumn: null,
          colorData: null,
        });
      }
    }
    if (resolved.colorBy?.color_scale) {
      set({ colorScaleName: resolved.colorBy.color_scale });
    }

    // Apply selection
    if (selType === "category") {
      const columnData = obsResults[selection.target];
      if (!columnData) return;
      const matchingIndices = findMatchingIndices(columnData, selection.values);
      set({ selectionGeometry: null, selectedPointIndices: matchingIndices });
    } else if (selType === "rectangle" || selType === "lasso") {
      const currentObsm = get().obsmData;
      if (!currentObsm?.data || !currentObsm?.shape) return;

      const { points } = buildScatterplotPoints({
        data: currentObsm.data,
        shape: currentObsm.shape,
      });

      const indices = [];
      if (selType === "rectangle") {
        const [minX, minY, maxX, maxY] = selection.bounds;
        for (const pt of points) {
          const [px, py] = pt.position;
          if (px >= minX && px <= maxX && py >= minY && py <= maxY) {
            indices.push(pt.index);
          }
        }
        set({ selectionGeometry: { type: "rectangle", bounds: selection.bounds } });
      } else {
        for (const pt of points) {
          const [px, py] = pt.position;
          if (pointInPolygon(px, py, selection.polygon)) {
            indices.push(pt.index);
          }
        }
        set({ selectionGeometry: { type: "lasso", polygon: selection.polygon } });
      }

      set({ selectedPointIndices: indices });
    }
  },

  // Validate raw config, resolve initial view, apply it, and populate selections
  applyFilterConfig: async (raw) => {
    try {
      const { adata, loading } = get();
      if (!adata || loading) {
        console.debug("[CZL:postMessage] Store not ready, queuing config for after initialization");
        set({ pendingFilterConfig: raw });
        return { success: true, queued: true };
      }

      const result = FilterSchema.safeParse(raw);
      if (!result.success) {
        const errorMsg = result.error.issues.map(i => i.message).join("; ");
        return { success: false, error: errorMsg };
      }

      const { defaults: parsedDefaults = {}, initial_view, saved_views } = result.data;
      set({ viewConfigDefaults: parsedDefaults, filterJson: JSON.stringify(raw, null, 2) });

      const initialMatch = resolveInitialView(initial_view, saved_views);
      if (!initialMatch) {
        const msg = typeof initial_view === "number"
          ? `Index ${initial_view} out of range (${saved_views.length} saved views, use 0-${saved_views.length - 1})`
          : `View "${initial_view}" not found in saved_views`;
        return { success: false, error: msg };
      }

      await get().applyView(initialMatch);

      set({
        appliedSelections: saved_views,
        activeSelectionIndex: saved_views.indexOf(initialMatch),
      });

      return { success: true };
    } catch (err) {
      return { success: false, error: err.message };
    }
  },

  setActiveSelectionIndex: (index) => set({ activeSelectionIndex: index }),

  setAppliedSelections: (selections) => set({ appliedSelections: selections }),

  setFilterJson: (json) => set({ filterJson: json }),
}));

export default useAppStore;
