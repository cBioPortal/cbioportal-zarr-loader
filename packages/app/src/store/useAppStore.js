import { create } from "zustand";
import { AnnDataStore } from "@cbioportal-cell-explorer/zarrstore";
import { fetchFeatureFlags } from "../utils/featureFlags";

import { createColumnsSlice } from "./slices/columnsSlice";
import { createEmbeddingSlice } from "./slices/embeddingSlice";
import { createPlotsSlice } from "./slices/plotsSlice";
import { createDotplotSlice } from "./slices/dotplotSlice";
import { createConfigSlice } from "./slices/configSlice";

const createCoreSlice = (set, get) => ({
  // Embedded mode (iframe or ?embedded query param) — evaluated once at store creation
  isEmbedded: window.self !== window.top ||
    new URLSearchParams(window.location.search).has("embedded"),

  // Feature flags (fetched from GitHub, local fallback)
  featureFlags: {},

  // Core data
  url: null,
  adata: null,
  metadata: null,
  loading: true,
  error: null,
  pendingFilterConfig: null,

  // Cached indices
  obsIndex: null,
  varIndex: null,

  initialize: async (url) => {
    // Reset all dataset-dependent state before loading new data
    set({
      loading: true,
      error: null,
      // Embedding slice
      selectedObsm: null,
      obsmData: null,
      obsmLoading: false,
      obsmTime: null,
      obsmStreamingData: null,
      obsmStreamingLoading: false,
      obsmStreamingTime: null,
      obsmStreamingProgress: null,
      selectedGene: null,
      geneExpression: null,
      geneLoading: false,
      colorColumn: null,
      colorData: null,
      colorLoading: false,
      tooltipColumns: [],
      tooltipData: {},
      tooltipColumnLoading: null,
      selectedPointIndices: [],
      selectionGeometry: null,
      // Columns slice
      obsColumnsSelected: [],
      obsColumnsData: {},
      obsColumnLoading: null,
      obsColumnTime: null,
      varColumnsSelected: [],
      varColumnsData: {},
      varColumnLoading: null,
      varColumnTime: null,
      // Plots slice
      plotGene: null,
      plotGeneExpression: null,
      plotGeneLoading: false,
      plotObsColumn: null,
      plotObsData: null,
      plotObsLoading: false,
      // Dotplot slice
      dotplotGenes: [],
      dotplotGeneExpressions: {},
      dotplotGeneLoading: null,
      dotplotObsColumn: null,
      dotplotObsData: null,
      dotplotObsLoading: false,
      // Config slice
      filterJson: "",
      viewConfigDefaults: {},
      appliedSelections: [],
      activeSelectionIndex: undefined,
    });

    const timings = {};

    try {
      // Phase 1: Critical path — minimum needed to render the scatterplot
      let start = performance.now();
      const [store, featureFlags] = await Promise.all([
        AnnDataStore.open(url),
        fetchFeatureFlags(),
      ]);
      timings.open = performance.now() - start;

      const obsmKeys = store.obsmKeys();
      const layerKeys = store.layerKeys();

      // Reveal the UI — scatterplot can render as soon as fetchObsm fires
      set({
        url,
        adata: store,
        metadata: { obsColumns: [], varColumns: [], obsmKeys, layerKeys, varNames: [], geneNames: [], timings, chunks: null },
        featureFlags,
        loading: false,
        error: null,
      });

      // Phase 2: Background metadata — populates sidebars, gene search, etc.
      start = performance.now();
      const [obsColumns, varColumns, varNames] = await Promise.all([
        store.obsColumns(),
        store.varColumns(),
        store.varNames(),
      ]);
      timings.obsColumns = performance.now() - start;
      timings.varColumns = timings.obsColumns;
      timings.varNames = timings.obsColumns;

      // Gene symbol lookup (depends on varNames)
      let geneNames = varNames;
      const geneSymbolCandidates = ["gene_symbol", "GeneSymbol", "feature_name", "var_name", "gene_name"];
      for (const col of geneSymbolCandidates) {
        try {
          const names = await store.varColumn(col);
          if (names && names.length === varNames.length) {
            geneNames = names;
            break;
          }
        } catch {
          // column doesn't exist, try next
        }
      }

      // X chunk info (for profiler)
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

      // Update metadata with full results
      set({
        metadata: { obsColumns, varColumns, obsmKeys, layerKeys, varNames, geneNames, timings, chunks },
        varIndex: varNames,
      });

      // Phase 3: obsNames — large fetch, deferred until the initial obsm fetch settles
      // so we don't compete for bandwidth/decompression with the scatterplot data.
      if (get().obsmLoading) {
        await new Promise((resolve) => {
          const unsub = useAppStore.subscribe((state) => {
            if (!state.obsmLoading) { unsub(); resolve(); }
          });
        });
      }
      start = performance.now();
      const obsNames = await store.obsNames();
      timings.obsNames = performance.now() - start;
      set((prev) => ({
        obsIndex: obsNames,
        metadata: { ...prev.metadata, timings: { ...prev.metadata.timings, obsNames: timings.obsNames } },
      }));

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
});

const useAppStore = create((...a) => ({
  ...createCoreSlice(...a),
  ...createColumnsSlice(...a),
  ...createEmbeddingSlice(...a),
  ...createPlotsSlice(...a),
  ...createDotplotSlice(...a),
  ...createConfigSlice(...a),
}));

export default useAppStore;
