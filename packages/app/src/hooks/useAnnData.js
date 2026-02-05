import { useState, useEffect, useCallback } from "react";
import { AnnDataStore } from "@anndata-zarr/zarrstore";

export default function useAnnData(url) {
  const [adata, setAdata] = useState(null);
  const [metadata, setMetadata] = useState(null);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(true);

  // Cached indices
  const [obsIndex, setObsIndex] = useState(null);
  const [varIndex, setVarIndex] = useState(null);

  // Obs column state
  const [selectedObsColumn, setSelectedObsColumn] = useState(null);
  const [obsColumnData, setObsColumnData] = useState(null);
  const [obsColumnLoading, setObsColumnLoading] = useState(false);
  const [obsColumnTime, setObsColumnTime] = useState(null);

  // Var column state
  const [selectedVarColumn, setSelectedVarColumn] = useState(null);
  const [varColumnData, setVarColumnData] = useState(null);
  const [varColumnLoading, setVarColumnLoading] = useState(false);
  const [varColumnTime, setVarColumnTime] = useState(null);

  // Obsm state
  const [selectedObsm, setSelectedObsm] = useState(null);
  const [obsmData, setObsmData] = useState(null);
  const [obsmLoading, setObsmLoading] = useState(false);
  const [obsmTime, setObsmTime] = useState(null);

  // Initial load
  useEffect(() => {
    async function load() {
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

        setAdata(store);
        setMetadata({ obsColumns, varColumns, obsmKeys, layerKeys, timings, chunks });
      } catch (err) {
        setError(err.message);
        console.error(err);
      } finally {
        setLoading(false);
      }
    }
    load();
  }, [url]);

  // Fetch obs column
  const fetchObsColumn = useCallback(async (colName) => {
    if (!adata) return;
    setSelectedObsColumn(colName);
    setObsColumnLoading(true);
    setObsColumnData(null);
    try {
      const start = performance.now();
      let index = obsIndex;
      if (!index) {
        index = await adata.obsNames();
        setObsIndex(index);
      }
      const values = await adata.obsColumn(colName);
      setObsColumnTime(performance.now() - start);
      setObsColumnData({ values, index });
    } catch (err) {
      console.error(err);
      setObsColumnData({ error: err.message });
    } finally {
      setObsColumnLoading(false);
    }
  }, [adata, obsIndex]);

  // Fetch var column
  const fetchVarColumn = useCallback(async (colName) => {
    if (!adata) return;
    setSelectedVarColumn(colName);
    setVarColumnLoading(true);
    setVarColumnData(null);
    try {
      const start = performance.now();
      let index = varIndex;
      if (!index) {
        index = await adata.varNames();
        setVarIndex(index);
      }
      const values = await adata.varColumn(colName);
      setVarColumnTime(performance.now() - start);
      setVarColumnData({ values, index });
    } catch (err) {
      console.error(err);
      setVarColumnData({ error: err.message });
    } finally {
      setVarColumnLoading(false);
    }
  }, [adata, varIndex]);

  // Fetch obsm
  const fetchObsm = useCallback(async (key) => {
    if (!adata) return;
    setSelectedObsm(key);
    setObsmLoading(true);
    setObsmData(null);
    try {
      const start = performance.now();
      let index = obsIndex;
      if (!index) {
        index = await adata.obsNames();
        setObsIndex(index);
      }
      const result = await adata.obsm(key);
      setObsmTime(performance.now() - start);
      setObsmData({ ...result, index });
    } catch (err) {
      console.error(err);
      setObsmData({ error: err.message });
    } finally {
      setObsmLoading(false);
    }
  }, [adata, obsIndex]);

  // Fetch color data for scatterplot (doesn't update state, just returns values)
  const fetchColorData = useCallback(async (colName) => {
    if (!adata) return null;
    return adata.obsColumn(colName);
  }, [adata]);

  return {
    // Store
    adata,
    loading,
    error,
    metadata,

    // Obs columns
    obsColumn: {
      selected: selectedObsColumn,
      data: obsColumnData,
      loading: obsColumnLoading,
      time: obsColumnTime,
      fetch: fetchObsColumn,
    },

    // Var columns
    varColumn: {
      selected: selectedVarColumn,
      data: varColumnData,
      loading: varColumnLoading,
      time: varColumnTime,
      fetch: fetchVarColumn,
    },

    // Obsm
    obsm: {
      selected: selectedObsm,
      data: obsmData,
      loading: obsmLoading,
      time: obsmTime,
      fetch: fetchObsm,
    },

    // Utilities
    fetchColorData,
  };
}
