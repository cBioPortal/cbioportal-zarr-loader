import * as zarr from "zarrita";
import { ZarrStore } from "./ZarrStore.js";
import {
  readArray,
  decodeDataframe,
  decodeColumn,
  decodeCategorical,
  decodeSparseMatrix,
  decodeNode,
  toStringArray,
} from "./decoders.js";

export class AnnDataStore {
  #zarrStore;
  #shape;
  #attrs;
  #consolidatedMetadata;
  #cache = new Map();

  constructor(zarrStore, shape, consolidatedMetadata = null) {
    this.#zarrStore = zarrStore;
    this.#attrs = zarrStore.attrs;
    this.#shape = shape;
    this.#consolidatedMetadata = consolidatedMetadata;
  }

  #cached(key, fn) {
    if (!this.#cache.has(key)) {
      this.#cache.set(key, fn());
    }
    return this.#cache.get(key);
  }

  clearCache() {
    this.#cache.clear();
  }

  static async open(url) {
    const zarrStore = await ZarrStore.open(url);
    const attrs = zarrStore.attrs;

    if (attrs["encoding-type"] !== "anndata") {
      throw new Error(
        `Expected encoding-type "anndata", got "${attrs["encoding-type"]}"`,
      );
    }

    const shape = await AnnDataStore.#resolveShape(zarrStore);
    const consolidatedMetadata = await AnnDataStore.#loadConsolidatedMetadata(
      zarrStore,
    );
    return new AnnDataStore(zarrStore, shape, consolidatedMetadata);
  }

  static async fromZarrStore(zarrStore) {
    const attrs = zarrStore.attrs;

    if (attrs["encoding-type"] !== "anndata") {
      throw new Error(
        `Expected encoding-type "anndata", got "${attrs["encoding-type"]}"`,
      );
    }

    const shape = await AnnDataStore.#resolveShape(zarrStore);
    const consolidatedMetadata = await AnnDataStore.#loadConsolidatedMetadata(
      zarrStore,
    );
    return new AnnDataStore(zarrStore, shape, consolidatedMetadata);
  }

  static async #loadConsolidatedMetadata(zarrStore) {
    try {
      const response = await fetch(
        zarrStore.store.url.replace(/\/$/, "") + "/.zmetadata",
      );
      if (!response.ok) return null;
      const data = await response.json();
      return data.metadata || null;
    } catch {
      return null;
    }
  }

  static async #resolveShape(zarrStore) {
    // Try opening X as an array first (dense), fall back to group (sparse)
    try {
      const xArr = await zarrStore.openArray("X");
      return xArr.shape;
    } catch {
      const xGroup = await zarrStore.openGroup("X");
      return xGroup.attrs.shape;
    }
  }

  // --- Metadata (synchronous) ---

  get shape() {
    return this.#shape;
  }

  get nObs() {
    return this.#shape[0];
  }

  get nVar() {
    return this.#shape[1];
  }

  get attrs() {
    return this.#attrs;
  }

  get zarrStore() {
    return this.#zarrStore;
  }

  // --- X matrix ---

  async X(sliceRange) {
    let node;
    try {
      node = await this.#zarrStore.openArray("X");
    } catch {
      node = await this.#zarrStore.openGroup("X");
    }

    if (node.attrs?.["encoding-type"]?.endsWith("_matrix")) {
      return decodeSparseMatrix(node);
    }

    // Dense array
    if (sliceRange) {
      const [start, end] = sliceRange;
      const chunk = await zarr.get(node, [zarr.slice(start, end), null]);
      return { data: chunk.data, shape: chunk.shape };
    }
    return readArray(node);
  }

  async geneExpression(geneName) {
    // Get gene index from var names
    const varNames = await this.varNames();
    const geneIndex = varNames.indexOf(geneName);
    if (geneIndex === -1) {
      throw new Error(`Gene "${geneName}" not found`);
    }

    // Try to open X as dense array
    let node;
    try {
      node = await this.#zarrStore.openArray("X");
    } catch {
      node = await this.#zarrStore.openGroup("X");
    }

    if (node.attrs?.["encoding-type"]?.endsWith("_matrix")) {
      // Sparse matrix - need to decode and extract column
      const sparse = await decodeSparseMatrix(node);
      // For CSR matrix, we need to iterate through all rows
      // This is less efficient but works for any sparse format
      const result = new Float32Array(this.#shape[0]);
      const { data, indices, indptr } = sparse;
      for (let row = 0; row < this.#shape[0]; row++) {
        const rowStart = indptr[row];
        const rowEnd = indptr[row + 1];
        for (let j = rowStart; j < rowEnd; j++) {
          if (indices[j] === geneIndex) {
            result[row] = data[j];
            break;
          }
        }
      }
      return result;
    }

    // Dense array - slice the column
    const chunk = await zarr.get(node, [null, geneIndex]);
    return chunk.data;
  }

  // --- obs / var dataframes ---

  obs() {
    return this.#cached("obs", async () => {
      const group = await this.#zarrStore.openGroup("obs");
      return decodeDataframe(group);
    });
  }

  obsColumn(name) {
    return this.#cached(`obs:${name}`, async () => {
      const group = await this.#zarrStore.openGroup("obs");
      return decodeColumn(group, name);
    });
  }

  obsColumns() {
    return this.#cached("obsColumns", async () => {
      const group = await this.#zarrStore.openGroup("obs");
      return Array.from(group.attrs["column-order"]);
    });
  }

  obsNames() {
    return this.#cached("obsNames", async () => {
      const group = await this.#zarrStore.openGroup("obs");
      const indexKey = group.attrs["_index"];
      // Index can be an array or a categorical group
      try {
        const arr = await zarr.open(group.resolve(indexKey), { kind: "array" });
        const result = await readArray(arr);
        return toStringArray(result.data);
      } catch {
        // It's a categorical group
        const catGroup = await zarr.open(group.resolve(indexKey), {
          kind: "group",
        });
        const decoded = await decodeCategorical(catGroup);
        return decoded.values;
      }
    });
  }

  var() {
    return this.#cached("var", async () => {
      const group = await this.#zarrStore.openGroup("var");
      return decodeDataframe(group);
    });
  }

  varColumn(name) {
    return this.#cached(`var:${name}`, async () => {
      const group = await this.#zarrStore.openGroup("var");
      return decodeColumn(group, name);
    });
  }

  varColumns() {
    return this.#cached("varColumns", async () => {
      const group = await this.#zarrStore.openGroup("var");
      return Array.from(group.attrs["column-order"]);
    });
  }

  varNames() {
    return this.#cached("varNames", async () => {
      const group = await this.#zarrStore.openGroup("var");
      const indexKey = group.attrs["_index"];
      // Index can be an array or a categorical group
      try {
        const arr = await zarr.open(group.resolve(indexKey), { kind: "array" });
        const result = await readArray(arr);
        return toStringArray(result.data);
      } catch {
        // It's a categorical group
        const catGroup = await zarr.open(group.resolve(indexKey), {
          kind: "group",
        });
        const decoded = await decodeCategorical(catGroup);
        return decoded.values;
      }
    });
  }

  // --- Dict-of-matrices slots ---

  #slotKeys(path) {
    if (!this.#consolidatedMetadata) {
      return [];
    }
    const prefix = path + "/";
    const keys = new Set();
    for (const key of Object.keys(this.#consolidatedMetadata)) {
      if (key.startsWith(prefix)) {
        const rest = key.slice(prefix.length);
        const slashIndex = rest.indexOf("/");
        if (slashIndex > 0) {
          keys.add(rest.slice(0, slashIndex));
        }
      }
    }
    return Array.from(keys);
  }

  #slotNode(path, key) {
    return this.#cached(`${path}:${key}`, async () => {
      const node = await this.#zarrStore.openGroup(`${path}/${key}`);
      return decodeNode(node);
    });
  }

  #slotArray(path, key) {
    return this.#cached(`${path}:${key}`, async () => {
      try {
        const arr = await this.#zarrStore.openArray(`${path}/${key}`);
        return readArray(arr);
      } catch {
        const node = await this.#zarrStore.openGroup(`${path}/${key}`);
        return decodeNode(node);
      }
    });
  }

  obsm(key) {
    return this.#slotArray("obsm", key);
  }

  obsmKeys() {
    return this.#slotKeys("obsm");
  }

  varm(key) {
    return this.#slotArray("varm", key);
  }

  varmKeys() {
    return this.#slotKeys("varm");
  }

  obsp(key) {
    return this.#slotNode("obsp", key);
  }

  obspKeys() {
    return this.#slotKeys("obsp");
  }

  varp(key) {
    return this.#slotNode("varp", key);
  }

  varpKeys() {
    return this.#slotKeys("varp");
  }

  // --- Layers ---

  layer(key) {
    return this.#slotArray("layers", key);
  }

  layerKeys() {
    return this.#slotKeys("layers");
  }

  // --- Unstructured (uns) ---

  uns(key) {
    return this.#slotArray("uns", key);
  }

  unsKeys() {
    return this.#slotKeys("uns");
  }
}
