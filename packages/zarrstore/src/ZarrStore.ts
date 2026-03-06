import * as zarr from "zarrita";
import type { Readable } from "zarrita";
import { InstrumentedStore } from "./InstrumentedStore";
import { ConsolidatedStore, buildV3Cache, buildV2Cache } from "./ConsolidatedStore";
import type { FetchStats } from "./InstrumentedStore";

interface ConsolidatedMetadataMap {
  [key: string]: unknown;
}

export class ZarrStore {
  store: InstrumentedStore;
  root: zarr.Group<Readable>;
  attrs: Record<string, unknown>;
  consolidatedMetadata: ConsolidatedMetadataMap | null;
  #effectiveStore: InstrumentedStore | ConsolidatedStore;

  constructor(
    store: InstrumentedStore,
    root: zarr.Group<Readable>,
    consolidatedMetadata: ConsolidatedMetadataMap | null = null,
    effectiveStore?: ConsolidatedStore,
  ) {
    this.store = store;
    this.root = root;
    this.attrs = root.attrs;
    this.consolidatedMetadata = consolidatedMetadata;
    this.#effectiveStore = effectiveStore ?? store;
  }

  static async open(url: string): Promise<ZarrStore> {
    const fetchStore = new zarr.FetchStore(url);
    const instrumented = new InstrumentedStore(fetchStore);

    let effectiveStore: ConsolidatedStore | undefined;
    let consolidatedMetadata: ConsolidatedMetadataMap | null = null;

    // Try v3: fetch root zarr.json, check for consolidated_metadata
    const rootBytes = await instrumented.get("/zarr.json" as any);
    if (rootBytes) {
      const rootJson = JSON.parse(new TextDecoder().decode(rootBytes)) as Record<string, unknown>;
      const consolidated = rootJson.consolidated_metadata as
        | { metadata?: ConsolidatedMetadataMap }
        | undefined;
      if (consolidated?.metadata) {
        consolidatedMetadata = consolidated.metadata;
        const cache = buildV3Cache(rootJson, consolidatedMetadata);
        effectiveStore = new ConsolidatedStore(instrumented, cache);
      } else {
        // v3 without consolidated — cache just the root zarr.json we already fetched
        const cache = new Map<string, Uint8Array>();
        cache.set("/zarr.json", rootBytes);
        effectiveStore = new ConsolidatedStore(instrumented, cache);
      }
    } else {
      // Try v2: fetch .zmetadata
      const zmetaBytes = await instrumented.get("/.zmetadata" as any);
      if (zmetaBytes) {
        const zmeta = JSON.parse(new TextDecoder().decode(zmetaBytes)) as Record<string, unknown>;
        if (zmeta.zarr_consolidated_format === 1 && zmeta.metadata) {
          consolidatedMetadata = zmeta.metadata as ConsolidatedMetadataMap;
          const cache = buildV2Cache(consolidatedMetadata);
          effectiveStore = new ConsolidatedStore(instrumented, cache);
        }
      }
    }

    const root = await zarr.open(effectiveStore ?? instrumented, { kind: "group" });
    return new ZarrStore(instrumented, root, consolidatedMetadata, effectiveStore);
  }

  async openArray(path: string): Promise<zarr.Array<zarr.DataType, Readable>> {
    return zarr.open(this.root.resolve(path), { kind: "array" });
  }

  async openGroup(path: string): Promise<zarr.Group<Readable>> {
    return zarr.open(this.root.resolve(path), { kind: "group" });
  }

  snapshotFetchStats(): FetchStats {
    return this.#effectiveStore.snapshot();
  }
}
