import { describe, it, expect, vi } from "vitest";
import { ConsolidatedStore, buildV3Cache, buildV2Cache } from "./ConsolidatedStore";

function mockInnerStore(responses: Record<string, Uint8Array | undefined> = {}) {
  return {
    url: "http://example.com/test.zarr",
    get: vi.fn(async (key: string) => responses[key]),
    getRange: vi.fn(async () => new Uint8Array([1, 2, 3])),
    snapshot: vi.fn(() => ({ requests: 0, bytes: 0 })),
  };
}

const encoder = new TextEncoder();

describe("ConsolidatedStore", () => {
  it("serves cached metadata from get() without hitting inner store", async () => {
    const inner = mockInnerStore();
    const cache = new Map<string, Uint8Array>();
    cache.set("/obs/zarr.json", encoder.encode('{"node_type":"group"}'));
    const store = new ConsolidatedStore(inner as any, cache);

    const result = await store.get("/obs/zarr.json" as any);
    expect(result).toEqual(encoder.encode('{"node_type":"group"}'));
    expect(inner.get).not.toHaveBeenCalled();
  });

  it("falls through to inner store for uncached paths", async () => {
    const responseBytes = new Uint8Array([10, 20, 30]);
    const inner = mockInnerStore({ "/some/chunk": responseBytes });
    const store = new ConsolidatedStore(inner as any, new Map());

    const result = await store.get("/some/chunk" as any);
    expect(result).toEqual(responseBytes);
    expect(inner.get).toHaveBeenCalledWith("/some/chunk", undefined);
  });

  it("always delegates getRange to inner store", async () => {
    const inner = mockInnerStore();
    const cache = new Map<string, Uint8Array>();
    cache.set("/obs/zarr.json", encoder.encode("cached"));
    const store = new ConsolidatedStore(inner as any, cache);

    await store.getRange("/obs/zarr.json" as any, { suffixLength: 10 });
    expect(inner.getRange).toHaveBeenCalled();
  });

  it("exposes url from inner store", () => {
    const inner = mockInnerStore();
    const store = new ConsolidatedStore(inner as any, new Map());
    expect(store.url).toBe("http://example.com/test.zarr");
  });

  it("delegates snapshot() to inner store", () => {
    const inner = mockInnerStore();
    const store = new ConsolidatedStore(inner as any, new Map());
    store.snapshot();
    expect(inner.snapshot).toHaveBeenCalled();
  });
});

describe("buildV3Cache", () => {
  it("maps consolidated metadata keys to /key/zarr.json paths", () => {
    const rootJson = { zarr_format: 3, node_type: "group" };
    const meta = {
      obs: { zarr_format: 3, node_type: "group", attributes: {} },
      "obs/cell_type": { zarr_format: 3, node_type: "array", shape: [100] },
    };
    const cache = buildV3Cache(rootJson, meta);

    expect(cache.has("/zarr.json")).toBe(true);
    expect(cache.has("/obs/zarr.json")).toBe(true);
    expect(cache.has("/obs/cell_type/zarr.json")).toBe(true);

    const parsed = JSON.parse(new TextDecoder().decode(cache.get("/obs/zarr.json")));
    expect(parsed.node_type).toBe("group");
  });

  it("caches root zarr.json so zarrita doesn't re-fetch it", () => {
    const rootJson = { zarr_format: 3, node_type: "group", attributes: { "encoding-type": "anndata" } };
    const cache = buildV3Cache(rootJson, {});

    const parsed = JSON.parse(new TextDecoder().decode(cache.get("/zarr.json")));
    expect(parsed.attributes["encoding-type"]).toBe("anndata");
  });
});

describe("buildV2Cache", () => {
  it("maps v2 metadata keys to /key paths", () => {
    const meta = {
      ".zgroup": { zarr_format: 2 },
      ".zattrs": { "encoding-type": "anndata" },
      "obs/.zattrs": { _index: "index", "column-order": ["a"] },
      "obs/.zgroup": { zarr_format: 2 },
      "X/.zarray": { shape: [100, 50], chunks: [100, 10], dtype: "<f4" },
    };
    const cache = buildV2Cache(meta);

    expect(cache.has("/.zgroup")).toBe(true);
    expect(cache.has("/.zattrs")).toBe(true);
    expect(cache.has("/obs/.zattrs")).toBe(true);
    expect(cache.has("/X/.zarray")).toBe(true);

    const parsed = JSON.parse(new TextDecoder().decode(cache.get("/obs/.zattrs")));
    expect(parsed._index).toBe("index");
  });
});
