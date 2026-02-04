import * as zarr from "zarrita";

export class ZarrStore {
  constructor(store, root) {
    this.store = store;
    this.root = root;
    this.attrs = root.attrs;
  }

  static async open(url) {
    const store = new zarr.FetchStore(url);
    const root = await zarr.open(store, { kind: "group" });
    return new ZarrStore(store, root);
  }
}
