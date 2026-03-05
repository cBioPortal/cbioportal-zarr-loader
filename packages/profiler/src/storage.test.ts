import { describe, it, expect, beforeEach, vi } from "vitest";
import {
  getProfileHistory,
  saveProfileSession,
  exportProfileHistory,
  importProfileHistory,
} from "./storage";
import type { ProfileSession } from "./types";

function makeSession(url: string, timestamp: number, entries: unknown[] = []): ProfileSession {
  return { url, timestamp, nObs: 100, nVar: 50, entries: entries as ProfileSession["entries"] };
}

describe("profileStorage", () => {
  beforeEach(() => {
    localStorage.clear();
  });

  describe("exportProfileHistory", () => {
    it("triggers a download with correct filename and JSON content", () => {
      saveProfileSession("https://a.zarr", 10, 5, [{ method: "get", duration: 1 }] as any);

      const clicked = vi.fn();
      const created: Record<string, string> = {};
      vi.spyOn(document, "createElement").mockReturnValueOnce({
        set href(v: string) { created.href = v; },
        set download(v: string) { created.download = v; },
        click: clicked,
      } as unknown as HTMLElement);
      const appendChild = vi.spyOn(document.body, "appendChild").mockImplementation(() => null as any);
      const removeChild = vi.spyOn(document.body, "removeChild").mockImplementation(() => null as any);
      const revokeURL = vi.spyOn(URL, "revokeObjectURL").mockImplementation(() => {});

      exportProfileHistory();

      expect(clicked).toHaveBeenCalledOnce();
      expect(created.download).toMatch(/^czl-profile-\d{4}-\d{2}-\d{2}\.json$/);
      expect(appendChild).toHaveBeenCalledOnce();
      expect(removeChild).toHaveBeenCalledOnce();
      expect(revokeURL).toHaveBeenCalledOnce();
    });
  });

  describe("importProfileHistory", () => {
    function fileFromJSON(data: unknown) {
      const json = JSON.stringify(data);
      return new File([json], "test.json", { type: "application/json" });
    }

    it("imports sessions into empty history", async () => {
      const sessions = [makeSession("https://a.zarr", 1), makeSession("https://b.zarr", 2)];
      const result = await importProfileHistory(fileFromJSON(sessions));

      expect(result).toEqual({ success: true, count: 2 });
      expect(getProfileHistory()).toHaveLength(2);
    });

    it("deduplicates by timestamp + url", async () => {
      saveProfileSession("https://a.zarr", 10, 5, []);
      const existing = getProfileHistory();

      const file = fileFromJSON([
        makeSession("https://a.zarr", existing[0].timestamp),
        makeSession("https://b.zarr", 999),
      ]);
      const result = await importProfileHistory(file);

      expect(result).toEqual({ success: true, count: 1 });
      expect(getProfileHistory()).toHaveLength(2);
    });

    it("places new sessions before existing ones", async () => {
      saveProfileSession("https://existing.zarr", 10, 5, []);
      const file = fileFromJSON([makeSession("https://new.zarr", 999)]);
      await importProfileHistory(file);

      const history = getProfileHistory();
      expect(history[0].url).toBe("https://new.zarr");
      expect(history[1].url).toBe("https://existing.zarr");
    });

    it("caps at 50 sessions", async () => {
      const sessions = Array.from({ length: 55 }, (_, i) =>
        makeSession(`https://${i}.zarr`, i)
      );
      const result = await importProfileHistory(fileFromJSON(sessions));

      expect(result.success).toBe(true);
      expect(getProfileHistory()).toHaveLength(50);
    });

    it("returns error for non-array JSON", async () => {
      const result = await importProfileHistory(fileFromJSON({ not: "array" }));
      expect(result).toEqual({ success: false, error: "File does not contain a JSON array" });
    });

    it("returns error for invalid JSON", async () => {
      const file = new File(["not json{{{"], "bad.json", { type: "application/json" });
      const result = await importProfileHistory(file);
      expect(result).toEqual({ success: false, error: "Invalid JSON file" });
    });

    it("importing same file twice produces no duplicates", async () => {
      const sessions = [makeSession("https://a.zarr", 100), makeSession("https://b.zarr", 200)];
      await importProfileHistory(fileFromJSON(sessions));
      const result = await importProfileHistory(fileFromJSON(sessions));

      expect(result).toEqual({ success: true, count: 0 });
      expect(getProfileHistory()).toHaveLength(2);
    });
  });
});
