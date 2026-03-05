import { describe, it, expect } from "vitest";
import { getMethodColor, formatBytes, formatShape, DEFAULT_METHOD_COLOR } from "./constants";

describe("getMethodColor", () => {
  it("returns the correct color for known methods", () => {
    expect(getMethodColor("obsm")).toBe("#1f77b4");
    expect(getMethodColor("obs")).toBe("#ff7f0e");
    expect(getMethodColor("X")).toBe("#8c564b");
  });

  it("returns default color for unknown methods", () => {
    expect(getMethodColor("unknown")).toBe(DEFAULT_METHOD_COLOR);
    expect(getMethodColor("")).toBe(DEFAULT_METHOD_COLOR);
  });
});

describe("formatBytes", () => {
  it("formats zero bytes", () => {
    expect(formatBytes(0)).toBe("0 B");
  });

  it("formats bytes under 1 KB", () => {
    expect(formatBytes(512)).toBe("512 B");
    expect(formatBytes(1)).toBe("1 B");
  });

  it("formats kilobytes", () => {
    expect(formatBytes(1024)).toBe("1.0 KB");
    expect(formatBytes(1536)).toBe("1.5 KB");
  });

  it("formats megabytes", () => {
    expect(formatBytes(1024 * 1024)).toBe("1.0 MB");
    expect(formatBytes(2.5 * 1024 * 1024)).toBe("2.5 MB");
  });
});

describe("formatShape", () => {
  it("returns empty string for null/undefined/empty", () => {
    expect(formatShape(null)).toBe("");
    expect(formatShape(undefined)).toBe("");
    expect(formatShape([])).toBe("");
  });

  it("formats single dimension", () => {
    expect(formatShape([100])).toBe("[100]");
  });

  it("formats multiple dimensions with multiplication sign", () => {
    expect(formatShape([100, 200])).toBe("[100 \u00d7 200]");
    expect(formatShape([10, 20, 30])).toBe("[10 \u00d7 20 \u00d7 30]");
  });
});
