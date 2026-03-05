import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import type { ProfileEntry } from "../types";
import MethodBreakdownChart from "./MethodBreakdownChart";
import CacheEfficiencyChart from "./CacheEfficiencyChart";
import BytesByMethodChart from "./BytesByMethodChart";
import RequestsByMethodChart from "./RequestsByMethodChart";
import BytesVsDurationChart from "./BytesVsDurationChart";
import SessionWaterfallChart from "./SessionWaterfallChart";

const sampleEntries: ProfileEntry[] = [
  { id: 1, method: "obsm", key: "X_umap", cacheHit: false, startTime: 0, duration: 50, fetches: { requests: 2, bytes: 1024 } },
  { id: 2, method: "obs", key: "cell_type", cacheHit: true, startTime: 10, duration: 20, fetches: { requests: 1, bytes: 512 } },
  { id: 3, method: "geneExpression", key: "TP53", cacheHit: false, startTime: 30, duration: 100, fetches: { requests: 4, bytes: 4096 } },
];

const sampleSessions = [
  { label: "Session 1", entries: sampleEntries },
  { label: "Session 2", entries: [sampleEntries[0], sampleEntries[2]] },
];

describe("MethodBreakdownChart", () => {
  it("renders without crashing", () => {
    const { container } = render(<MethodBreakdownChart sessions={sampleSessions} />);
    expect(container.querySelector("svg")).toBeTruthy();
  });

  it("returns null for empty sessions", () => {
    const { container } = render(<MethodBreakdownChart sessions={[]} />);
    expect(container.innerHTML).toBe("");
  });
});

describe("CacheEfficiencyChart", () => {
  it("renders without crashing", () => {
    const { container } = render(<CacheEfficiencyChart sessions={sampleSessions} />);
    expect(container.querySelector("svg")).toBeTruthy();
  });

  it("returns null for empty sessions", () => {
    const { container } = render(<CacheEfficiencyChart sessions={[]} />);
    expect(container.innerHTML).toBe("");
  });
});

describe("BytesByMethodChart", () => {
  it("renders without crashing", () => {
    const { container } = render(<BytesByMethodChart entries={sampleEntries} />);
    expect(container.querySelector("svg")).toBeTruthy();
  });

  it("returns null for empty entries", () => {
    const { container } = render(<BytesByMethodChart entries={[]} />);
    expect(container.innerHTML).toBe("");
  });
});

describe("RequestsByMethodChart", () => {
  it("renders without crashing", () => {
    const { container } = render(<RequestsByMethodChart entries={sampleEntries} />);
    expect(container.querySelector("svg")).toBeTruthy();
  });

  it("returns null for empty entries", () => {
    const { container } = render(<RequestsByMethodChart entries={[]} />);
    expect(container.innerHTML).toBe("");
  });
});

describe("BytesVsDurationChart", () => {
  it("renders without crashing", () => {
    const { container } = render(<BytesVsDurationChart entries={sampleEntries} />);
    expect(container.querySelector("svg")).toBeTruthy();
  });

  it("returns null for empty entries", () => {
    const { container } = render(<BytesVsDurationChart entries={[]} />);
    expect(container.innerHTML).toBe("");
  });
});

describe("SessionWaterfallChart", () => {
  it("renders without crashing", () => {
    const { container } = render(<SessionWaterfallChart entries={sampleEntries} />);
    expect(container.querySelector("svg")).toBeTruthy();
  });

  it("returns null for empty entries", () => {
    const { container } = render(<SessionWaterfallChart entries={[]} />);
    expect(container.innerHTML).toBe("");
  });
});
