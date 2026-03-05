import { describe, it, expect, vi, beforeAll, afterEach } from "vitest";
import { render, screen, cleanup } from "@testing-library/react";
import ProfileBar from "./ProfileBar";
import type { ProfileEntry } from "../types";

beforeAll(() => {
  globalThis.ResizeObserver = class {
    observe() {}
    unobserve() {}
    disconnect() {}
  } as unknown as typeof ResizeObserver;
});

afterEach(cleanup);

const sampleEntries: ProfileEntry[] = [
  { id: 1, method: "obsm", key: "X_umap", cacheHit: false, startTime: 0, duration: 50, fetches: { requests: 2, bytes: 1024 } },
  { id: 2, method: "obs", key: "cell_type", cacheHit: true, startTime: 10, duration: 20 },
];

function createMockProfiler(entries: ProfileEntry[] = sampleEntries) {
  let subscriber: (() => void) | null = null;
  return {
    entries,
    version: 1,
    subscribe: (cb: () => void) => {
      subscriber = cb;
      return () => { subscriber = null; };
    },
    clear: vi.fn(),
    toJSON: () => entries,
    _notify: () => subscriber?.(),
  };
}

describe("ProfileBar", () => {
  it("renders without crashing when no profiler provided", () => {
    render(<ProfileBar />);
    expect(screen.getByText("Query Profiler")).toBeTruthy();
  });

  it("displays entry count and stats", () => {
    const profiler = createMockProfiler();
    render(<ProfileBar profiler={profiler} />);
    expect(screen.getByText("Queries: 2")).toBeTruthy();
  });

  it("calls onSave when Save button is clicked", () => {
    const profiler = createMockProfiler();
    const onSave = vi.fn();
    const { container } = render(<ProfileBar profiler={profiler} onSave={onSave} />);

    // Find the actual <button> element containing "Save"
    const buttons = container.querySelectorAll("button");
    const saveBtn = Array.from(buttons).find((b) => b.textContent?.includes("Save"));
    expect(saveBtn).toBeTruthy();

    // Simulate click with stopPropagation to avoid parent handler
    const event = new MouseEvent("click", { bubbles: true, cancelable: true });
    saveBtn!.dispatchEvent(event);

    expect(onSave).toHaveBeenCalledWith(sampleEntries);
  });

  it("calls profiler.clear when Clear button is clicked", () => {
    const profiler = createMockProfiler();
    const { container } = render(<ProfileBar profiler={profiler} />);

    const buttons = container.querySelectorAll("button");
    const clearBtn = Array.from(buttons).find((b) => b.textContent?.includes("Clear"));
    expect(clearBtn).toBeTruthy();

    const event = new MouseEvent("click", { bubbles: true, cancelable: true });
    clearBtn!.dispatchEvent(event);

    expect(profiler.clear).toHaveBeenCalledOnce();
  });

  it("renders renderLink content when provided", () => {
    const profiler = createMockProfiler();
    render(
      <ProfileBar
        profiler={profiler}
        renderLink={(children) => <a href="/profile" data-testid="profile-link">{children}</a>}
      />
    );

    expect(screen.getByTestId("profile-link")).toBeTruthy();
  });

  it("does not render link when renderLink is not provided", () => {
    const profiler = createMockProfiler();
    render(<ProfileBar profiler={profiler} />);
    expect(screen.queryByTestId("profile-link")).toBeNull();
  });

  it("disables Save and Clear when no entries", () => {
    const profiler = createMockProfiler([]);
    const { container } = render(<ProfileBar profiler={profiler} />);

    const buttons = container.querySelectorAll("button");
    const saveBtn = Array.from(buttons).find((b) => b.textContent?.includes("Save"));
    const clearBtn = Array.from(buttons).find((b) => b.textContent?.includes("Clear"));

    expect(saveBtn?.disabled).toBe(true);
    expect(clearBtn?.disabled).toBe(true);
  });
});
