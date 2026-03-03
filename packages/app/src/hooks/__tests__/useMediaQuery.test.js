import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, act } from "@testing-library/react";
import useMediaQuery from "../useMediaQuery";

describe("useMediaQuery", () => {
  let listeners;
  let mockMatchMedia;

  beforeEach(() => {
    listeners = [];
    mockMatchMedia = vi.fn((query) => {
      const mql = {
        matches: false,
        media: query,
        addEventListener: vi.fn((event, cb) => listeners.push(cb)),
        removeEventListener: vi.fn((event, cb) => {
          listeners = listeners.filter((l) => l !== cb);
        }),
      };
      return mql;
    });
    window.matchMedia = mockMatchMedia;
  });

  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("returns false when media query does not match", () => {
    const { result } = renderHook(() => useMediaQuery("(min-width: 768px)"));
    expect(result.current).toBe(false);
  });

  it("returns true when media query matches initially", () => {
    mockMatchMedia.mockImplementation((query) => ({
      matches: true,
      media: query,
      addEventListener: vi.fn(),
      removeEventListener: vi.fn(),
    }));
    const { result } = renderHook(() => useMediaQuery("(min-width: 768px)"));
    expect(result.current).toBe(true);
  });

  it("updates when the media query match changes", () => {
    const { result } = renderHook(() => useMediaQuery("(min-width: 768px)"));
    expect(result.current).toBe(false);

    act(() => {
      listeners.forEach((cb) => cb({ matches: true }));
    });
    expect(result.current).toBe(true);
  });

  it("cleans up listener on unmount", () => {
    const { unmount } = renderHook(() => useMediaQuery("(min-width: 768px)"));
    const mql = mockMatchMedia.mock.results[0].value;
    unmount();
    expect(mql.removeEventListener).toHaveBeenCalledWith("change", expect.any(Function));
  });
});
