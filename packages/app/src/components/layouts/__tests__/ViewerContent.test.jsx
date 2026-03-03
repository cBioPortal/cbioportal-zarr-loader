import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router";

// Must mock before importing ViewerContent
vi.mock("../../../store/useAppStore", () => {
  const store = vi.fn(() => ({
    featureFlags: {},
    isEmbedded: false,
    metadata: { obsmKeys: ["X_umap"] },
    selectedObsm: "X_umap",
    obsmData: { data: new Float32Array([1, 2, 3, 4]), shape: [2, 2] },
    obsmLoading: false,
    fetchObsm: vi.fn(),
  }));
  store.getState = vi.fn(() => store());
  store.setState = vi.fn();
  return { default: store };
});

vi.mock("../ExplorerLayout", () => ({
  default: () => <div data-testid="explorer-layout">ExplorerLayout</div>,
}));

vi.mock("../../views/ObsmTab", () => ({ default: () => <div>ObsmTab</div> }));
vi.mock("../../views/ColumnsTab", () => ({ default: () => <div>ColumnsTab</div> }));
vi.mock("../../views/PlotsTab", () => ({ default: () => <div>PlotsTab</div> }));
vi.mock("../../views/DotplotTab", () => ({ default: () => <div>DotplotTab</div> }));
vi.mock("../../views/InfoTab", () => ({ default: () => <div>InfoTab</div> }));

// Import after mocks
const { ViewerContent } = await import("../../../App.jsx");

describe("ViewerContent", () => {
  it("renders ViewerTabs by default (no layout param)", () => {
    render(
      <MemoryRouter initialEntries={["/"]}>
        <ViewerContent />
      </MemoryRouter>,
    );
    expect(screen.queryByTestId("explorer-layout")).not.toBeInTheDocument();
    // Tabs should render (look for tab text)
    expect(screen.getByText("Explore")).toBeInTheDocument();
  });

  it("renders ExplorerLayout when ?layout=v3 is set", () => {
    render(
      <MemoryRouter initialEntries={["/?layout=v3"]}>
        <ViewerContent />
      </MemoryRouter>,
    );
    expect(screen.getByTestId("explorer-layout")).toBeInTheDocument();
  });
});
