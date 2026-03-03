import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen, cleanup } from "@testing-library/react";

// Mock deck.gl
vi.mock("@deck.gl/react", () => ({
  default: vi.fn((props) => <div data-testid="deckgl" />),
}));
vi.mock("@deck.gl/layers", () => ({
  ScatterplotLayer: vi.fn(),
}));
vi.mock("@deck.gl/core", () => ({
  OrthographicView: vi.fn(),
}));

// Mock child components
vi.mock("../ui/SelectionOverlay", () => ({
  default: () => <div data-testid="selection-overlay" />,
}));
vi.mock("../../hooks/useSelectionInteraction", () => ({
  default: () => ({
    selectMode: "pan",
    setSelectMode: vi.fn(),
    selectionRectRef: { current: null },
    lassoSvgRef: { current: null },
    handleMouseDown: vi.fn(),
    handleMouseMove: vi.fn(),
    handleMouseUp: vi.fn(),
  }),
}));

// Override ResizeObserver with a trackable mock
const mockObserve = vi.fn();
const mockDisconnect = vi.fn();
globalThis.ResizeObserver = class {
  observe(...args) { mockObserve(...args); }
  disconnect(...args) { mockDisconnect(...args); }
  unobserve() {}
};

const { default: EmbeddingScatterplotGL } = await import(
  "./EmbeddingScatterplotGL"
);

const defaultProps = {
  data: new Float32Array([0, 0, 1, 1, 2, 2]),
  shape: [3, 2],
  label: "X_umap",
  selectedPointIndices: [],
  setSelectedPoints: vi.fn(),
  clearSelectedPoints: vi.fn(),
  selectionGeometry: null,
  setSelectionGeometry: vi.fn(),
};

describe("EmbeddingScatterplotGL", () => {
  afterEach(cleanup);

  it("renders DeckGL", () => {
    render(<EmbeddingScatterplotGL {...defaultProps} />);
    expect(screen.getByTestId("deckgl")).toBeInTheDocument();
  });

  it("renders axis labels from the label prop", () => {
    render(<EmbeddingScatterplotGL {...defaultProps} />);
    expect(screen.getByText("X_umap_1")).toBeInTheDocument();
    expect(screen.getByText("X_umap_2")).toBeInTheDocument();
  });

  it("renders rectangle and lasso select buttons", () => {
    render(<EmbeddingScatterplotGL {...defaultProps} />);
    expect(screen.getByTitle("Rectangle select")).toBeInTheDocument();
    expect(screen.getByTitle("Lasso select")).toBeInTheDocument();
  });

  it("shows selection count badge when points are selected", () => {
    render(
      <EmbeddingScatterplotGL {...defaultProps} selectedPointIndices={[0, 1]} />,
    );
    expect(screen.getByText("2 selected")).toBeInTheDocument();
  });

  it("does not show selection count when nothing is selected", () => {
    render(<EmbeddingScatterplotGL {...defaultProps} />);
    expect(screen.queryByText(/selected/)).not.toBeInTheDocument();
  });

  it("shows clear selection button when points are selected", () => {
    render(
      <EmbeddingScatterplotGL {...defaultProps} selectedPointIndices={[0]} />,
    );
    expect(screen.getByTitle("Clear selection")).toBeInTheDocument();
  });

  it("does not render hexbin, legend, tooltip, or color elements", () => {
    render(<EmbeddingScatterplotGL {...defaultProps} />);
    expect(screen.queryByTitle("Switch to hexbin")).not.toBeInTheDocument();
    expect(screen.queryByTestId("expression-legend")).not.toBeInTheDocument();
    expect(screen.queryByTestId("collapsible-legend")).not.toBeInTheDocument();
    expect(screen.queryByTestId("hover-tooltip")).not.toBeInTheDocument();
    expect(screen.queryByTitle("Plot settings")).not.toBeInTheDocument();
  });

  it("renders selection overlay", () => {
    render(<EmbeddingScatterplotGL {...defaultProps} />);
    expect(screen.getByTestId("selection-overlay")).toBeInTheDocument();
  });

  it("attaches ResizeObserver to fill available width", () => {
    render(<EmbeddingScatterplotGL {...defaultProps} />);
    expect(mockObserve).toHaveBeenCalled();
  });
});
