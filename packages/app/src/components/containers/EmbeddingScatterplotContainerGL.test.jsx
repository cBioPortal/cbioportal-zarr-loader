import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, cleanup } from "@testing-library/react";

// Mock the GL chart to capture props
let capturedProps = {};
vi.mock("../charts/EmbeddingScatterplotGL", () => ({
  default: (props) => {
    capturedProps = props;
    return <div data-testid="scatterplot-gl">GL: {props.label}</div>;
  },
}));

// Mock the store
let storeState = {};
vi.mock("../../store/useAppStore", () => {
  const store = vi.fn(() => storeState);
  store.getState = vi.fn(() => storeState);
  store.setState = vi.fn();
  return { default: store };
});

const { default: EmbeddingScatterplotContainerGL } = await import(
  "./EmbeddingScatterplotContainerGL"
);

const mockSetSelectedPoints = vi.fn();
const mockClearSelectedPoints = vi.fn();
const mockSetSelectionGeometry = vi.fn();

const baseStoreState = {
  selectedPointIndices: [],
  setSelectedPoints: mockSetSelectedPoints,
  clearSelectedPoints: mockClearSelectedPoints,
  selectionGeometry: null,
  setSelectionGeometry: mockSetSelectionGeometry,
};

describe("EmbeddingScatterplotContainerGL", () => {
  beforeEach(() => {
    storeState = { ...baseStoreState };
    capturedProps = {};
  });
  afterEach(cleanup);

  it("passes Float32Array data directly to chart", () => {
    const data = new Float32Array([1, 2, 3, 4, 5, 6]);
    render(
      <EmbeddingScatterplotContainerGL data={data} shape={[3, 2]} label="X_umap" />,
    );
    expect(screen.getByTestId("scatterplot-gl")).toBeInTheDocument();
    expect(capturedProps.data).toBe(data); // same reference — no copy
    expect(capturedProps.shape).toEqual([3, 2]);
    expect(capturedProps.label).toBe("X_umap");
  });

  it("converts Float64Array to Float32Array", () => {
    const data = new Float64Array([1.5, 2.5]);
    render(
      <EmbeddingScatterplotContainerGL data={data} shape={[1, 2]} label="X_umap" />,
    );
    expect(capturedProps.data).toBeInstanceOf(Float32Array);
  });

  it("passes selection state from store", () => {
    storeState = { ...baseStoreState, selectedPointIndices: [0, 2] };
    const data = new Float32Array([1, 2, 3, 4, 5, 6]);
    render(
      <EmbeddingScatterplotContainerGL data={data} shape={[3, 2]} label="X_umap" />,
    );
    expect(capturedProps.selectedPointIndices).toEqual([0, 2]);
    expect(capturedProps.setSelectedPoints).toBe(mockSetSelectedPoints);
    expect(capturedProps.clearSelectedPoints).toBe(mockClearSelectedPoints);
  });
});
