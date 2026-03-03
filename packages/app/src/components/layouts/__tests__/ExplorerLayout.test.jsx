import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, cleanup } from "@testing-library/react";
import ExplorerLayout from "../ExplorerLayout";

// Mock the child components
vi.mock("../ExplorerLeftSidebar", () => ({
  default: () => <div data-testid="left-sidebar">Left</div>,
}));
vi.mock("../ExplorerRightSidebar", () => ({
  default: () => <div data-testid="right-sidebar">Right</div>,
}));
vi.mock("../../containers/EmbeddingScatterplotContainerGL", () => ({
  default: (props) => <div data-testid="scatterplot">Scatterplot: {props.label}</div>,
}));

// Mock the store
vi.mock("../../../store/useAppStore", () => {
  const store = vi.fn(() => ({
    metadata: { obsmKeys: ["X_umap"] },
    selectedObsm: "X_umap",
    obsmData: { data: new Float32Array([1, 2, 3, 4]), shape: [2, 2] },
    obsmLoading: false,
    fetchObsm: vi.fn(),
    featureFlags: {},
  }));
  store.getState = vi.fn(() => store());
  store.setState = vi.fn();
  return { default: store };
});

// Mock useMediaQuery
let mockIsDesktop = true;
vi.mock("../../../hooks/useMediaQuery", () => ({
  default: () => mockIsDesktop,
}));

describe("ExplorerLayout", () => {
  beforeEach(() => {
    mockIsDesktop = true;
  });

  afterEach(() => {
    cleanup();
  });

  it("renders three columns on desktop", () => {
    render(<ExplorerLayout />);
    expect(screen.getByTestId("left-sidebar")).toBeInTheDocument();
    expect(screen.getByTestId("scatterplot")).toBeInTheDocument();
    expect(screen.getByTestId("right-sidebar")).toBeInTheDocument();
  });

  it("shows drawer toggle buttons on mobile", () => {
    mockIsDesktop = false;
    render(<ExplorerLayout />);
    // Scatterplot should still render
    expect(screen.getByTestId("scatterplot")).toBeInTheDocument();
    // Drawer toggle buttons should exist
    expect(screen.getByLabelText("Open left sidebar")).toBeInTheDocument();
    expect(screen.getByLabelText("Open right sidebar")).toBeInTheDocument();
  });

  it("opens left drawer when left button is clicked on mobile", () => {
    mockIsDesktop = false;
    render(<ExplorerLayout />);

    // Click the left sidebar button
    fireEvent.click(screen.getByLabelText("Open left sidebar"));

    // After clicking, left sidebar content should appear in the drawer
    expect(screen.getByTestId("left-sidebar")).toBeInTheDocument();
  });
});
