import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router";
import useAppStore from "../../../store/useAppStore";
import ExplorerLeftSidebar from "../../layouts/ExplorerLeftSidebar";

// Mock recentUrls utility
vi.mock("../../../utils/recentUrls", () => ({
  getRecentUrls: () => [
    { url: "https://example.com/data/cells.zarr", lastLoaded: 1000 },
  ],
}));

beforeEach(() => {
  useAppStore.setState({
    url: "https://example.com/data/cells.zarr",
    isEmbedded: false,
    metadata: {
      obsmKeys: ["X_umap", "X_pca"],
      obsColumns: ["cell_type", "batch"],
      geneNames: ["TP53", "BRCA1"],
    },
    selectedObsm: "X_umap",
    obsmLoading: false,
    colorColumn: null,
    selectedGene: null,
    colorLoading: false,
    featureFlags: {},
    fetchObsm: vi.fn(),
    setColorColumn: vi.fn(),
    setSelectedGene: vi.fn(),
    clearGeneSelection: vi.fn(),
  });
});

function renderWithRouter(ui) {
  return render(<MemoryRouter>{ui}</MemoryRouter>);
}

describe("ExplorerLeftSidebar", () => {
  it("renders all three sidebar sections", () => {
    renderWithRouter(<ExplorerLeftSidebar />);
    expect(screen.getAllByText("Dataset").length).toBeGreaterThan(0);
    expect(screen.getAllByText("Embedding").length).toBeGreaterThan(0);
    expect(screen.getAllByText("Color by").length).toBeGreaterThan(0);
  });

  it("renders the ZExplorer header", () => {
    renderWithRouter(<ExplorerLeftSidebar />);
    expect(screen.getAllByText("ZExplorer").length).toBeGreaterThan(0);
  });
});
