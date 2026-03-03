import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import useAppStore from "../../../store/useAppStore";
import SidebarEmbeddingPicker from "../SidebarEmbeddingPicker";

// Provide minimal store state for tests
beforeEach(() => {
  useAppStore.setState({
    metadata: { obsmKeys: ["X_umap", "X_pca", "X_tsne"] },
    selectedObsm: "X_umap",
    obsmLoading: false,
    fetchObsm: vi.fn(),
  });
});

describe("SidebarEmbeddingPicker", () => {
  it("renders the section label", () => {
    render(<SidebarEmbeddingPicker />);
    expect(screen.getByText("Embedding")).toBeInTheDocument();
  });

  it("shows the selected embedding in the dropdown", () => {
    render(<SidebarEmbeddingPicker />);
    // Ant Select renders the selected value; may appear in multiple elements
    const matches = screen.getAllByTitle("X_umap");
    expect(matches.length).toBeGreaterThan(0);
  });
});
