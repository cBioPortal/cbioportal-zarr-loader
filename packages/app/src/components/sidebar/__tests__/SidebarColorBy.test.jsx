import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import useAppStore from "../../../store/useAppStore";
import SidebarColorBy from "../SidebarColorBy";

beforeEach(() => {
  useAppStore.setState({
    metadata: {
      obsColumns: ["cell_type", "batch", "sample_id"],
      geneNames: ["TP53", "BRCA1", "EGFR"],
    },
    colorColumn: null,
    selectedGene: null,
    colorLoading: false,
    geneLoading: false,
    setColorColumn: vi.fn(),
    setSelectedGene: vi.fn(),
    clearGeneSelection: vi.fn(),
  });
});

describe("SidebarColorBy", () => {
  it("renders the section label and mode selector", () => {
    render(<SidebarColorBy />);
    expect(screen.getAllByText("Color by").length).toBeGreaterThan(0);
  });

  it("renders the searchable select with placeholder", () => {
    render(<SidebarColorBy />);
    expect(screen.getAllByText("Select column...").length).toBeGreaterThan(0);
  });

  it("shows the selected column when set", () => {
    useAppStore.setState({ colorColumn: "cell_type" });
    render(<SidebarColorBy />);
    expect(screen.getAllByText("cell_type").length).toBeGreaterThan(0);
  });
});
