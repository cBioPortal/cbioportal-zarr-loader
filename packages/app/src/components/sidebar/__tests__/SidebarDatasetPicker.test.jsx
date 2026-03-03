import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router";
import useAppStore from "../../../store/useAppStore";
import SidebarDatasetPicker from "../SidebarDatasetPicker";

// Mock recentUrls utility
vi.mock("../../../utils/recentUrls", () => ({
  getRecentUrls: () => [
    { url: "https://example.com/data/cells.zarr", lastLoaded: 1000 },
    { url: "https://example.com/other/neurons.zarr", lastLoaded: 900 },
  ],
}));

beforeEach(() => {
  useAppStore.setState({
    url: "https://example.com/data/cells.zarr",
    isEmbedded: false,
  });
});

function renderWithRouter(ui) {
  return render(<MemoryRouter>{ui}</MemoryRouter>);
}

describe("SidebarDatasetPicker", () => {
  it("renders the section label", () => {
    renderWithRouter(<SidebarDatasetPicker />);
    expect(screen.getAllByText("Dataset").length).toBeGreaterThan(0);
  });

  it("shows filename of current URL in the select", () => {
    renderWithRouter(<SidebarDatasetPicker />);
    // The select displays the filename label
    expect(screen.getAllByText("cells.zarr").length).toBeGreaterThan(0);
  });

  it("renders nothing when embedded", () => {
    useAppStore.setState({ isEmbedded: true });
    const { container } = renderWithRouter(<SidebarDatasetPicker />);
    expect(container.querySelector("[style]")).toBeNull();
  });

  it("renders the paste URL input", () => {
    renderWithRouter(<SidebarDatasetPicker />);
    const inputs = screen.getAllByPlaceholderText("Paste URL...");
    expect(inputs.length).toBeGreaterThan(0);
  });
});
