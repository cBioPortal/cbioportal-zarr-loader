import { describe, it, expect, afterEach } from "vitest";
import { render, screen, fireEvent, cleanup } from "@testing-library/react";
import CollapsibleLegend from "./CollapsibleLegend";

afterEach(cleanup);

const makeCategories = (n) =>
  Array.from({ length: n }, (_, i) => [`cat${i}`, [i * 10, i * 5, 200]]);

describe("CollapsibleLegend", () => {
  it("renders all categories when count is within limit", () => {
    const categories = makeCategories(5);
    render(<CollapsibleLegend categories={categories} maxHeight={400} />);
    for (let i = 0; i < 5; i++) {
      expect(screen.getByText(`cat${i}`)).toBeInTheDocument();
    }
    expect(screen.queryByText(/Show all/)).not.toBeInTheDocument();
  });

  it("renders color swatches with correct background colors", () => {
    const categories = [["TypeA", [31, 119, 180]]];
    const { container } = render(
      <CollapsibleLegend categories={categories} maxHeight={400} />
    );
    const swatch = container.querySelector('[style*="background-color"]');
    expect(swatch).toBeInTheDocument();
    expect(swatch.style.backgroundColor).toBe("rgb(31, 119, 180)");
  });

  it("shows only first 20 items when categories exceed limit", () => {
    const categories = makeCategories(25);
    render(<CollapsibleLegend categories={categories} maxHeight={400} />);
    for (let i = 0; i < 20; i++) {
      expect(screen.getByText(`cat${i}`)).toBeInTheDocument();
    }
    expect(screen.queryByText("cat20")).not.toBeInTheDocument();
    expect(screen.getByText("Show all (25)")).toBeInTheDocument();
  });

  it("expands to show all items on click", () => {
    const categories = makeCategories(25);
    render(<CollapsibleLegend categories={categories} maxHeight={400} />);
    fireEvent.click(screen.getByText("Show all (25)"));
    for (let i = 0; i < 25; i++) {
      expect(screen.getByText(`cat${i}`)).toBeInTheDocument();
    }
    expect(screen.getByText("Show less")).toBeInTheDocument();
  });

  it("collapses back to limit on second click", () => {
    const categories = makeCategories(25);
    render(<CollapsibleLegend categories={categories} maxHeight={400} />);
    fireEvent.click(screen.getByText("Show all (25)"));
    fireEvent.click(screen.getByText("Show less"));
    expect(screen.queryByText("cat20")).not.toBeInTheDocument();
    expect(screen.getByText("Show all (25)")).toBeInTheDocument();
  });
});
