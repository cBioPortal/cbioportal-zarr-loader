import { describe, it, expect, afterEach } from "vitest";
import { render, screen, cleanup } from "@testing-library/react";
import HoverTooltip from "./HoverTooltip";

const baseHoverInfo = {
  x: 100,
  y: 200,
  object: {
    position: [1.23456789, -0.98765432],
    index: 0,
    category: "ClusterA",
    expression: 2.71828,
  },
};

const baseProps = {
  hoverInfo: baseHoverInfo,
  colorColumn: "celltype",
  selectedGene: "TP53",
  tooltipData: {},
  hasColorData: false,
  hasGeneExpression: false,
};

afterEach(cleanup);

describe("HoverTooltip", () => {
  it("renders x/y coordinates formatted to 4 decimals", () => {
    render(<HoverTooltip {...baseProps} />);
    expect(screen.getByText("x: 1.2346")).toBeInTheDocument();
    expect(screen.getByText("y: -0.9877")).toBeInTheDocument();
  });

  it("shows category when hasColorData is true", () => {
    render(<HoverTooltip {...baseProps} hasColorData={true} />);
    expect(screen.getByText("celltype: ClusterA")).toBeInTheDocument();
  });

  it("hides category when hasColorData is false", () => {
    render(<HoverTooltip {...baseProps} hasColorData={false} />);
    expect(screen.queryByText("celltype: ClusterA")).not.toBeInTheDocument();
  });

  it("shows gene expression when hasGeneExpression is true", () => {
    render(<HoverTooltip {...baseProps} hasGeneExpression={true} />);
    expect(screen.getByText("TP53: 2.7183")).toBeInTheDocument();
  });

  it("hides gene expression when hasGeneExpression is false", () => {
    render(<HoverTooltip {...baseProps} hasGeneExpression={false} />);
    expect(screen.queryByText(/TP53:/)).not.toBeInTheDocument();
  });

  it("renders all tooltip columns from tooltipData", () => {
    const tooltipData = {
      tissue: ["brain", "liver", "lung"],
      stage: ["I", "II", "III"],
    };
    render(
      <HoverTooltip
        {...baseProps}
        tooltipData={tooltipData}
        hoverInfo={{
          ...baseHoverInfo,
          object: { ...baseHoverInfo.object, index: 1 },
        }}
      />
    );
    expect(screen.getByText("tissue: liver")).toBeInTheDocument();
    expect(screen.getByText("stage: II")).toBeInTheDocument();
  });

  it("renders nothing extra with empty tooltipData", () => {
    const { container } = render(
      <HoverTooltip {...baseProps} tooltipData={{}} />
    );
    // Should have exactly 2 divs for x and y (plus the wrapper)
    const innerDivs = container.firstChild.querySelectorAll(":scope > div");
    expect(innerDivs).toHaveLength(2); // x and y only
  });
});
