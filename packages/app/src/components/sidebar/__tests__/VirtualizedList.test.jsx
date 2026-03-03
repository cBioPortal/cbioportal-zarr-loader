import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import VirtualizedList from "../VirtualizedList";

describe("VirtualizedList", () => {
  const items = ["alpha", "beta", "gamma", "delta", "epsilon"];

  it("renders visible items", () => {
    render(
      <VirtualizedList items={items} selected={null} onSelect={() => {}} height={200} />
    );
    expect(screen.getByText("alpha")).toBeInTheDocument();
  });

  it("highlights the selected item", () => {
    render(
      <VirtualizedList items={items} selected="beta" onSelect={() => {}} height={200} />
    );
    // react-window v2 may render items more than once; find the one with the highlight
    const rows = screen.getAllByTitle("beta");
    const highlighted = rows.find(
      (el) => el.style.backgroundColor === "rgb(230, 244, 255)"
    );
    expect(highlighted).toBeTruthy();
  });

  it("calls onSelect when an item is clicked", async () => {
    const onSelect = vi.fn();
    render(
      <VirtualizedList items={items} selected={null} onSelect={onSelect} height={200} />
    );
    // react-window v2 may render items more than once; click the last (interactive) one
    const targets = screen.getAllByText("gamma");
    await userEvent.click(targets[targets.length - 1]);
    expect(onSelect).toHaveBeenCalledWith("gamma");
  });

  it("applies loading styles to loading item", () => {
    render(
      <VirtualizedList items={items} selected={null} onSelect={() => {}} height={200} loading="alpha" />
    );
    // react-window v2 may render items more than once; find the one with loading styles
    const rows = screen.getAllByTitle("alpha");
    const loadingRow = rows.find((el) => el.style.opacity === "0.5");
    expect(loadingRow).toBeTruthy();
    expect(loadingRow.style.cursor).toBe("wait");
  });

  it("renders nothing when items is empty", () => {
    const { container } = render(
      <VirtualizedList items={[]} selected={null} onSelect={() => {}} height={200} />
    );
    expect(container.firstChild).toBeNull();
  });
});
