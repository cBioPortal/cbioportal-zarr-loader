import { describe, it, expect, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import ProfilePage from "./ProfilePage";

describe("ProfilePage", () => {
  beforeEach(() => {
    localStorage.clear();
  });

  it("renders empty state when no history exists", () => {
    render(<ProfilePage />);
    expect(screen.getByText(/no profile sessions saved/i)).toBeTruthy();
  });

  it("renders the page title", () => {
    render(<ProfilePage />);
    expect(screen.getAllByText("Profile History").length).toBeGreaterThan(0);
  });

  it("renders action buttons", () => {
    render(<ProfilePage />);
    const buttons = screen.getAllByRole("button");
    const buttonTexts = buttons.map((b) => b.textContent);
    expect(buttonTexts.some((t) => t?.includes("Download"))).toBe(true);
    expect(buttonTexts.some((t) => t?.includes("Import"))).toBe(true);
  });
});
