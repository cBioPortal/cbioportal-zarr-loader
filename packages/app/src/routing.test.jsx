import { describe, it, expect, afterEach, vi } from "vitest";
import { render, screen, cleanup } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { createMemoryRouter, RouterProvider } from "react-router";
import App, { ViewerLayout, ViewerTabs } from "./App";

// Stub heavy child components so tests stay fast and focused on routing
vi.mock("./components/views/ObsmTab", () => ({ default: () => <div>ObsmTab</div> }));
vi.mock("./components/views/ColumnsTab", () => ({ default: () => <div>ColumnsTab</div> }));
vi.mock("./components/views/PlotsTab", () => ({ default: () => <div>PlotsTab</div> }));
vi.mock("./components/views/DotplotTab", () => ({ default: () => <div>DotplotTab</div> }));
vi.mock("./components/views/InfoTab", () => ({ default: () => <div>InfoTab</div> }));
vi.mock("./hooks/usePostMessage", () => ({ default: () => {} }));
vi.mock("./hooks/useIframeResize", () => ({ default: () => {} }));
vi.mock("./store/useAppStore", () => {
  const store = () => ({
    loading: false,
    error: null,
    isEmbedded: false,
    featureFlags: { loadDataset: true, profile: true, dotplot: false },
    initialize: vi.fn().mockResolvedValue(undefined),
    adata: null,
    url: null,
  });
  store.getState = store;
  store.setState = vi.fn();
  store.subscribe = vi.fn(() => vi.fn());
  return { default: store };
});
vi.mock("@cbioportal-zarr-loader/profiler", () => ({
  ProfilePage: () => <div data-testid="profile-page">Profile Page</div>,
  ProfileBar: ({ renderLink }) => renderLink ? renderLink("Profile History") : null,
  PROFILE_BAR_HEIGHT: 0,
  saveProfileSession: vi.fn(),
}));
vi.mock("./pages/LoadPage", () => ({
  default: () => <div data-testid="load-page">Load Page</div>,
}));
vi.mock("./constants", () => ({ DEFAULT_URL: "http://test.example/data.zarr" }));

afterEach(cleanup);

function renderWithRouter(initialPath = "/") {
  const router = createMemoryRouter(
    [
      {
        element: <App />,
        children: [
          { path: "load", element: <LoadPage /> },
          { path: "profile", element: <ProfilePage /> },
          {
            element: <ViewerLayout />,
            children: [
              { index: true, element: <ViewerTabs /> },
              { path: "*", element: <ViewerTabs /> },
            ],
          },
        ],
      },
    ],
    { initialEntries: [initialPath] },
  );
  return render(<RouterProvider router={router} />);
}

// Pull mocked components into scope for the route definitions
import LoadPage from "./pages/LoadPage";
import { ProfilePage } from "@cbioportal-zarr-loader/profiler";

describe("Routing", () => {
  it("renders ViewerTabs at the root path", async () => {
    renderWithRouter("/");
    expect(screen.getByText(/Loading AnnData|ObsmTab/)).toBeInTheDocument();
  });

  it("renders LoadPage at /load", () => {
    renderWithRouter("/load");
    expect(screen.getByTestId("load-page")).toBeInTheDocument();
  });

  it("renders ProfilePage at /profile", () => {
    renderWithRouter("/profile");
    expect(screen.getByTestId("profile-page")).toBeInTheDocument();
  });

  it("navigates from root to /load when Load Dataset is clicked", async () => {
    const user = userEvent.setup();
    renderWithRouter("/?url=http://test.example/data.zarr");

    const loadLink = screen.getByRole("link", { name: /load dataset/i });
    await user.click(loadLink);

    expect(screen.getByTestId("load-page")).toBeInTheDocument();
  });

  it("navigates from root to /profile when Profile History is clicked", async () => {
    const user = userEvent.setup();
    renderWithRouter("/?url=http://test.example/data.zarr");

    const profileLink = screen.getByRole("link", { name: /profile history/i });
    await user.click(profileLink);

    expect(screen.getByTestId("profile-page")).toBeInTheDocument();
  });

  it("navigates back to root from /load when logo is clicked", async () => {
    const user = userEvent.setup();
    renderWithRouter("/load?url=http://test.example/data.zarr");

    expect(screen.getByTestId("load-page")).toBeInTheDocument();

    const logo = screen.getByRole("link", { name: /zexplorer/i });
    await user.click(logo);

    expect(screen.queryByTestId("load-page")).not.toBeInTheDocument();
  });
});
