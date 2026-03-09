import React from "react";
import ReactDOM from "react-dom/client";
import { createBrowserRouter, RouterProvider } from "react-router";
import App, { ViewerLayout, ViewerContent } from "./App.jsx";
import LoadPage from "./pages/LoadPage";
import { ProfilePage } from "@cbioportal-cell-explorer/profiler";

const router = createBrowserRouter([
  {
    path: "/",
    Component: App,
    children: [
      { path: "load", Component: LoadPage },
      { path: "profile", Component: ProfilePage },
      {
        Component: ViewerLayout,
        children: [
          { index: true, Component: ViewerContent },
          { path: "*", Component: ViewerContent },
        ],
      },
    ],
  },
], {
  basename: import.meta.env.BASE_URL,
});

ReactDOM.createRoot(document.getElementById("app")).render(
  <React.StrictMode>
    <RouterProvider router={router} />
  </React.StrictMode>,
);
