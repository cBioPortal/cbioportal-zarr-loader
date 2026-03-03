import { useEffect, useLayoutEffect, useMemo, useRef } from "react";
import { Navigate, Outlet, useSearchParams, Link } from "react-router";
import {
  Layout,
  Spin,
  Alert,
  Tabs,
  Button,
} from "antd";
import { GithubOutlined, UploadOutlined } from "@ant-design/icons";
import ColumnsTab from "./components/views/ColumnsTab";
import InfoTab from "./components/views/InfoTab";
import ObsmTab from "./components/views/ObsmTab";
import PlotsTab from "./components/views/PlotsTab";
import DotplotTab from "./components/views/DotplotTab";
import ExplorerLayout from "./components/layouts/ExplorerLayout";
import { ProfileBar, PROFILE_BAR_HEIGHT, saveProfileSession } from "@cbioportal-zarr-loader/profiler";

import useAppStore from "./store/useAppStore";
import usePostMessage from "./hooks/usePostMessage";
import useIframeResize from "./hooks/useIframeResize";
import useLinkWithParams from "./hooks/useLinkWithParams";
import { saveRecentUrl } from "./utils/recentUrls";
import { DEFAULT_URL } from "./constants";

const { Header, Content } = Layout;

/**
 * ViewerLayout — pathless nested layout route.
 * Handles URL extraction, store initialization, loading/error guards.
 * Renders <Outlet /> (ViewerTabs) when data is ready.
 */
export function ViewerLayout() {
  const [searchParams] = useSearchParams();
  const url = searchParams.get("url") || DEFAULT_URL;

  const {
    loading,
    error,
    initialize,
  } = useAppStore();

  const initUrlRef = useRef(null);

  useEffect(() => {
    if (!url) return;
    if (initUrlRef.current === url) return;
    initUrlRef.current = url;
    initialize(url).then(() => {
      const { error } = useAppStore.getState();
      if (!error) saveRecentUrl(url);
    });
  }, [initialize, url]);

  if (!url) {
    return <Navigate to="/load" replace />;
  }

  if (loading) {
    return (
      <div style={{ padding: 24, textAlign: "center" }}>
        <Spin size="large" />
        <p style={{ marginTop: 16 }}>Loading AnnData from {url}...</p>
      </div>
    );
  }

  if (error) {
    return (
      <div style={{ padding: 24 }}>
        <Alert
          type="error"
          message="Error loading AnnData"
          description={
            <>
              <p>{error}</p>
              <p>Make sure the Zarr store is being served at {url}</p>
            </>
          }
        />
      </div>
    );
  }

  return <Outlet />;
}

/**
 * ViewerTabs — pure tabs UI.
 * Reads featureFlags and isEmbedded from the store, renders tab items.
 */
export function ViewerTabs() {
  const { featureFlags, isEmbedded } = useAppStore();

  const tabItems = [
    {
      key: "explorer",
      label: "Explore",
      children: <ObsmTab />,
    },
    {
      key: "columns",
      label: "Data",
      children: <ColumnsTab />,
    },
    {
      key: "plots",
      label: "Plots",
      children: <PlotsTab />,
    },
    ...(featureFlags.dotplot ? [{ key: "dotplot", label: "Dotplot", children: <DotplotTab /> }] : []),
    {
      key: "info",
      label: "Info",
      children: <InfoTab />,
    },
  ];

  return (
    <div style={{ padding: isEmbedded ? "0 24px 24px" : 24 }}>
      <Tabs items={tabItems} defaultActiveKey={import.meta.env.VITE_DEFAULT_TAB || "explorer"} />
    </div>
  );
}

/**
 * ViewerContent — layout switcher.
 * Reads ?layout=v3 to choose between tabbed layout and explorer layout.
 */
export function ViewerContent() {
  const [searchParams] = useSearchParams();
  const layoutVersion = searchParams.get("layout");

  if (layoutVersion === "v3") {
    return <ExplorerLayout />;
  }

  return <ViewerTabs />;
}

/**
 * App — root layout.
 * Owns app-level hooks (postMessage, iframe resize) so they survive route navigation.
 */
export default function App() {
  const { featureFlags, isEmbedded } = useAppStore();
  const adata = useAppStore((s) => s.adata);
  const url = useAppStore((s) => s.url);
  const linkTo = useLinkWithParams();
  const [searchParams] = useSearchParams();
  const isV3 = searchParams.get("layout") === "v3";

  const postMessageHandlers = useMemo(() => ({
    applyConfig: async (payload) => {
      const result = await useAppStore.getState().applyFilterConfig(payload);
      if (!result.success) console.error("[CZL:postMessage] applyConfig failed:", result.error);
    },
  }), []);

  usePostMessage(postMessageHandlers, import.meta.env.VITE_POSTMESSAGE_ORIGIN || "*");
  useIframeResize();

  // Apply viewport height constraint on #app only for v3 layout.
  // The old tab layout needs a scrollable page — hard-coding height: 100vh
  // on #app causes scrollbar toggle loops with the scatterplot resize handler.
  useLayoutEffect(() => {
    const el = document.getElementById("app");
    if (!el) return;
    if (isV3) {
      el.style.height = "100vh";
      el.style.overflow = "hidden";
    } else {
      el.style.height = "";
      el.style.overflow = "";
    }
  }, [isV3]);

  return (
    <Layout style={{ minHeight: "100vh", height: isV3 ? "100vh" : undefined }}>
      {!isEmbedded && !isV3 && (
        <Header
          style={{
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
            padding: "0 24px",
            background: "#fff",
            borderBottom: "1px solid #f0f0f0",
          }}
        >
          <Link to={linkTo("/")} style={{ fontSize: 18, fontWeight: 600, color: "inherit", textDecoration: "none" }}>
            cBioportal ZExplorer
          </Link>
          <nav style={{ display: "flex", gap: 16, alignItems: "center" }}>
            {featureFlags.loadDataset && (
              <Link to={linkTo("/load")}>
                <Button type="text" icon={<UploadOutlined />}>Load Dataset</Button>
              </Link>
            )}
            {featureFlags.profile && (
              <Link to={linkTo("/profile")}>
                <Button type="text">Profile History</Button>
              </Link>
            )}
            <a
              href="https://github.com/cbioportal/cbioportal-zarr-loader"
              target="_blank"
              rel="noopener noreferrer"
            >
              <GithubOutlined style={{ fontSize: 20 }} />
            </a>
          </nav>
        </Header>
      )}
      <Content style={{ background: "#fff", paddingBottom: featureFlags.profile ? PROFILE_BAR_HEIGHT : 0, ...(isV3 ? { flex: 1, overflow: "hidden" } : {}) }}>
        <Outlet />
      </Content>
      {featureFlags.profile && (
        <ProfileBar
          profiler={adata?.profiler}
          onSave={(entries) => {
            if (adata) saveProfileSession(url, adata.nObs, adata.nVar, entries);
          }}
        />
      )}
    </Layout>
  );
}
