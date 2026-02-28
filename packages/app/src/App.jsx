import { useEffect, useMemo } from "react";
import { Routes, Route, Navigate, useSearchParams, Link } from "react-router";
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
import LoadPage from "./pages/LoadPage";
import ProfilePage from "./pages/ProfilePage";
import ProfileBar, { PROFILE_BAR_HEIGHT } from "./components/ui/ProfileBar";

import useAppStore from "./store/useAppStore";
import usePostMessage from "./hooks/usePostMessage";
import useIframeResize from "./hooks/useIframeResize";
import useLinkWithParams from "./hooks/useLinkWithParams";
import { saveRecentUrl } from "./utils/recentUrls";
import { DEFAULT_URL } from "./constants";

const isEmbedded = window.self !== window.top || new URLSearchParams(window.location.search).has("embedded");

const { Header, Content } = Layout;

function ViewerContent() {
  const [searchParams] = useSearchParams();
  const url = searchParams.get("url") || DEFAULT_URL;

  const {
    loading,
    error,
    metadata,
    featureFlags,
    initialize,
  } = useAppStore();

  useEffect(() => {
    if (!url) return;
    // Skip re-initialization if we already have data for this URL
    const { url: currentUrl, adata } = useAppStore.getState();
    if (currentUrl === url && adata) return;
    initialize(url).then(() => {
      const { error } = useAppStore.getState();
      if (!error) saveRecentUrl(url);
    });
  }, [initialize, url]);

  if (!url) {
    return <Navigate to="/load" replace />;
  }

  const postMessageHandlers = useMemo(() => ({
    applyConfig: async (payload) => {
      const result = await useAppStore.getState().applyFilterConfig(payload);
      if (!result.success) console.error("[CZL:postMessage] applyConfig failed:", result.error);
    },
  }), []);

  usePostMessage(postMessageHandlers, import.meta.env.VITE_POSTMESSAGE_ORIGIN || "*");
  useIframeResize();

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

  const { obsColumns, varColumns } = metadata;

  const tabItems = [
    {
      key: "explorer",
      label: "Explore",
      children: <ObsmTab />,
    },
    {
      key: "columns",
      label: `Data (${obsColumns.length + varColumns.length})`,
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

export default function App() {
  const { featureFlags } = useAppStore();
  const linkTo = useLinkWithParams();

  return (
    <Layout style={{ minHeight: "100vh" }}>
      {!isEmbedded && (
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
      <Content style={{ background: "#fff", paddingBottom: featureFlags.profile ? PROFILE_BAR_HEIGHT : 0 }}>
        <Routes>
          <Route path="/load" element={<LoadPage />} />
          <Route path="/profile" element={<ProfilePage />} />
          <Route path="/*" element={<ViewerContent />} />
        </Routes>
      </Content>
      {featureFlags.profile && <ProfileBar />}
    </Layout>
  );
}
