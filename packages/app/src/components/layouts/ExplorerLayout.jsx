import { useEffect, useState } from "react";
import { Button, Drawer, Spin, Alert } from "antd";
import { MenuFoldOutlined, MenuUnfoldOutlined } from "@ant-design/icons";
import useMediaQuery from "../../hooks/useMediaQuery";
import useAppStore from "../../store/useAppStore";
import EmbeddingScatterplotContainer from "../containers/EmbeddingScatterplotContainer";
import ExplorerLeftSidebar from "./ExplorerLeftSidebar";
import ExplorerRightSidebar from "./ExplorerRightSidebar";

const SIDEBAR_WIDTH = 280;
const DESKTOP_QUERY = "(min-width: 768px)";

export default function ExplorerLayout() {
  const isDesktop = useMediaQuery(DESKTOP_QUERY);
  const [leftDrawerOpen, setLeftDrawerOpen] = useState(false);
  const [rightDrawerOpen, setRightDrawerOpen] = useState(false);

  const {
    metadata,
    selectedObsm,
    obsmData,
    obsmLoading,
    fetchObsm,
    featureFlags,
  } = useAppStore();

  const { obsmKeys } = metadata;
  const isEmbedding = selectedObsm && obsmData?.shape?.[1] >= 2;

  // Auto-fetch an embedding on mount (same logic as ObsmTab)
  useEffect(() => {
    if (!selectedObsm && obsmKeys.length > 0) {
      const defaultKey =
        obsmKeys.find((k) => /umap/i.test(k)) || obsmKeys[0];
      fetchObsm(defaultKey);
    }
  }, [obsmKeys, selectedObsm, fetchObsm]);

  const leftSidebar = <ExplorerLeftSidebar />;
  const rightSidebar = <ExplorerRightSidebar />;

  const centerContent = (
    <Spin spinning={obsmLoading} tip="Loading coordinates...">
      {obsmData?.error ? (
        <Alert type="error" message={obsmData.error} />
      ) : isEmbedding ? (
        <EmbeddingScatterplotContainer
          data={obsmData.data}
          shape={obsmData.shape}
          label={selectedObsm}
          showHexbinToggle={!!featureFlags.hexbin}
          debugMode={!!featureFlags.deckglDebug}
        />
      ) : (
        <div
          style={{
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            height: "100%",
            color: "#999",
          }}
        >
          {selectedObsm
            ? "Selected key is not a 2D embedding"
            : "Loading embedding..."}
        </div>
      )}
    </Spin>
  );

  if (isDesktop) {
    return (
      <div style={{ display: "flex", height: "100%" }}>
        <div
          style={{
            width: SIDEBAR_WIDTH,
            flexShrink: 0,
            borderRight: "1px solid #f0f0f0",
            overflow: "auto",
          }}
        >
          {leftSidebar}
        </div>
        <div style={{ flex: 1, minWidth: 0 }}>{centerContent}</div>
        <div
          style={{
            width: SIDEBAR_WIDTH,
            flexShrink: 0,
            borderLeft: "1px solid #f0f0f0",
            overflow: "auto",
          }}
        >
          {rightSidebar}
        </div>
      </div>
    );
  }

  // Mobile: full-width center with drawer sidebars
  return (
    <div style={{ height: "100%" }}>
      <div
        style={{
          display: "flex",
          justifyContent: "space-between",
          padding: "8px 16px",
          borderBottom: "1px solid #f0f0f0",
        }}
      >
        <Button
          icon={<MenuUnfoldOutlined />}
          onClick={() => setLeftDrawerOpen(true)}
          aria-label="Open left sidebar"
        />
        <Button
          icon={<MenuFoldOutlined />}
          onClick={() => setRightDrawerOpen(true)}
          aria-label="Open right sidebar"
        />
      </div>
      <div style={{ flex: 1 }}>{centerContent}</div>
      <Drawer
        placement="left"
        open={leftDrawerOpen}
        onClose={() => setLeftDrawerOpen(false)}
        width={SIDEBAR_WIDTH}
        styles={{ body: { padding: 0 } }}
      >
        {leftSidebar}
      </Drawer>
      <Drawer
        placement="right"
        open={rightDrawerOpen}
        onClose={() => setRightDrawerOpen(false)}
        width={SIDEBAR_WIDTH}
        styles={{ body: { padding: 0 } }}
      >
        {rightSidebar}
      </Drawer>
    </div>
  );
}
