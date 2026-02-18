import { useEffect, useMemo } from "react";
import { Routes, Route } from "react-router";
import {
  Typography,
  Spin,
  Alert,
  Tabs,
} from "antd";
import ColumnsTab from "./components/ColumnsTab";
import InfoTab from "./components/InfoTab";
import ObsmTab from "./components/ObsmTab";
import PlotsTab from "./components/PlotsTab";

import useAppStore from "./store/useAppStore";
import usePostMessage from "./hooks/usePostMessage";

const { Title } = Typography;

const URL = "https://cbioportal-public-imaging.assets.cbioportal.org/msk_spectrum_tme_2022/zarr/spectrum_all_cells.zarr";

export default function App() {
  const {
    loading,
    error,
    metadata,
    initialize,
  } = useAppStore();

  useEffect(() => {
    initialize(URL);
  }, [initialize]);

  const postMessageHandlers = useMemo(() => ({
    applyConfig: async (payload) => {
      const result = await useAppStore.getState().applyFilterConfig(payload);
      if (!result.success) console.error("[CZL:postMessage] applyConfig failed:", result.error);
    },
  }), []);

  usePostMessage(postMessageHandlers, import.meta.env.VITE_POSTMESSAGE_ORIGIN || "*");

  if (loading) {
    return (
      <div style={{ padding: 24, textAlign: "center" }}>
        <Spin size="large" />
        <p style={{ marginTop: 16 }}>Loading AnnData from {URL}...</p>
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
              <p>Make sure the Zarr store is being served at {URL}</p>
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
    {
      key: "info",
      label: "Info",
      children: <InfoTab />,
    },
  ];

  return (
    <Routes>
      <Route
        path="/*"
        element={
          <div style={{ padding: 24 }}>
            <Title level={3}>cBioportal ZExplorer</Title>
            <Tabs items={tabItems} defaultActiveKey="explorer" />
          </div>
        }
      />
    </Routes>
  );
}
