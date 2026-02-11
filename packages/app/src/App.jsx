import { useEffect } from "react";
import {
  Typography,
  Spin,
  Alert,
  Tabs,
} from "antd";
import ObsPage from "./pages/ObsPage";
import VarPage from "./pages/VarPage";
import ObsmPage from "./pages/ObsmPage";
import LayersPage from "./pages/LayersPage";
import PlotsPage from "./pages/PlotsPage";
import InfoPage from "./pages/InfoPage";
import useAppStore from "./store/useAppStore";

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

  const { obsColumns, varColumns, obsmKeys, layerKeys } = metadata;

  const tabItems = [
    {
      key: "obs",
      label: `obs (${obsColumns.length})`,
      children: <ObsPage />,
    },
    {
      key: "var",
      label: `var (${varColumns.length})`,
      children: <VarPage />,
    },
    {
      key: "obsm",
      label: `obsm (${obsmKeys.length})`,
      children: <ObsmPage />,
    },
    {
      key: "layers",
      label: `layers (${layerKeys.length})`,
      children: <LayersPage />,
    },
    {
      key: "plots",
      label: "Plots",
      children: <PlotsPage />,
    },
    {
      key: "info",
      label: "Info",
      children: <InfoPage />,
    },
  ];

  return (
    <div style={{ padding: 24 }}>
      <Title level={3}>AnnData Zarr Loader</Title>
      <Tabs items={tabItems} defaultActiveKey="obsm" />
    </div>
  );
}
