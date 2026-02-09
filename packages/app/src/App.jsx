import { useEffect } from "react";
import {
  Card,
  Typography,
  Spin,
  Alert,
  Tabs,
  Descriptions,
} from "antd";
import ColumnExplorer from "./components/ColumnExplorer";
import ObsmTab from "./components/ObsmTab";
import PlotsTab from "./components/PlotsTab";

import useAppStore from "./store/useAppStore";

const { Title, Text } = Typography;

const URL = "https://cbioportal-public-imaging.assets.cbioportal.org/msk_spectrum_tme_2022/zarr/spectrum_all_cells.zarr";

export default function App() {
  const {
    adata,
    loading,
    error,
    metadata,
    initialize,
    // Obs column (multi-select)
    obsColumnsSelected,
    obsColumnsData,
    obsColumnLoading,
    obsColumnTime,
    obsIndex,
    toggleObsColumn,
    clearObsColumns,
    // Var column (multi-select)
    varColumnsSelected,
    varColumnsData,
    varColumnLoading,
    varColumnTime,
    varIndex,
    toggleVarColumn,
    clearVarColumns,
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

  const { obsColumns, varColumns, obsmKeys, layerKeys, chunks } = metadata;

  const tabItems = [
    {
      key: "obs",
      label: `obs (${obsColumns.length})`,
      children: (
        <ColumnExplorer
          columns={obsColumns}
          selectedColumns={obsColumnsSelected}
          columnsData={obsColumnsData}
          index={obsIndex}
          loading={obsColumnLoading}
          time={obsColumnTime}
          onToggleColumn={toggleObsColumn}
          onClearAll={clearObsColumns}
        />
      ),
    },
    {
      key: "var",
      label: `var (${varColumns.length})`,
      children: (
        <ColumnExplorer
          columns={varColumns}
          selectedColumns={varColumnsSelected}
          columnsData={varColumnsData}
          index={varIndex}
          loading={varColumnLoading}
          time={varColumnTime}
          onToggleColumn={toggleVarColumn}
          onClearAll={clearVarColumns}
        />
      ),
    },
    {
      key: "obsm",
      label: `obsm (${obsmKeys.length})`,
      children: <ObsmTab />,
    },
    {
      key: "layers",
      label: `layers (${layerKeys.length})`,
      children: (
        <Card title="Layers" size="small">
          {layerKeys.length ? (
            <div>
              {layerKeys.map((k) => (
                <div key={k} style={{ padding: "4px 0" }}>{k}</div>
              ))}
            </div>
          ) : (
            <Text type="secondary">(none)</Text>
          )}
        </Card>
      ),
    },
    {
      key: "plots",
      label: "Plots",
      children: <PlotsTab />,
    },
    {
      key: "info",
      label: "Info",
      children: (
        <Card title="Dataset" size="small">
          <Descriptions column={1} size="small">
            <Descriptions.Item label="Shape">
              {adata.nObs.toLocaleString()} obs × {adata.nVar.toLocaleString()} var
            </Descriptions.Item>
            <Descriptions.Item label="Chunk size">
              {chunks ? chunks.join(" × ") : "N/A"}
            </Descriptions.Item>
            <Descriptions.Item label="Encoding">
              {adata.attrs["encoding-type"]} v{adata.attrs["encoding-version"]}
            </Descriptions.Item>
            <Descriptions.Item label="URL">
              <Text copyable style={{ fontSize: 12 }}>{URL}</Text>
            </Descriptions.Item>
          </Descriptions>
        </Card>
      ),
    },
  ];

  return (
    <div style={{ padding: 24 }}>
      <Title level={3}>AnnData Zarr Loader</Title>
      <Tabs items={tabItems} defaultActiveKey="obsm" />
    </div>
  );
}
