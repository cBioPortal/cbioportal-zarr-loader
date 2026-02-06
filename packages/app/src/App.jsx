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
import useAppStore from "./store/useAppStore";

const { Title, Text } = Typography;

// const URL = "http://localhost:3000/pbmc3k.zarr";
// const URL = "http://localhost:3000/spectrum_all_cells.zarr";
const URL = "https://cbioportal-public-imaging.assets.cbioportal.org/msk_spectrum_tme_2022/zarr/spectrum_all_cells.zarr";

export default function App() {
  const {
    adata,
    loading,
    error,
    metadata,
    initialize,
    // Obs column
    selectedObsColumn,
    obsColumnData,
    obsColumnLoading,
    obsColumnTime,
    fetchObsColumn,
    // Var column
    selectedVarColumn,
    varColumnData,
    varColumnLoading,
    varColumnTime,
    fetchVarColumn,
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
          selectedColumn={selectedObsColumn}
          columnData={obsColumnData}
          loading={obsColumnLoading}
          time={obsColumnTime}
          onSelectColumn={fetchObsColumn}
        />
      ),
    },
    {
      key: "var",
      label: `var (${varColumns.length})`,
      children: (
        <ColumnExplorer
          columns={varColumns}
          selectedColumn={selectedVarColumn}
          columnData={varColumnData}
          loading={varColumnLoading}
          time={varColumnTime}
          onSelectColumn={fetchVarColumn}
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
