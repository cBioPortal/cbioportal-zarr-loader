import {
  Row,
  Col,
  Card,
  Typography,
  Spin,
  Alert,
  Table,
  List,
  Tabs,
  Descriptions,
} from "antd";
import ColumnExplorer from "./components/ColumnExplorer";
import ObsmTab from "./components/ObsmTab";
import useAnnData from "./hooks/useAnnData";

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
    obsColumn,
    varColumn,
    obsm,
    fetchColorData,
  } = useAnnData(URL);

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

  const { obsColumns, varColumns, obsmKeys, layerKeys, timings, chunks } = metadata;

  const timingsData = Object.entries(timings).map(([key, ms]) => ({
    key,
    operation: key,
    time: `${ms.toFixed(1)} ms`,
  }));
  timingsData.push({
    key: "total",
    operation: "Total",
    time: `${Object.values(timings).reduce((a, b) => a + b, 0).toFixed(1)} ms`,
  });

  const tabItems = [
    {
      key: "obs",
      label: `obs (${obsColumns.length})`,
      children: (
        <ColumnExplorer
          columns={obsColumns}
          selectedColumn={obsColumn.selected}
          columnData={obsColumn.data}
          loading={obsColumn.loading}
          time={obsColumn.time}
          onSelectColumn={obsColumn.fetch}
        />
      ),
    },
    {
      key: "var",
      label: `var (${varColumns.length})`,
      children: (
        <ColumnExplorer
          columns={varColumns}
          selectedColumn={varColumn.selected}
          columnData={varColumn.data}
          loading={varColumn.loading}
          time={varColumn.time}
          onSelectColumn={varColumn.fetch}
        />
      ),
    },
    {
      key: "obsm",
      label: `obsm (${obsmKeys.length})`,
      children: (
        <ObsmTab
          obsmKeys={obsmKeys}
          obsm={obsm}
          obsColumns={obsColumns}
          fetchColorData={fetchColorData}
        />
      ),
    },
    {
      key: "layers",
      label: `layers (${layerKeys.length})`,
      children: (
        <Card title="Layers">
          {layerKeys.length ? (
            <List
              size="small"
              dataSource={layerKeys}
              renderItem={(k) => <List.Item>{k}</List.Item>}
            />
          ) : (
            <Text type="secondary">(none)</Text>
          )}
        </Card>
      ),
    },
  ];

  return (
    <div style={{ padding: 24 }}>
      <Title level={2}>AnnData Store</Title>
      <Text type="secondary" copyable>
        {URL}
      </Text>

      <Row gutter={[16, 16]} style={{ marginTop: 24 }}>
        <Col xs={24} md={12}>
          <Card title="Dataset" size="small">
            <Descriptions column={1} size="small">
              <Descriptions.Item label="Shape">
                {adata.nObs.toLocaleString()} × {adata.nVar.toLocaleString()}
              </Descriptions.Item>
              <Descriptions.Item label="Chunk size">
                {chunks ? chunks.join(" × ") : "N/A"}
              </Descriptions.Item>
              <Descriptions.Item label="Encoding">
                {adata.attrs["encoding-type"]} v{adata.attrs["encoding-version"]}
              </Descriptions.Item>
            </Descriptions>
          </Card>
        </Col>
        <Col xs={24} md={12}>
          <Card title="Fetch Keys Timings" size="small">
            <Table
              size="small"
              pagination={false}
              dataSource={timingsData}
              columns={[
                { title: "Operation", dataIndex: "operation", key: "operation" },
                { title: "Time", dataIndex: "time", key: "time" },
              ]}
            />
          </Card>
        </Col>
      </Row>

      <Tabs items={tabItems} style={{ marginTop: 24 }} />
    </div>
  );
}
