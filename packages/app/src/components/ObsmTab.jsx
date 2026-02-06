import { useEffect } from "react";
import {
  Row,
  Col,
  Card,
  Typography,
  Spin,
  Alert,
  Table,
  List,
  Space,
  Button,
} from "antd";
import EmbeddingScatterplot from "./EmbeddingScatterplot";
import useAppStore from "../store/useAppStore";

const { Title, Text } = Typography;

export default function ObsmTab() {
  const {
    metadata,
    selectedObsm,
    obsmData,
    obsmLoading,
    obsmTime,
    fetchObsm,
  } = useAppStore();

  const { obsmKeys, obsColumns, geneNames } = metadata;
  const isEmbedding = selectedObsm && /umap|tsne|pca/i.test(selectedObsm) && obsmData?.shape?.[1] >= 2;

  // Auto-fetch UMAP embedding on mount
  useEffect(() => {
    if (!selectedObsm && obsmKeys.length > 0) {
      const umapKey = obsmKeys.find(k => /umap/i.test(k));
      if (umapKey) {
        fetchObsm(umapKey);
      }
    }
  }, [obsmKeys, selectedObsm, fetchObsm]);

  return (
    <Row gutter={[16, 16]}>
      <Col xs={24} md={6}>
        <Card title="Keys" size="small">
          {obsmKeys.length ? (
            <List
              size="small"
              dataSource={obsmKeys}
              renderItem={(k) => (
                <List.Item style={{ padding: "4px 0" }}>
                  <Button
                    type={selectedObsm === k ? "primary" : "text"}
                    size="small"
                    onClick={() => fetchObsm(k)}
                  >
                    {k}
                  </Button>
                </List.Item>
              )}
            />
          ) : (
            <Text type="secondary">(none)</Text>
          )}
        </Card>
      </Col>
      <Col xs={24} md={18}>
        {selectedObsm ? (
          <Card title={`obsm: ${selectedObsm}`} size="small">
            {obsmLoading ? (
              <Spin />
            ) : obsmData?.error ? (
              <Alert type="error" message={obsmData.error} />
            ) : (
              <>
                <Space style={{ marginBottom: 16 }} wrap>
                  <Text>Fetched in {obsmTime?.toFixed(1)} ms</Text>
                  <Text type="secondary">|</Text>
                  <Text>Shape: {obsmData?.shape?.join(" × ")}</Text>
                </Space>

                {isEmbedding && (
                  <EmbeddingScatterplot
                    data={obsmData.data}
                    shape={obsmData.shape}
                    label={selectedObsm}
                  />
                )}

                <Title level={5}>First 10 rows</Title>
                <Table
                  size="small"
                  pagination={false}
                  scroll={{ x: true }}
                  dataSource={obsmData?.index?.slice(0, 10).map((id, i) => {
                    const row = { key: id, index: id };
                    for (let j = 0; j < obsmData.shape[1]; j++) {
                      row[`col${j}`] = obsmData.data[i * obsmData.shape[1] + j]?.toFixed(4);
                    }
                    return row;
                  })}
                  columns={[
                    { title: "Index", dataIndex: "index", key: "index", fixed: "left" },
                    ...Array.from({ length: obsmData?.shape?.[1] || 0 }, (_, j) => ({
                      title: String(j),
                      dataIndex: `col${j}`,
                      key: `col${j}`,
                    })),
                  ]}
                />
              </>
            )}
          </Card>
        ) : (
          <Card size="small">
            <Text type="secondary">Select a key to view its data</Text>
          </Card>
        )}
      </Col>
    </Row>
  );
}
