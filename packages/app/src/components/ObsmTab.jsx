import { useEffect } from "react";
import {
  Row,
  Col,
  Card,
  Typography,
  Spin,
  Alert,
  Space,
  Button,
} from "antd";
import { ReloadOutlined } from "@ant-design/icons";
import EmbeddingScatterplot from "./EmbeddingScatterplot";
import StreamingEmbeddingScatterplot from "./StreamingEmbeddingScatterplot";
import GeneList from "./GeneList";
import useAppStore from "../store/useAppStore";

const { Text } = Typography;

export default function ObsmTab() {
  const {
    metadata,
    selectedObsm,
    obsmData,
    obsmLoading,
    obsmTime,
    fetchObsm,
    obsmStreamingData,
    obsmStreamingLoading,
  } = useAppStore();

  const { obsmKeys, obsColumns, geneNames } = metadata;
  const hasEmbeddingShape = obsmData?.shape?.[1] >= 2 || obsmStreamingData?.shape?.[1] >= 2;
  const isEmbedding = selectedObsm && /umap|tsne|pca/i.test(selectedObsm) && (obsmLoading || obsmStreamingLoading || hasEmbeddingShape);

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
            <div>
              {obsmKeys.map((k) => (
                <div key={k} style={{ padding: "4px 0" }}>
                  <Button
                    type={selectedObsm === k ? "primary" : "text"}
                    size="small"
                    onClick={() => fetchObsm(k)}
                  >
                    {k}
                  </Button>
                </div>
              ))}
            </div>
          ) : (
            <Text type="secondary">(none)</Text>
          )}
        </Card>
        <GeneList height={300} style={{ marginTop: 16 }} />
      </Col>
      <Col xs={24} md={18}>
        {selectedObsm ? (
          <Card
            title={`obsm: ${selectedObsm}`}
            size="small"
            extra={
              <Button
                size="small"
                icon={<ReloadOutlined />}
                onClick={() => fetchObsm(selectedObsm)}
                loading={obsmLoading || obsmStreamingLoading}
              >
                Reload
              </Button>
            }
          >
            {isEmbedding && (
              <div style={{ display: "flex", gap: 16 }}>
                <div style={{ flex: 1 }}>
                  <Text strong>Non-streaming</Text>
                  {obsmLoading ? (
                    <Spin style={{ display: "block", marginTop: 16 }} />
                  ) : obsmData?.error ? (
                    <Alert type="error" message={obsmData.error} />
                  ) : (
                    <>
                      <Space style={{ marginBottom: 8, marginTop: 8 }} wrap>
                        <Text>Fetched in {obsmTime?.toFixed(1)} ms</Text>
                        <Text type="secondary">|</Text>
                        <Text>Shape: {obsmData?.shape?.join(" × ")}</Text>
                      </Space>
                      <EmbeddingScatterplot
                        data={obsmData.data}
                        shape={obsmData.shape}
                        label={selectedObsm}
                      />
                    </>
                  )}
                </div>
                <div style={{ flex: 1 }}>
                  <Text strong>Streaming</Text>
                  {obsmStreamingLoading ? (
                    <Spin style={{ display: "block", marginTop: 16 }} />
                  ) : (
                    <div style={{ marginTop: 8 }}>
                      <StreamingEmbeddingScatterplot label={selectedObsm} />
                    </div>
                  )}
                </div>
              </div>
            )}
            {!isEmbedding && (
              obsmLoading ? (
                <Spin />
              ) : obsmData?.error ? (
                <Alert type="error" message={obsmData.error} />
              ) : (
                <Space style={{ marginBottom: 16 }} wrap>
                  <Text>Fetched in {obsmTime?.toFixed(1)} ms</Text>
                  <Text type="secondary">|</Text>
                  <Text>Shape: {obsmData?.shape?.join(" × ")}</Text>
                </Space>
              )
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
