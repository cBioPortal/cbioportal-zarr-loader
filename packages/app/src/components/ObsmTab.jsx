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

const { Title, Text } = Typography;

export default function ObsmTab({ obsmKeys, obsm, obsColumns, fetchColorData }) {
  const { selected, data, loading, time, fetch } = obsm;
  const isEmbedding = selected && /umap|tsne|pca/i.test(selected) && data?.shape?.[1] >= 2;

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
                    type={selected === k ? "primary" : "text"}
                    size="small"
                    onClick={() => fetch(k)}
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
        {selected ? (
          <Card title={`obsm: ${selected}`} size="small">
            {loading ? (
              <Spin />
            ) : data?.error ? (
              <Alert type="error" message={data.error} />
            ) : (
              <>
                <Space style={{ marginBottom: 16 }} wrap>
                  <Text>Fetched in {time?.toFixed(1)} ms</Text>
                  <Text type="secondary">|</Text>
                  <Text>Shape: {data?.shape?.join(" × ")}</Text>
                </Space>

                {isEmbedding && (
                  <EmbeddingScatterplot
                    data={data.data}
                    shape={data.shape}
                    label={selected}
                    obsColumns={obsColumns}
                    fetchColorData={fetchColorData}
                  />
                )}

                <Title level={5}>First 10 rows</Title>
                <Table
                  size="small"
                  pagination={false}
                  scroll={{ x: true }}
                  dataSource={data?.index?.slice(0, 10).map((id, i) => {
                    const row = { key: id, index: id };
                    for (let j = 0; j < data.shape[1]; j++) {
                      row[`col${j}`] = data.data[i * data.shape[1] + j]?.toFixed(4);
                    }
                    return row;
                  })}
                  columns={[
                    { title: "Index", dataIndex: "index", key: "index", fixed: "left" },
                    ...Array.from({ length: data?.shape?.[1] || 0 }, (_, j) => ({
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
