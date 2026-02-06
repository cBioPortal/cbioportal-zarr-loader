import { useMemo } from "react";
import {
  Row,
  Col,
  Card,
  Typography,
  Button,
  Spin,
  Alert,
  Table,
  List,
  Space,
} from "antd";

const { Text } = Typography;

export default function ColumnExplorer({
  columns,
  selectedColumn,
  columnData,
  loading,
  time,
  onSelectColumn,
}) {
  const valueCounts = useMemo(() => {
    if (!columnData?.values) return [];
    const counts = {};
    for (const v of columnData.values) {
      const key = String(v);
      counts[key] = (counts[key] || 0) + 1;
    }
    return Object.entries(counts)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 15)
      .map(([value, count]) => ({ key: value, value, count }));
  }, [columnData?.values]);

  const tableData = useMemo(() => {
    if (!columnData?.index) return [];
    return columnData.index.map((id, i) => ({
      key: id,
      index: id,
      value: String(columnData.values[i]),
    }));
  }, [columnData]);

  return (
    <Row gutter={[16, 16]}>
      <Col xs={24} md={6}>
        <Card title="Columns" size="small">
          <List
            size="small"
            dataSource={columns}
            style={{ maxHeight: 200, overflow: "auto" }}
            renderItem={(c) => (
              <List.Item style={{ padding: "4px 0" }}>
                <Button
                  type={selectedColumn === c ? "primary" : "text"}
                  size="small"
                  onClick={() => onSelectColumn(c)}
                >
                  {c}
                </Button>
              </List.Item>
            )}
          />
        </Card>
        {selectedColumn && (
          <Card
            title={`Value Counts: ${selectedColumn}`}
            size="small"
            style={{ marginTop: 16 }}
          >
            {loading ? (
              <Spin />
            ) : columnData?.error ? (
              <Alert type="error" message={columnData.error} />
            ) : (
              <Table
                size="small"
                pagination={false}
                style={{ maxHeight: 250, overflow: "auto" }}
                dataSource={valueCounts}
                columns={[
                  { title: "Value", dataIndex: "value", key: "value" },
                  { title: "Count", dataIndex: "count", key: "count" },
                ]}
              />
            )}
          </Card>
        )}
      </Col>
      <Col xs={24} md={18}>
        {selectedColumn ? (
          <Card title={`Column: ${selectedColumn}`} size="small">
            {loading ? (
              <Spin />
            ) : columnData?.error ? (
              <Alert type="error" message={columnData.error} />
            ) : (
              <>
                <Space style={{ marginBottom: 16 }}>
                  <Text>Fetched in {time?.toFixed(1)} ms</Text>
                  <Text type="secondary">|</Text>
                  <Text>Length: {columnData?.values?.length?.toLocaleString()}</Text>
                </Space>
                <Table
                  size="small"
                  pagination={{
                    defaultPageSize: 10,
                    showSizeChanger: true,
                    pageSizeOptions: [10, 25, 50, 100],
                    showTotal: (total, range) => `${range[0]}-${range[1]} of ${total}`,
                  }}
                  dataSource={tableData}
                  columns={[
                    { title: "Index", dataIndex: "index", key: "index" },
                    { title: selectedColumn, dataIndex: "value", key: "value" },
                  ]}
                />
              </>
            )}
          </Card>
        ) : (
          <Card size="small">
            <Text type="secondary">Select a column to view its data</Text>
          </Card>
        )}
      </Col>
    </Row>
  );
}
