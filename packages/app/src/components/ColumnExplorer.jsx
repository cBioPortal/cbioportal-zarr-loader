import { useMemo } from "react";
import {
  Card,
  Typography,
  Table,
  Space,
} from "antd";
import TabLayout from "./TabLayout";
import SearchableList from "./SearchableList";

const { Text } = Typography;

export default function ColumnExplorer({
  columns,
  selectedColumns,
  columnsData,
  index,
  loading,
  time,
  onToggleColumn,
  onClearAll,
}) {
  const lastSelected = selectedColumns.length
    ? selectedColumns[selectedColumns.length - 1]
    : null;

  const valueCounts = useMemo(() => {
    if (!lastSelected || !columnsData[lastSelected]) return [];
    const counts = {};
    for (const v of columnsData[lastSelected]) {
      const key = String(v);
      counts[key] = (counts[key] || 0) + 1;
    }
    return Object.entries(counts)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 15)
      .map(([value, count]) => ({ key: value, value, count }));
  }, [lastSelected, columnsData]);

  const tableColumns = useMemo(() => {
    const cols = [{ title: "Index", dataIndex: "index", key: "index" }];
    for (const col of selectedColumns) {
      cols.push({ title: col, dataIndex: col, key: col });
    }
    return cols;
  }, [selectedColumns]);

  const tableData = useMemo(() => {
    if (!index) return [];
    return index.map((id, i) => {
      const row = { key: id, index: id };
      for (const col of selectedColumns) {
        row[col] = columnsData[col] ? String(columnsData[col][i]) : "";
      }
      return row;
    });
  }, [index, selectedColumns, columnsData]);

  return (
    <TabLayout
      sidebar={
        <>
          <SearchableList
            title="Columns"
            items={columns}
            selected={selectedColumns}
            onSelect={onToggleColumn}
            onClear={onClearAll}
            loading={loading}
            multiSelect
            placeholder="Search columns..."
            height={300}
          />
          {lastSelected && (
            <Card
              title={`Value Counts: ${lastSelected}`}
              size="small"
              style={{ marginTop: 16 }}
            >
              {loading === lastSelected ? null : (
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
        </>
      }
    >
      <Card title="Data" size="small">
        <Space style={{ marginBottom: 16 }}>
          <Text>Rows: {index?.length?.toLocaleString() ?? 0}</Text>
          {time != null && (
            <>
              <Text type="secondary">|</Text>
              <Text>Last fetch: {time.toFixed(1)} ms</Text>
            </>
          )}
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
          columns={tableColumns}
        />
      </Card>
    </TabLayout>
  );
}
