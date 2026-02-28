import { useState, useMemo } from "react";
import { Table, Tag, Button, Card, Space, Typography, Popconfirm, Statistic, Row, Col } from "antd";
import { DeleteOutlined, ClearOutlined } from "@ant-design/icons";
import {
  getProfileHistory,
  removeProfileSession,
  clearProfileHistory,
} from "../utils/profileStorage";
import MethodBreakdownChart from "../components/charts/MethodBreakdownChart";
import CacheEfficiencyChart from "../components/charts/CacheEfficiencyChart";

const { Text, Title } = Typography;

const entryColumns = [
  { title: "#", dataIndex: "id", key: "id", width: 50 },
  { title: "Method", dataIndex: "method", key: "method", width: 120 },
  { title: "Key", dataIndex: "key", key: "key", ellipsis: true },
  {
    title: "Cache",
    dataIndex: "cacheHit",
    key: "cacheHit",
    width: 80,
    render: (hit) => <Tag color={hit ? "green" : "red"}>{hit ? "HIT" : "MISS"}</Tag>,
  },
  {
    title: "Duration",
    dataIndex: "duration",
    key: "duration",
    width: 100,
    render: (ms) => `${ms.toFixed(1)} ms`,
    sorter: (a, b) => a.duration - b.duration,
  },
];

function SummaryStats({ history }) {
  const stats = useMemo(() => {
    let totalQueries = 0;
    let totalHits = 0;
    let slowestQuery = 0;

    for (const session of history) {
      for (const e of session.entries) {
        totalQueries++;
        if (e.cacheHit) totalHits++;
        if (e.duration > slowestQuery) slowestQuery = e.duration;
      }
    }

    const avgHitRate = totalQueries > 0 ? ((totalHits / totalQueries) * 100).toFixed(1) : 0;

    return { totalSessions: history.length, totalQueries, avgHitRate, slowestQuery };
  }, [history]);

  return (
    <Row gutter={16} style={{ marginBottom: 24 }}>
      <Col span={6}>
        <Card size="small">
          <Statistic title="Total Sessions" value={stats.totalSessions} />
        </Card>
      </Col>
      <Col span={6}>
        <Card size="small">
          <Statistic title="Total Queries" value={stats.totalQueries} />
        </Card>
      </Col>
      <Col span={6}>
        <Card size="small">
          <Statistic title="Avg Cache Hit Rate" value={stats.avgHitRate} suffix="%" />
        </Card>
      </Col>
      <Col span={6}>
        <Card size="small">
          <Statistic title="Slowest Query" value={stats.slowestQuery.toFixed(1)} suffix="ms" />
        </Card>
      </Col>
    </Row>
  );
}

function Charts({ history }) {
  const sessions = useMemo(
    () =>
      history.map((s, i) => ({
        label: s.url ? new URL(s.url).pathname.split("/").pop() || `Session ${i + 1}` : `Session ${i + 1}`,
        entries: s.entries,
      })),
    [history],
  );

  return (
    <Row gutter={24} style={{ marginBottom: 24 }}>
      <Col span={12}>
        <Card size="small" title="Duration by Method">
          <MethodBreakdownChart sessions={sessions} width={440} />
        </Card>
      </Col>
      <Col span={12}>
        <Card size="small" title="Cache Efficiency">
          <CacheEfficiencyChart sessions={sessions} width={440} />
        </Card>
      </Col>
    </Row>
  );
}

export default function ProfilePage() {
  const [history, setHistory] = useState(() => getProfileHistory());

  const handleRemove = (index) => {
    removeProfileSession(index);
    setHistory(getProfileHistory());
  };

  const handleClearAll = () => {
    clearProfileHistory();
    setHistory([]);
  };

  const sessionColumns = [
    {
      title: "URL",
      dataIndex: "url",
      key: "url",
      ellipsis: true,
      render: (url) => <Text style={{ maxWidth: 200 }} ellipsis={{ tooltip: url }}>{url}</Text>,
    },
    {
      title: "Date",
      dataIndex: "timestamp",
      key: "timestamp",
      width: 180,
      render: (ts) => new Date(ts).toLocaleString(),
      sorter: (a, b) => a.timestamp - b.timestamp,
      defaultSortOrder: "descend",
    },
    {
      title: "Shape",
      key: "shape",
      width: 140,
      render: (_, r) => `${r.nObs?.toLocaleString()} x ${r.nVar?.toLocaleString()}`,
    },
    {
      title: "Queries",
      key: "queries",
      width: 80,
      render: (_, r) => r.entries.length,
    },
    {
      title: "Total Time",
      key: "totalTime",
      width: 110,
      render: (_, r) => `${r.entries.reduce((s, e) => s + e.duration, 0).toFixed(1)} ms`,
      sorter: (a, b) =>
        a.entries.reduce((s, e) => s + e.duration, 0) -
        b.entries.reduce((s, e) => s + e.duration, 0),
    },
    {
      title: "Hit Rate",
      key: "hitRate",
      width: 90,
      render: (_, r) => {
        const hits = r.entries.filter((e) => e.cacheHit).length;
        return r.entries.length > 0 ? `${((hits / r.entries.length) * 100).toFixed(0)}%` : "-";
      },
    },
    {
      title: "",
      key: "actions",
      width: 50,
      render: (_, __, index) => (
        <Button
          type="text"
          danger
          size="small"
          icon={<DeleteOutlined />}
          onClick={() => handleRemove(index)}
        />
      ),
    },
  ];

  return (
    <div style={{ padding: 24 }}>
      <Space style={{ marginBottom: 16, width: "100%", justifyContent: "space-between" }}>
        <Title level={4} style={{ margin: 0 }}>Profile History</Title>
        {history.length > 0 && (
          <Popconfirm title="Clear all profile history?" onConfirm={handleClearAll}>
            <Button icon={<ClearOutlined />} danger>Clear All</Button>
          </Popconfirm>
        )}
      </Space>
      {history.length === 0 ? (
        <Card>
          <Text type="secondary">No profile sessions saved yet. Use the profiler drawer to record and save sessions.</Text>
        </Card>
      ) : (
        <>
          <SummaryStats history={history} />
          <Charts history={history} />
          <Table
            size="small"
            dataSource={history.map((s, i) => ({ ...s, key: i }))}
            columns={sessionColumns}
            pagination={false}
            expandable={{
              expandedRowRender: (record) => (
                <Table
                  size="small"
                  dataSource={record.entries}
                  columns={entryColumns}
                  rowKey="id"
                  pagination={{ defaultPageSize: 25, showSizeChanger: true }}
                />
              ),
            }}
          />
        </>
      )}
    </div>
  );
}
