import { useState, useMemo } from "react";
import { Table, Tag, Button, Card, Space, Typography, Popconfirm, Statistic, Row, Col, Input } from "antd";
import { DeleteOutlined, ClearOutlined, SearchOutlined } from "@ant-design/icons";
import {
  getProfileHistory,
  removeProfileSession,
  clearProfileHistory,
} from "../utils/profileStorage";
import MethodBreakdownChart from "../components/charts/MethodBreakdownChart";
import CacheEfficiencyChart from "../components/charts/CacheEfficiencyChart";
import BytesByMethodChart from "../components/charts/BytesByMethodChart";
import RequestsByMethodChart from "../components/charts/RequestsByMethodChart";
import BytesVsDurationChart from "../components/charts/BytesVsDurationChart";
import SessionWaterfallChart from "../components/charts/SessionWaterfallChart";

const { Text, Title } = Typography;

function formatBytes(bytes) {
  if (bytes === 0) return "0 B";
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

function formatShape(shape) {
  if (!shape || shape.length === 0) return "-";
  return shape.map((d) => d.toLocaleString()).join(" × ");
}

/** Build per-method I/O summary rows from a list of profile entries. */
function buildMethodSummary(entries) {
  const byMethod = {};
  for (const e of entries) {
    if (!byMethod[e.method]) {
      byMethod[e.method] = {
        method: e.method,
        calls: 0,
        cacheMisses: 0,
        totalRequests: 0,
        totalBytes: 0,
        totalDuration: 0,
        // track chunk info from first entry that has it
        arrayShape: null,
        chunkShape: null,
        dtype: null,
        sharded: null,
      };
    }
    const m = byMethod[e.method];
    m.calls++;
    if (!e.cacheHit) m.cacheMisses++;
    m.totalRequests += e.fetches?.requests ?? 0;
    m.totalBytes += e.fetches?.bytes ?? 0;
    m.totalDuration += e.duration;
    if (e.chunks && !m.arrayShape) {
      m.arrayShape = e.chunks.arrayShape;
      m.chunkShape = e.chunks.chunkShape;
      m.dtype = e.chunks.dtype;
      m.sharded = e.chunks.sharded;
    }
  }
  return Object.values(byMethod);
}

const methodSummaryColumns = [
  { title: "Method", dataIndex: "method", key: "method", width: 130 },
  {
    title: "Storage",
    key: "storage",
    width: 90,
    render: (_, r) =>
      r.sharded != null ? (
        <Tag color={r.sharded ? "purple" : "blue"}>{r.sharded ? "Sharded" : "Chunked"}</Tag>
      ) : (
        <Text type="secondary">-</Text>
      ),
  },
  {
    title: "Array Shape",
    key: "arrayShape",
    width: 140,
    render: (_, r) => r.arrayShape ? formatShape(r.arrayShape) : <Text type="secondary">-</Text>,
  },
  {
    title: "Chunk Shape",
    key: "chunkShape",
    width: 140,
    render: (_, r) => r.chunkShape ? formatShape(r.chunkShape) : <Text type="secondary">-</Text>,
  },
  {
    title: "dtype",
    key: "dtype",
    width: 80,
    render: (_, r) => r.dtype || <Text type="secondary">-</Text>,
  },
  {
    title: "Calls",
    dataIndex: "calls",
    key: "calls",
    width: 70,
    sorter: (a, b) => a.calls - b.calls,
  },
  {
    title: "Fetches",
    dataIndex: "totalRequests",
    key: "totalRequests",
    width: 80,
    sorter: (a, b) => a.totalRequests - b.totalRequests,
  },
  {
    title: "Avg / Call",
    key: "avgReqPerCall",
    width: 90,
    render: (_, r) => r.cacheMisses > 0 ? (r.totalRequests / r.cacheMisses).toFixed(1) : "-",
  },
  {
    title: "Total Bytes",
    dataIndex: "totalBytes",
    key: "totalBytes",
    width: 100,
    render: (b) => formatBytes(b),
    sorter: (a, b) => a.totalBytes - b.totalBytes,
  },
  {
    title: "Avg Bytes / Fetch",
    key: "avgBytesPerReq",
    width: 120,
    render: (_, r) => r.totalRequests > 0 ? formatBytes(Math.round(r.totalBytes / r.totalRequests)) : "-",
  },
  {
    title: "Total Time",
    dataIndex: "totalDuration",
    key: "totalDuration",
    width: 100,
    render: (ms) => `${ms.toFixed(1)} ms`,
    sorter: (a, b) => a.totalDuration - b.totalDuration,
  },
];

function keySearchFilter() {
  return {
    filterDropdown: ({ setSelectedKeys, selectedKeys, confirm, clearFilters }) => (
      <div style={{ padding: 8 }}>
        <Input
          placeholder="Search key…"
          value={selectedKeys[0]}
          onChange={(e) => setSelectedKeys(e.target.value ? [e.target.value] : [])}
          onPressEnter={confirm}
          style={{ width: 200, marginBottom: 8, display: "block" }}
          size="small"
        />
        <Space>
          <Button type="primary" onClick={confirm} icon={<SearchOutlined />} size="small" style={{ width: 90 }}>
            Search
          </Button>
          <Button onClick={() => { clearFilters?.(); confirm(); }} size="small" style={{ width: 90 }}>
            Reset
          </Button>
        </Space>
      </div>
    ),
    filterIcon: (filtered) => <SearchOutlined style={{ color: filtered ? "#1677ff" : undefined }} />,
    onFilter: (value, record) => {
      const v = value.toLowerCase();
      return record.key.toLowerCase().includes(v) || (record.label?.toLowerCase().includes(v) ?? false);
    },
  };
}

function buildEntryColumns(entries) {
  const methods = [...new Set(entries.map((e) => e.method))].sort();
  return [
    { title: "#", dataIndex: "id", key: "id", width: 50 },
    {
      title: "Method",
      dataIndex: "method",
      key: "method",
      width: 120,
      filters: methods.map((m) => ({ text: m, value: m })),
      onFilter: (value, record) => record.method === value,
    },
    {
      title: "Key",
      dataIndex: "key",
      key: "key",
      ellipsis: true,
      ...keySearchFilter(),
      render: (key, record) =>
        record.label ? (
          <span>
            {key} <Text type="secondary" style={{ fontSize: 11 }}>({record.label})</Text>
          </span>
        ) : (
          key
        ),
    },
    {
      title: "Status",
      key: "status",
      width: 90,
      filters: [
        { text: "HIT", value: "hit" },
        { text: "MISS", value: "miss" },
        { text: "ABORTED", value: "aborted" },
      ],
      onFilter: (value, record) => {
        if (value === "aborted") return !!record.aborted;
        if (value === "hit") return record.cacheHit && !record.aborted;
        return !record.cacheHit && !record.aborted;
      },
      render: (_, record) => {
        if (record.aborted) return <Tag color="orange">ABORTED</Tag>;
        return <Tag color={record.cacheHit ? "green" : "red"}>{record.cacheHit ? "HIT" : "MISS"}</Tag>;
      },
    },
    {
      title: "Storage",
      key: "storage",
      width: 90,
      filters: [
        { text: "Sharded", value: "sharded" },
        { text: "Chunked", value: "chunked" },
      ],
      onFilter: (value, record) => {
        if (!record.chunks) return false;
        return value === "sharded" ? record.chunks.sharded : !record.chunks.sharded;
      },
      render: (_, r) =>
        r.chunks ? (
          <Tag color={r.chunks.sharded ? "purple" : "blue"}>
            {r.chunks.sharded ? "Shard" : "Chunk"}
          </Tag>
        ) : null,
    },
    {
      title: "Chunk Shape",
      key: "chunkShape",
      width: 130,
      render: (_, r) => r.chunks ? formatShape(r.chunks.chunkShape) : <Text type="secondary">-</Text>,
    },
    {
      title: "Requests",
      key: "requests",
      width: 80,
      render: (_, r) => r.fetches?.requests ?? <Text type="secondary">-</Text>,
      sorter: (a, b) => (a.fetches?.requests ?? 0) - (b.fetches?.requests ?? 0),
    },
    {
      title: "Bytes",
      key: "bytes",
      width: 90,
      render: (_, r) => r.fetches ? formatBytes(r.fetches.bytes) : <Text type="secondary">-</Text>,
      sorter: (a, b) => (a.fetches?.bytes ?? 0) - (b.fetches?.bytes ?? 0),
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
}

function SummaryStats({ history }) {
  const stats = useMemo(() => {
    let totalQueries = 0;
    let totalHits = 0;
    let slowestQuery = 0;
    let totalBytes = 0;

    for (const session of history) {
      for (const e of session.entries) {
        totalQueries++;
        if (e.cacheHit) totalHits++;
        if (e.duration > slowestQuery) slowestQuery = e.duration;
        totalBytes += e.fetches?.bytes ?? 0;
      }
    }

    const avgHitRate = totalQueries > 0 ? ((totalHits / totalQueries) * 100).toFixed(1) : 0;

    return { totalSessions: history.length, totalQueries, avgHitRate, slowestQuery, totalBytes };
  }, [history]);

  return (
    <Row gutter={16} style={{ marginBottom: 24 }}>
      <Col span={5}>
        <Card size="small">
          <Statistic title="Total Sessions" value={stats.totalSessions} />
        </Card>
      </Col>
      <Col span={5}>
        <Card size="small">
          <Statistic title="Total Queries" value={stats.totalQueries} />
        </Card>
      </Col>
      <Col span={5}>
        <Card size="small">
          <Statistic title="Avg Cache Hit Rate" value={stats.avgHitRate} suffix="%" />
        </Card>
      </Col>
      <Col span={5}>
        <Card size="small">
          <Statistic title="Slowest Query" value={stats.slowestQuery.toFixed(1)} suffix="ms" />
        </Card>
      </Col>
      <Col span={4}>
        <Card size="small">
          <Statistic title="Total Transfer" value={formatBytes(stats.totalBytes)} />
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

function SessionDetail({ record }) {
  const methodRows = useMemo(() => buildMethodSummary(record.entries), [record.entries]);
  const entryColumns = useMemo(() => buildEntryColumns(record.entries), [record.entries]);

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 16 }}>
      {/* Row 1: three charts side by side */}
      <Row gutter={16}>
        <Col span={8}>
          <Card size="small" title="Bytes by Method">
            <BytesByMethodChart entries={record.entries} width={320} />
          </Card>
        </Col>
        <Col span={8}>
          <Card size="small" title="Requests by Method">
            <RequestsByMethodChart entries={record.entries} width={320} />
          </Card>
        </Col>
        <Col span={8}>
          <Card size="small" title="Bytes vs Duration">
            <BytesVsDurationChart entries={record.entries} width={320} />
          </Card>
        </Col>
      </Row>

      {/* Row 2: full-width waterfall */}
      <Card size="small" title="Session Waterfall Timeline">
        <SessionWaterfallChart entries={record.entries} width={1000} />
      </Card>

      <div>
        <Text strong style={{ fontSize: 12, display: "block", marginBottom: 8 }}>
          Method I/O Summary
        </Text>
        <Table
          size="small"
          dataSource={methodRows}
          columns={methodSummaryColumns}
          rowKey="method"
          pagination={false}
        />
      </div>
      <div>
        <Text strong style={{ fontSize: 12, display: "block", marginBottom: 8 }}>
          All Queries
        </Text>
        <Table
          size="small"
          dataSource={record.entries}
          columns={entryColumns}
          rowKey="id"
          pagination={{ defaultPageSize: 25, showSizeChanger: true }}
        />
      </div>
    </div>
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
      title: "Host",
      dataIndex: "url",
      key: "host",
      width: 140,
      render: (url) => {
        try { return new URL(url).host; } catch { return "-"; }
      },
    },
    {
      title: "Dataset",
      dataIndex: "url",
      key: "dataset",
      render: (url) => {
        let name = url;
        try {
          name = decodeURIComponent(new URL(url).pathname).replace(/\/$/, "").split("/").pop() || url;
        } catch { /* use raw url */ }
        return name;
      },
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
      title: "Transfer",
      key: "totalBytes",
      width: 100,
      render: (_, r) => formatBytes(r.entries.reduce((s, e) => s + (e.fetches?.bytes ?? 0), 0)),
      sorter: (a, b) =>
        a.entries.reduce((s, e) => s + (e.fetches?.bytes ?? 0), 0) -
        b.entries.reduce((s, e) => s + (e.fetches?.bytes ?? 0), 0),
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
              expandedRowRender: (record) => <SessionDetail record={record} />,
            }}
          />
        </>
      )}
    </div>
  );
}
