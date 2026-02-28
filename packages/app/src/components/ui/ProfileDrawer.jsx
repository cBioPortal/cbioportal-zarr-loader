import { useSyncExternalStore, useCallback, useRef, useEffect, useState, useMemo } from "react";
import { Drawer, Table, Tag, Space, Button, Typography, message, Collapse } from "antd";
import { DeleteOutlined, SaveOutlined } from "@ant-design/icons";
import { Group } from "@visx/group";
import { scaleLinear } from "@visx/scale";
import { AxisBottom } from "@visx/axis";
import { useTooltip, TooltipWithBounds, defaultStyles } from "@visx/tooltip";
import useAppStore from "../../store/useAppStore";
import { saveProfileSession } from "../../utils/profileStorage";

const { Text } = Typography;

// Method → color mapping (consistent across waterfall and charts)
const METHOD_COLORS = {
  obsm: "#1f77b4",
  obs: "#ff7f0e",
  geneExpression: "#2ca02c",
  var: "#d62728",
  obsColumns: "#9467bd",
  X: "#8c564b",
  uns: "#e377c2",
};
const DEFAULT_METHOD_COLOR = "#7f7f7f";

function getMethodColor(method) {
  return METHOD_COLORS[method] || DEFAULT_METHOD_COLOR;
}

const BAR_HEIGHT = 18;
const BAR_GAP = 3;
const MARGIN = { top: 10, right: 16, bottom: 30, left: 16 };
const WATERFALL_MAX_VISIBLE = 60;

const tooltipStyles = {
  ...defaultStyles,
  fontSize: 12,
  padding: "6px 10px",
};

function WaterfallTimeline({ entries, width }) {
  const { showTooltip, hideTooltip, tooltipOpen, tooltipData, tooltipLeft, tooltipTop } =
    useTooltip();
  const scrollRef = useRef(null);

  // Auto-scroll to bottom when new entries arrive
  useEffect(() => {
    if (scrollRef.current) {
      scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
    }
  }, [entries.length]);

  if (entries.length === 0) {
    return (
      <div style={{ padding: "16px 0", textAlign: "center", color: "#999" }}>
        <Text type="secondary">No queries recorded yet. Interact with the viewer to see profiling data.</Text>
      </div>
    );
  }

  const innerWidth = width - MARGIN.left - MARGIN.right;
  const svgHeight = entries.length * (BAR_HEIGHT + BAR_GAP) + MARGIN.top + MARGIN.bottom;

  // Time axis: origin at the earliest entry's startTime
  const t0 = entries.length > 0 ? entries[0].startTime : 0;
  const tMax = entries.reduce((max, e) => Math.max(max, e.startTime + e.duration - t0), 1);

  const xScale = scaleLinear({
    domain: [0, tMax],
    range: [0, innerWidth],
  });

  return (
    <div style={{ position: "relative" }}>
      <div
        ref={scrollRef}
        style={{
          maxHeight: WATERFALL_MAX_VISIBLE * (BAR_HEIGHT + BAR_GAP) + MARGIN.top + MARGIN.bottom,
          overflowY: "auto",
          overflowX: "hidden",
        }}
      >
        <svg width={width} height={svgHeight}>
          <Group left={MARGIN.left} top={MARGIN.top}>
            {entries.map((entry, i) => {
              const x = xScale(entry.startTime - t0);
              const barWidth = Math.max(xScale(entry.duration) - xScale(0), 2);
              const y = i * (BAR_HEIGHT + BAR_GAP);
              const color = getMethodColor(entry.method);

              return (
                <g key={entry.id}>
                  <rect
                    x={x}
                    y={y}
                    width={barWidth}
                    height={BAR_HEIGHT}
                    fill={entry.cacheHit ? "transparent" : color}
                    stroke={color}
                    strokeWidth={entry.cacheHit ? 1.5 : 0}
                    rx={2}
                    style={{ cursor: "pointer" }}
                    onMouseEnter={(e) => {
                      const svg = e.currentTarget.ownerSVGElement;
                      const point = svg.createSVGPoint();
                      point.x = e.clientX;
                      point.y = e.clientY;
                      const svgPoint = point.matrixTransform(svg.getScreenCTM().inverse());
                      showTooltip({
                        tooltipData: entry,
                        tooltipLeft: svgPoint.x,
                        tooltipTop: svgPoint.y - 10,
                      });
                    }}
                    onMouseLeave={hideTooltip}
                  />
                  {/* Method label inside bar if wide enough */}
                  {barWidth > 40 && (
                    <text
                      x={x + 4}
                      y={y + BAR_HEIGHT / 2}
                      fontSize={10}
                      fill={entry.cacheHit ? color : "#fff"}
                      dominantBaseline="central"
                      pointerEvents="none"
                    >
                      {entry.method}
                    </text>
                  )}
                </g>
              );
            })}
            <AxisBottom
              scale={xScale}
              top={entries.length * (BAR_HEIGHT + BAR_GAP)}
              numTicks={5}
              tickFormat={(v) => `${v.toFixed(0)} ms`}
              tickLabelProps={() => ({ fontSize: 10, textAnchor: "middle", fill: "#999" })}
              stroke="#ddd"
              tickStroke="#ddd"
            />
          </Group>
        </svg>
      </div>

      {tooltipOpen && tooltipData && (
        <TooltipWithBounds left={tooltipLeft} top={tooltipTop} style={tooltipStyles}>
          <div><strong>{tooltipData.method}</strong></div>
          <div style={{ fontSize: 11, color: "#666" }}>{tooltipData.key}</div>
          <div>{tooltipData.duration.toFixed(1)} ms</div>
          <div>
            <Tag
              color={tooltipData.cacheHit ? "green" : "red"}
              style={{ margin: 0, marginTop: 2 }}
            >
              {tooltipData.cacheHit ? "CACHE HIT" : "CACHE MISS"}
            </Tag>
          </div>
        </TooltipWithBounds>
      )}

      {/* Legend */}
      <div style={{ display: "flex", flexWrap: "wrap", gap: "8px 12px", padding: "8px 0 0", fontSize: 11 }}>
        {Object.entries(METHOD_COLORS).map(([method, color]) => (
          <span key={method} style={{ display: "flex", alignItems: "center", gap: 4 }}>
            <span style={{ width: 10, height: 10, borderRadius: 2, background: color, display: "inline-block" }} />
            {method}
          </span>
        ))}
        <span style={{ display: "flex", alignItems: "center", gap: 4 }}>
          <span style={{ width: 10, height: 10, borderRadius: 2, border: "1.5px solid #999", display: "inline-block" }} />
          cache hit (outline)
        </span>
      </div>
    </div>
  );
}

const tableColumns = [
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

// CSS for the pulse animation (injected once)
const pulseStyle = `
@keyframes czl-pulse {
  0% { opacity: 1; transform: scale(1); }
  50% { opacity: 0.4; transform: scale(1.3); }
  100% { opacity: 0; transform: scale(1); }
}
`;

export default function ProfileDrawer({ open, onClose }) {
  const adata = useAppStore((s) => s.adata);
  const url = useAppStore((s) => s.url);
  const profiler = adata?.profiler;

  const subscribe = useCallback(
    (cb) => (profiler ? profiler.subscribe(cb) : () => {}),
    [profiler],
  );
  const getSnapshot = useCallback(
    () => (profiler ? profiler.version : 0),
    [profiler],
  );

  const version = useSyncExternalStore(subscribe, getSnapshot);
  const entries = profiler?.entries ?? [];

  // Activity pulse: show for 1s after version changes
  const [pulseVisible, setPulseVisible] = useState(false);
  const prevVersionRef = useRef(version);
  useEffect(() => {
    if (version !== prevVersionRef.current && version > 0) {
      setPulseVisible(true);
      const timer = setTimeout(() => setPulseVisible(false), 1000);
      prevVersionRef.current = version;
      return () => clearTimeout(timer);
    }
    prevVersionRef.current = version;
  }, [version]);

  const totalTime = entries.reduce((sum, e) => sum + e.duration, 0);
  const cacheHits = entries.filter((e) => e.cacheHit).length;
  const hitRate = entries.length > 0 ? ((cacheHits / entries.length) * 100).toFixed(0) : 0;

  const handleClear = () => profiler?.clear();

  const handleSave = () => {
    if (!adata || entries.length === 0) return;
    saveProfileSession(url, adata.nObs, adata.nVar, profiler.toJSON());
    message.success("Profile session saved to history");
  };

  const dataSource = useMemo(() => [...entries].reverse(), [entries, version]);

  // Drawer body width estimate (drawer width minus padding)
  const waterfallWidth = 600 - 48;

  return (
    <Drawer
      title={
        <span style={{ display: "flex", alignItems: "center", gap: 8 }}>
          <style>{pulseStyle}</style>
          Query Profiler
          {pulseVisible && (
            <span
              style={{
                width: 8,
                height: 8,
                borderRadius: "50%",
                background: "#52c41a",
                display: "inline-block",
                animation: "czl-pulse 1s ease-out forwards",
              }}
            />
          )}
        </span>
      }
      placement="right"
      width={600}
      open={open}
      onClose={onClose}
      mask={false}
    >
      <Space style={{ marginBottom: 12 }} wrap>
        <Text strong>Queries: {entries.length}</Text>
        <Text type="secondary">|</Text>
        <Text>Cache hit: {hitRate}%</Text>
        <Text type="secondary">|</Text>
        <Text>Total: {totalTime.toFixed(1)} ms</Text>
        <Button icon={<SaveOutlined />} size="small" onClick={handleSave} disabled={entries.length === 0}>
          Save to history
        </Button>
        <Button icon={<DeleteOutlined />} size="small" danger onClick={handleClear} disabled={entries.length === 0}>
          Clear
        </Button>
      </Space>

      <WaterfallTimeline entries={entries} width={waterfallWidth} />

      <Collapse
        ghost
        items={[
          {
            key: "log",
            label: <Text type="secondary" style={{ fontSize: 12 }}>Query log table</Text>,
            children: (
              <Table
                size="small"
                pagination={{ defaultPageSize: 50, showSizeChanger: true }}
                dataSource={dataSource}
                columns={tableColumns}
                rowKey="id"
              />
            ),
          },
        ]}
        style={{ marginTop: 12 }}
      />
    </Drawer>
  );
}
