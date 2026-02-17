import { useState } from "react";
import { Typography } from "antd";

const { Text } = Typography;
const BREAKDOWN_LIMIT = 5;

function CollapsibleBreakdown({ col, breakdown, total, onHoverValue }) {
  const [expanded, setExpanded] = useState(false);
  const hasMore = breakdown.length > BREAKDOWN_LIMIT;
  const visible = expanded ? breakdown : breakdown.slice(0, BREAKDOWN_LIMIT);

  return (
    <div style={{ marginTop: 8 }}>
      <Text type="secondary" style={{ fontSize: 11 }}>{col}</Text>
      <table style={{ width: "100%", marginTop: 4, borderCollapse: "collapse" }}>
        <tbody>
          {visible.map(([val, count]) => (
            <tr
              key={val}
              style={{ cursor: "default" }}
              onMouseEnter={() => onHoverValue?.({ col, value: val })}
              onMouseLeave={() => onHoverValue?.(null)}
            >
              <td style={{ paddingRight: 8 }}>{val}</td>
              <td style={{ textAlign: "right", whiteSpace: "nowrap" }}>
                {count.toLocaleString()} ({((count / total) * 100).toFixed(1)}%)
              </td>
            </tr>
          ))}
        </tbody>
      </table>
      {hasMore && (
        <span
          onClick={() => setExpanded(!expanded)}
          style={{ color: "#1890ff", cursor: "pointer", fontSize: 11 }}
        >
          {expanded ? "Show less" : `Show all (${breakdown.length})`}
        </span>
      )}
    </div>
  );
}

export default function SelectionSummaryPanel({
  selectionSummary,
  selectedCount,
  categoryColorMap,
  colorColumn,
  selectedGene,
  maxHeight,
  onHoverCategory,
  onHoverTooltipValue,
}) {
  return (
    <div style={{ maxHeight, overflow: "auto", fontSize: 12, minWidth: 160, borderLeft: "1px solid #d9d9d9", paddingLeft: 16 }}>
      <Text strong style={{ fontSize: 12 }}>
        Selection ({selectedCount.toLocaleString()} points)
      </Text>

      {selectionSummary.categoryBreakdown && (
        <div style={{ marginTop: 8 }}>
          <Text type="secondary" style={{ fontSize: 11 }}>{colorColumn}</Text>
          <table style={{ width: "100%", marginTop: 4, borderCollapse: "collapse" }}>
            <tbody>
              {selectionSummary.categoryBreakdown.map(([cat, count]) => (
                <tr
                  key={cat}
                  style={{ cursor: "default" }}
                  onMouseEnter={() => onHoverCategory?.(cat)}
                  onMouseLeave={() => onHoverCategory?.(null)}
                >
                  <td style={{ paddingRight: 8, display: "flex", alignItems: "center", gap: 4 }}>
                    {categoryColorMap[cat] && (
                      <span style={{
                        display: "inline-block",
                        width: 8,
                        height: 8,
                        borderRadius: 2,
                        backgroundColor: `rgb(${categoryColorMap[cat].join(",")})`,
                        flexShrink: 0,
                      }} />
                    )}
                    <span>{cat}</span>
                  </td>
                  <td style={{ textAlign: "right", whiteSpace: "nowrap" }}>
                    {count.toLocaleString()} ({((count / selectedCount) * 100).toFixed(1)}%)
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      {selectionSummary.expressionStats && (
        <div style={{ marginTop: 8 }}>
          <Text type="secondary" style={{ fontSize: 11 }}>{selectedGene} expression</Text>
          <table style={{ width: "100%", marginTop: 4, borderCollapse: "collapse" }}>
            <tbody>
              <tr><td>Mean</td><td style={{ textAlign: "right" }}>{selectionSummary.expressionStats.mean.toFixed(4)}</td></tr>
              <tr><td>Min</td><td style={{ textAlign: "right" }}>{selectionSummary.expressionStats.min.toFixed(4)}</td></tr>
              <tr><td>Max</td><td style={{ textAlign: "right" }}>{selectionSummary.expressionStats.max.toFixed(4)}</td></tr>
            </tbody>
          </table>
        </div>
      )}

      {Object.entries(selectionSummary.tooltipBreakdowns).map(([col, breakdown]) => (
        <CollapsibleBreakdown key={col} col={col} breakdown={breakdown} total={selectedCount} onHoverValue={onHoverTooltipValue} />
      ))}
    </div>
  );
}
