import { useState } from "react";

const LEGEND_LIMIT = 20;

export default function CollapsibleLegend({ categories, maxHeight }) {
  const [expanded, setExpanded] = useState(false);
  const hasMore = categories.length > LEGEND_LIMIT;
  const visible = expanded ? categories : categories.slice(0, LEGEND_LIMIT);

  return (
    <div style={{ maxHeight, overflow: "auto", fontSize: 12 }}>
      {visible.map(([cat, color]) => (
        <div key={cat} style={{ display: "flex", alignItems: "center", gap: 4, marginBottom: 2 }}>
          <div
            style={{
              width: 12,
              height: 12,
              borderRadius: 2,
              backgroundColor: `rgb(${color.join(",")})`,
              flexShrink: 0,
            }}
          />
          <span>{cat}</span>
        </div>
      ))}
      {hasMore && (
        <span
          onClick={() => setExpanded(!expanded)}
          style={{ color: "#1890ff", cursor: "pointer", fontSize: 11 }}
        >
          {expanded ? "Show less" : `Show all (${categories.length})`}
        </span>
      )}
    </div>
  );
}
