import { Group } from "@visx/group";
import { scaleBand, scaleLinear, scaleOrdinal } from "@visx/scale";
import { AxisBottom, AxisLeft } from "@visx/axis";
import { useTooltip, TooltipWithBounds, defaultStyles } from "@visx/tooltip";

const METHOD_COLORS = {
  obsm: "#1f77b4",
  obs: "#ff7f0e",
  geneExpression: "#2ca02c",
  var: "#d62728",
  obsColumns: "#9467bd",
  X: "#8c564b",
  uns: "#e377c2",
};
const FALLBACK_COLOR = "#7f7f7f";

const tooltipStyles = {
  ...defaultStyles,
  fontSize: 12,
  padding: "6px 10px",
};

function formatBytes(bytes) {
  if (bytes === 0) return "0 B";
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

/**
 * Horizontal bar chart — total bytes transferred per method.
 * @param {{ entries: Array<{ method: string, fetches?: { bytes: number } }>, width?: number, height?: number }} props
 */
export default function BytesByMethodChart({ entries, width = 360, height: heightProp }) {
  const { showTooltip, hideTooltip, tooltipOpen, tooltipData, tooltipLeft, tooltipTop } =
    useTooltip();

  if (!entries || entries.length === 0) return null;

  // Aggregate bytes per method
  const byMethod = {};
  for (const e of entries) {
    const bytes = e.fetches?.bytes ?? 0;
    if (bytes > 0) {
      byMethod[e.method] = (byMethod[e.method] || 0) + bytes;
    }
  }

  const data = Object.entries(byMethod)
    .map(([method, bytes]) => ({ method, bytes }))
    .sort((a, b) => b.bytes - a.bytes);

  if (data.length === 0) return null;

  const MARGIN = { top: 10, right: 20, bottom: 40, left: 110 };
  const minBarHeight = 28;
  const height = heightProp || Math.max(160, data.length * minBarHeight + MARGIN.top + MARGIN.bottom + 20);
  const xMax = width - MARGIN.left - MARGIN.right;
  const yMax = height - MARGIN.top - MARGIN.bottom;

  const xScale = scaleLinear({
    domain: [0, Math.max(...data.map((d) => d.bytes))],
    range: [0, xMax],
    nice: true,
  });

  const yScale = scaleBand({
    domain: data.map((d) => d.method),
    range: [0, yMax],
    padding: 0.25,
  });

  const colorScale = scaleOrdinal({
    domain: data.map((d) => d.method),
    range: data.map((d) => METHOD_COLORS[d.method] || FALLBACK_COLOR),
  });

  return (
    <div style={{ position: "relative" }}>
      <svg width={width} height={height}>
        <Group left={MARGIN.left} top={MARGIN.top}>
          {data.map((d) => {
            const barWidth = xScale(d.bytes);
            const y = yScale(d.method);
            return (
              <rect
                key={d.method}
                x={0}
                y={y}
                width={barWidth}
                height={yScale.bandwidth()}
                fill={colorScale(d.method)}
                rx={2}
                style={{ cursor: "pointer" }}
                onMouseEnter={(e) => {
                  const svg = e.currentTarget.ownerSVGElement;
                  const point = svg.createSVGPoint();
                  point.x = e.clientX;
                  point.y = e.clientY;
                  const svgPoint = point.matrixTransform(svg.getScreenCTM().inverse());
                  showTooltip({
                    tooltipData: d,
                    tooltipLeft: svgPoint.x,
                    tooltipTop: svgPoint.y - 10,
                  });
                }}
                onMouseLeave={hideTooltip}
              />
            );
          })}
          <AxisBottom
            scale={xScale}
            top={yMax}
            numTicks={4}
            tickFormat={(v) => formatBytes(v)}
            tickLabelProps={() => ({ fontSize: 10, textAnchor: "middle", fill: "#666" })}
          />
          <AxisLeft
            scale={yScale}
            tickLabelProps={() => ({ fontSize: 11, textAnchor: "end", fill: "#666", dx: -4 })}
          />
        </Group>
      </svg>

      {tooltipOpen && tooltipData && (
        <TooltipWithBounds left={tooltipLeft} top={tooltipTop} style={tooltipStyles}>
          <div><strong>{tooltipData.method}</strong></div>
          <div>{formatBytes(tooltipData.bytes)}</div>
        </TooltipWithBounds>
      )}
    </div>
  );
}
