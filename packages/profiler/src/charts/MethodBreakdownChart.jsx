import { Group } from "@visx/group";
import { scaleBand, scaleLinear, scaleOrdinal } from "@visx/scale";
import { AxisBottom, AxisLeft } from "@visx/axis";
import { useTooltip, TooltipWithBounds, defaultStyles } from "@visx/tooltip";
import { BarStackHorizontal } from "@visx/shape";
import { METHOD_COLORS, DEFAULT_METHOD_COLOR } from "../constants";

const tooltipStyles = {
  ...defaultStyles,
  fontSize: 12,
  padding: "6px 10px",
};

/**
 * Horizontal stacked bar chart — total duration broken down by method per session.
 * @param {{ sessions: Array<{ label: string, entries: Array<{ method: string, duration: number }> }>, width?: number, height?: number }} props
 */
export default function MethodBreakdownChart({ sessions, width = 500, height: heightProp }) {
  const { showTooltip, hideTooltip, tooltipOpen, tooltipData, tooltipLeft, tooltipTop } =
    useTooltip();

  if (!sessions || sessions.length === 0) return null;

  // Aggregate duration per method per session
  const allMethods = new Set();
  const data = sessions.map((s, i) => {
    const byMethod = {};
    for (const e of s.entries) {
      byMethod[e.method] = (byMethod[e.method] || 0) + e.duration;
      allMethods.add(e.method);
    }
    return { label: s.label || `Session ${i + 1}`, ...byMethod };
  });

  const methods = [...allMethods];

  const MARGIN = { top: 10, right: 20, bottom: 40, left: 120 };
  const minBarHeight = 28;
  const height = heightProp || Math.max(200, sessions.length * minBarHeight + MARGIN.top + MARGIN.bottom + 20);
  const xMax = width - MARGIN.left - MARGIN.right;
  const yMax = height - MARGIN.top - MARGIN.bottom;

  const xScale = scaleLinear({
    domain: [0, Math.max(...data.map((d) => methods.reduce((sum, m) => sum + (d[m] || 0), 0)))],
    range: [0, xMax],
    nice: true,
  });

  const yScale = scaleBand({
    domain: data.map((d) => d.label),
    range: [0, yMax],
    padding: 0.25,
  });

  const colorScale = scaleOrdinal({
    domain: methods,
    range: methods.map((m) => METHOD_COLORS[m] || DEFAULT_METHOD_COLOR),
  });

  return (
    <div style={{ position: "relative" }}>
      <svg width={width} height={height}>
        <Group left={MARGIN.left} top={MARGIN.top}>
          <BarStackHorizontal
            data={data}
            keys={methods}
            y={(d) => d.label}
            xScale={xScale}
            yScale={yScale}
            color={colorScale}
          >
            {(barStacks) =>
              barStacks.map((barStack) =>
                barStack.bars.map((bar) => (
                  <rect
                    key={`${barStack.index}-${bar.index}`}
                    x={bar.x}
                    y={bar.y}
                    width={Math.max(0, bar.width)}
                    height={Math.max(0, bar.height)}
                    fill={bar.color}
                    rx={2}
                    style={{ cursor: "pointer" }}
                    onMouseEnter={(e) => {
                      const svg = e.currentTarget.ownerSVGElement;
                      const point = svg.createSVGPoint();
                      point.x = e.clientX;
                      point.y = e.clientY;
                      const svgPoint = point.matrixTransform(svg.getScreenCTM().inverse());
                      showTooltip({
                        tooltipData: { method: barStack.key, value: bar.bar.data[barStack.key], session: bar.bar.data.label },
                        tooltipLeft: svgPoint.x,
                        tooltipTop: svgPoint.y - 10,
                      });
                    }}
                    onMouseLeave={hideTooltip}
                  />
                )),
              )
            }
          </BarStackHorizontal>
          <AxisBottom
            scale={xScale}
            top={yMax}
            numTicks={5}
            tickFormat={(v) => `${v.toFixed(0)} ms`}
            tickLabelProps={() => ({ fontSize: 11, textAnchor: "middle", fill: "#666" })}
          />
          <AxisLeft
            scale={yScale}
            tickLabelProps={() => ({ fontSize: 11, textAnchor: "end", fill: "#666", dx: -4 })}
          />
        </Group>
      </svg>

      {/* Legend */}
      <div style={{ display: "flex", flexWrap: "wrap", gap: "6px 12px", padding: "4px 0 0", fontSize: 11 }}>
        {methods.map((m) => (
          <span key={m} style={{ display: "flex", alignItems: "center", gap: 4 }}>
            <span style={{ width: 10, height: 10, borderRadius: 2, background: colorScale(m), display: "inline-block" }} />
            {m}
          </span>
        ))}
      </div>

      {tooltipOpen && tooltipData && (
        <TooltipWithBounds left={tooltipLeft} top={tooltipTop} style={tooltipStyles}>
          <div><strong>{tooltipData.session}</strong></div>
          <div>{tooltipData.method}: {tooltipData.value?.toFixed(1)} ms</div>
        </TooltipWithBounds>
      )}
    </div>
  );
}
