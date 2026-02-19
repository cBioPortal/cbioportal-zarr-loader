import { useCallback } from "react";
import { Group } from "@visx/group";
import { scaleBand, scaleLinear, scaleSqrt } from "@visx/scale";
import { AxisBottom, AxisLeft } from "@visx/axis";
import { useTooltip, TooltipWithBounds, defaultStyles } from "@visx/tooltip";
import { COLOR_SCALES, colorScaleGradient, interpolateColorScale, rgbToString } from "../../utils/colors";

const MARGIN = { top: 16, right: 100, bottom: 40, left: 120 };
const MAX_RADIUS = 14;

const tooltipStyles = {
  ...defaultStyles,
  fontSize: 12,
  padding: "6px 10px",
};

export default function Dotplot({ genes, groups, data, width = 600, height = 400, showLabels = false }) {
  const { showTooltip, hideTooltip, tooltipOpen, tooltipData, tooltipLeft, tooltipTop } =
    useTooltip();

  const xMax = width - MARGIN.left - MARGIN.right;
  const yMax = height - MARGIN.top - MARGIN.bottom;

  // Map groups to integer labels for compact x-axis
  const groupIndices = groups.map((_, i) => i + 1);
  const groupToIndex = Object.fromEntries(groups.map((g, i) => [g, i + 1]));
  const indexToGroup = Object.fromEntries(groups.map((g, i) => [i + 1, g]));

  const xDomain = showLabels ? groups : groupIndices;
  const xScale = scaleBand({ domain: xDomain, range: [0, xMax], padding: 0.05 });
  const yScale = scaleBand({ domain: genes, range: [0, yMax], padding: 0.05 });

  const maxMean = Math.max(...data.map((d) => d.meanExpression), 0.01);

  const maxR = Math.min(MAX_RADIUS, xScale.bandwidth() / 2, yScale.bandwidth() / 2);
  const radiusScale = scaleSqrt({
    domain: [0, 1],
    range: [0, maxR],
  });

  const viridis = COLOR_SCALES.viridis;
  const colorScale = (val) => rgbToString(interpolateColorScale(val / maxMean, viridis));

  const handleMouseEnter = useCallback(
    (event, d) => {
      const svg = event.currentTarget.ownerSVGElement;
      const point = svg.createSVGPoint();
      point.x = event.clientX;
      point.y = event.clientY;
      const svgPoint = point.matrixTransform(svg.getScreenCTM().inverse());
      showTooltip({
        tooltipData: d,
        tooltipLeft: svgPoint.x,
        tooltipTop: svgPoint.y - 10,
      });
    },
    [showTooltip],
  );

  return (
    <div style={{ position: "relative" }}>
      <svg width={width} height={height}>
        <Group left={MARGIN.left} top={MARGIN.top}>
          {/* Horizontal gridlines */}
          {genes.map((gene) => (
            <line
              key={`h-${gene}`}
              x1={0}
              x2={xMax}
              y1={(yScale(gene) ?? 0) + yScale.bandwidth() / 2}
              y2={(yScale(gene) ?? 0) + yScale.bandwidth() / 2}
              stroke="#f0f0f0"
              strokeWidth={1}
            />
          ))}
          {/* Vertical gridlines */}
          {xDomain.map((val) => (
            <line
              key={`v-${val}`}
              x1={(xScale(val) ?? 0) + xScale.bandwidth() / 2}
              x2={(xScale(val) ?? 0) + xScale.bandwidth() / 2}
              y1={0}
              y2={yMax}
              stroke="#f0f0f0"
              strokeWidth={1}
            />
          ))}
          {data.map((d) => {
            const xKey = showLabels ? d.group : groupToIndex[d.group];
            const cx = (xScale(xKey) ?? 0) + xScale.bandwidth() / 2;
            const cy = (yScale(d.gene) ?? 0) + yScale.bandwidth() / 2;
            const isEmpty = d.fractionExpressing === 0;
            const r = isEmpty ? 3 : Math.max(radiusScale(d.fractionExpressing), 4);
            return (
              <circle
                key={`${d.gene}-${d.group}`}
                cx={cx}
                cy={cy}
                r={r}
                fill={isEmpty ? "#f0f0f0" : colorScale(d.meanExpression)}
                stroke={isEmpty ? "#ccc" : colorScale(d.meanExpression)}
                strokeWidth={isEmpty ? 1 : 0.5}
                style={{ cursor: "pointer" }}
                onMouseEnter={(e) => handleMouseEnter(e, d)}
                onMouseLeave={hideTooltip}
              />
            );
          })}
          <AxisBottom
            scale={xScale}
            top={yMax}
            tickComponent={({ x, y, formattedValue }) =>
              showLabels ? (
                <text
                  x={x}
                  y={y}
                  fontSize={11}
                  textAnchor="end"
                  dy={-4}
                  transform={`rotate(-45, ${x}, ${y})`}
                >
                  {formattedValue}
                </text>
              ) : (
                <text
                  x={x}
                  y={y}
                  fontSize={11}
                  textAnchor="middle"
                  dy="0.25em"
                  style={{ cursor: "pointer" }}
                  onMouseEnter={(e) => {
                    const svg = e.currentTarget.ownerSVGElement;
                    const pt = svg.createSVGPoint();
                    pt.x = e.clientX;
                    pt.y = e.clientY;
                    const svgPt = pt.matrixTransform(svg.getScreenCTM().inverse());
                    showTooltip({
                      tooltipData: { axisLabel: indexToGroup[formattedValue] },
                      tooltipLeft: svgPt.x,
                      tooltipTop: svgPt.y - 10,
                    });
                  }}
                  onMouseLeave={hideTooltip}
                >
                  {formattedValue}
                </text>
              )
            }
          />
          <AxisLeft
            scale={yScale}
            numTicks={genes.length}
            tickLabelProps={() => ({
              fontSize: 11,
              textAnchor: "end",
              dx: -4,
              dy: 4,
            })}
          />
          {/* Expression color legend */}
          <Group left={xMax + 16} top={0}>
            <text fontSize={10} fontWeight="bold" dy={-4}>Mean expr</text>
            <defs>
              <linearGradient id="dotplot-color-gradient" x1="0" y1="1" x2="0" y2="0">
                {viridis.map((c, i) => (
                  <stop
                    key={i}
                    offset={`${(i / (viridis.length - 1)) * 100}%`}
                    stopColor={rgbToString(c)}
                  />
                ))}
              </linearGradient>
            </defs>
            <rect
              width={14}
              height={Math.min(yMax, 100)}
              fill="url(#dotplot-color-gradient)"
              rx={2}
            />
            <text x={18} y={10} fontSize={9}>{maxMean.toFixed(2)}</text>
            <text x={18} y={Math.min(yMax, 100)} fontSize={9}>0</text>
            {/* Size legend */}
            <text fontSize={10} fontWeight="bold" y={Math.min(yMax, 100) + 24}>Fraction</text>
            {[0.25, 0.5, 0.75, 1.0].map((frac, i) => {
              const r = radiusScale(frac);
              const cy = Math.min(yMax, 100) + 46 + i * (maxR * 2 + 4);
              return (
                <g key={frac}>
                  <circle cx={maxR} cy={cy} r={r} fill="#888" />
                  <text x={maxR * 2 + 6} y={cy + 3} fontSize={9}>{(frac * 100)}%</text>
                </g>
              );
            })}
          </Group>
        </Group>
      </svg>
      {tooltipOpen && tooltipData && (
        <TooltipWithBounds
          left={tooltipLeft}
          top={tooltipTop}
          style={tooltipStyles}
        >
          {tooltipData.axisLabel ? (
            <div>{tooltipData.axisLabel}</div>
          ) : (
            <>
              <div><strong>{tooltipData.group}</strong></div>
              <div>{tooltipData.gene}</div>
              <div>Cells: {tooltipData.expressingCount.toLocaleString()} / {tooltipData.cellCount.toLocaleString()}</div>
              <div>Fraction: {(tooltipData.fractionExpressing * 100).toFixed(1)}%</div>
              <div>Mean expr: {tooltipData.meanExpression.toFixed(3)}</div>
            </>
          )}
        </TooltipWithBounds>
      )}
    </div>
  );
}
