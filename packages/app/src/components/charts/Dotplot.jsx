import { useCallback } from "react";
import { Group } from "@visx/group";
import { scaleBand, scaleLinear, scaleSqrt } from "@visx/scale";
import { AxisBottom, AxisLeft } from "@visx/axis";
import { useTooltip, TooltipWithBounds, defaultStyles } from "@visx/tooltip";
import { COLOR_SCALES, colorScaleGradient, interpolateColorScale, rgbToString } from "../../utils/colors";

const MAX_RADIUS = 14;

const tooltipStyles = {
  ...defaultStyles,
  fontSize: 12,
  padding: "6px 10px",
};

export default function Dotplot({ genes, groups, data, width = 600, height = 400, showLabels = false, swapAxes = false }) {
  const { showTooltip, hideTooltip, tooltipOpen, tooltipData, tooltipLeft, tooltipTop } =
    useTooltip();

  // When swapped: genes on x, groups on y. Default: groups on x, genes on y.
  const xItems = swapAxes ? genes : groups;
  const yItems = swapAxes ? groups : genes;

  // Map groups to integer labels for compact axis
  const groupIndices = groups.map((_, i) => i + 1);
  const groupToIndex = Object.fromEntries(groups.map((g, i) => [g, i + 1]));
  const indexToGroup = Object.fromEntries(groups.map((g, i) => [i + 1, g]));

  // Determine whether integer labels are used on each axis
  const useIntX = !swapAxes && !showLabels;
  const useIntY = swapAxes && !showLabels;

  const xDomain = swapAxes ? genes : (showLabels ? groups : groupIndices);
  const yDomain = swapAxes ? (showLabels ? groups : groupIndices) : genes;

  // Estimate bottom margin from longest x label when showing rotated labels
  const xLabelItems = showLabels ? xItems : [];
  const maxXLabelLen = xLabelItems.length > 0 ? Math.max(...xLabelItems.map((s) => s.length)) : 0;
  const bottomMargin = maxXLabelLen > 0 ? Math.max(40, maxXLabelLen * 6 + 20) : 40;

  // Estimate left margin from longest y label
  const yLabelItems = swapAxes ? (showLabels ? groups : groupIndices.map(String)) : genes;
  const maxYLabelLen = yLabelItems.length > 0 ? Math.max(...yLabelItems.map((s) => String(s).length)) : 0;
  const leftMargin = Math.max(40, maxYLabelLen * 7 + 16);

  const MARGIN = { top: 16, right: 100, bottom: bottomMargin, left: leftMargin };

  const xMax = width - MARGIN.left - MARGIN.right;
  const yMax = height - MARGIN.top - MARGIN.bottom;

  const xScale = scaleBand({ domain: xDomain, range: [0, xMax], padding: 0.05 });
  const yScale = scaleBand({ domain: yDomain, range: [0, yMax], padding: 0.05 });

  const maxMean = Math.max(...data.map((d) => d.meanExpression), 0.01);

  const maxR = Math.min(MAX_RADIUS, xScale.bandwidth() / 2, yScale.bandwidth() / 2);
  const radiusScale = scaleSqrt({
    domain: [0, 1],
    range: [0, maxR],
  });

  const viridis = COLOR_SCALES.viridis;
  const colorScale = (val) => rgbToString(interpolateColorScale(val / maxMean, viridis));

  // Helpers to map data point to x/y keys
  const getXKey = (d) => {
    if (swapAxes) return d.gene;
    return showLabels ? d.group : groupToIndex[d.group];
  };
  const getYKey = (d) => {
    if (swapAxes) return showLabels ? d.group : groupToIndex[d.group];
    return d.gene;
  };

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
          {yDomain.map((val) => (
            <line
              key={`h-${val}`}
              x1={0}
              x2={xMax}
              y1={(yScale(val) ?? 0) + yScale.bandwidth() / 2}
              y2={(yScale(val) ?? 0) + yScale.bandwidth() / 2}
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
            const cx = (xScale(getXKey(d)) ?? 0) + xScale.bandwidth() / 2;
            const cy = (yScale(getYKey(d)) ?? 0) + yScale.bandwidth() / 2;
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
            numTicks={xDomain.length}
            tickComponent={({ x, y, formattedValue }) =>
              useIntX ? (
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
              ) : (
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
              )
            }
          />
          <AxisLeft
            scale={yScale}
            numTicks={yDomain.length}
            tickComponent={({ x, y, formattedValue }) =>
              useIntY ? (
                <text
                  x={x}
                  y={y}
                  fontSize={11}
                  textAnchor="end"
                  dx={-4}
                  dy={4}
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
              ) : (
                <text
                  x={x}
                  y={y}
                  fontSize={11}
                  textAnchor="end"
                  dx={-4}
                  dy={4}
                >
                  {formattedValue}
                </text>
              )
            }
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
            {[0.2, 0.4, 0.6, 0.8, 1.0].map((frac, i) => {
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
