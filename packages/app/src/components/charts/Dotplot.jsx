import { useCallback } from "react";
import { Group } from "@visx/group";
import { scaleBand, scaleLinear } from "@visx/scale";
import { AxisBottom, AxisLeft } from "@visx/axis";
import { useTooltip, TooltipWithBounds, defaultStyles } from "@visx/tooltip";
import { scaleSequential } from "d3-scale";
import { interpolateViridis } from "d3-scale-chromatic";

const MARGIN = { top: 16, right: 16, bottom: 40, left: 120 };
const MAX_RADIUS = 14;

const tooltipStyles = {
  ...defaultStyles,
  fontSize: 12,
  padding: "6px 10px",
};

export default function Dotplot({ genes, groups, data, width = 600, height = 400 }) {
  const { showTooltip, hideTooltip, tooltipOpen, tooltipData, tooltipLeft, tooltipTop } =
    useTooltip();

  const xMax = width - MARGIN.left - MARGIN.right;
  const yMax = height - MARGIN.top - MARGIN.bottom;

  // Map groups to integer labels for compact x-axis
  const groupIndices = groups.map((_, i) => i + 1);
  const groupToIndex = Object.fromEntries(groups.map((g, i) => [g, i + 1]));
  const indexToGroup = Object.fromEntries(groups.map((g, i) => [i + 1, g]));

  const xScale = scaleBand({ domain: groupIndices, range: [0, xMax], padding: 0.05 });
  const yScale = scaleBand({ domain: genes, range: [0, yMax], padding: 0.05 });

  const maxMean = Math.max(...data.map((d) => d.meanExpression), 0.01);

  const radiusScale = scaleLinear({
    domain: [0, 1],
    range: [0, Math.min(MAX_RADIUS, xScale.bandwidth() / 2, yScale.bandwidth() / 2)],
  });

  const colorScale = scaleSequential(interpolateViridis).domain([0, maxMean]);

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
          {data.map((d) => {
            const cx = (xScale(groupToIndex[d.group]) ?? 0) + xScale.bandwidth() / 2;
            const cy = (yScale(d.gene) ?? 0) + yScale.bandwidth() / 2;
            const r = radiusScale(d.fractionExpressing);
            if (r === 0) return null;
            return (
              <circle
                key={`${d.gene}-${d.group}`}
                cx={cx}
                cy={cy}
                r={r}
                fill={colorScale(d.meanExpression)}
                stroke="#ccc"
                strokeWidth={0.5}
                style={{ cursor: "pointer" }}
                onMouseEnter={(e) => handleMouseEnter(e, d)}
                onMouseLeave={hideTooltip}
              />
            );
          })}
          <AxisBottom
            scale={xScale}
            top={yMax}
            tickComponent={({ x, y, formattedValue }) => (
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
            )}
          />
          <AxisLeft
            scale={yScale}
            tickLabelProps={() => ({
              fontSize: 11,
              textAnchor: "end",
              dx: -4,
              dy: 4,
            })}
          />
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
              <div>Fraction: {(tooltipData.fractionExpressing * 100).toFixed(1)}%</div>
              <div>Mean expr: {tooltipData.meanExpression.toFixed(3)}</div>
            </>
          )}
        </TooltipWithBounds>
      )}
    </div>
  );
}
