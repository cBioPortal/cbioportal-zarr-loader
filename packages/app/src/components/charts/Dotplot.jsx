import { Group } from "@visx/group";
import { scaleBand, scaleLinear } from "@visx/scale";
import { AxisBottom, AxisLeft } from "@visx/axis";
import { scaleSequential } from "d3-scale";
import { interpolateViridis } from "d3-scale-chromatic";

const MARGIN = { top: 16, right: 16, bottom: 80, left: 120 };
const MAX_RADIUS = 14;

export default function Dotplot({ genes, groups, data, width = 600, height = 400 }) {
  const xMax = width - MARGIN.left - MARGIN.right;
  const yMax = height - MARGIN.top - MARGIN.bottom;

  const xScale = scaleBand({ domain: genes, range: [0, xMax], padding: 0.2 });
  const yScale = scaleBand({ domain: groups, range: [0, yMax], padding: 0.2 });

  const maxMean = Math.max(...data.map((d) => d.meanExpression), 0.01);

  const radiusScale = scaleLinear({
    domain: [0, 1],
    range: [0, Math.min(MAX_RADIUS, xScale.bandwidth() / 2, yScale.bandwidth() / 2)],
  });

  const colorScale = scaleSequential(interpolateViridis).domain([0, maxMean]);

  return (
    <svg width={width} height={height}>
      <Group left={MARGIN.left} top={MARGIN.top}>
        {data.map((d) => {
          const cx = (xScale(d.gene) ?? 0) + xScale.bandwidth() / 2;
          const cy = (yScale(d.group) ?? 0) + yScale.bandwidth() / 2;
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
            />
          );
        })}
        <AxisBottom
          scale={xScale}
          top={yMax}
          tickLabelProps={() => ({
            fontSize: 11,
            textAnchor: "end",
            angle: -45,
            dx: -4,
            dy: -2,
          })}
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
  );
}
