import { COLOR_SCALES, colorScaleGradient } from "../utils/colors";

export default function ExpressionLegend({
  selectedGene,
  expressionRange,
  colorScaleName,
}) {
  return (
    <div style={{ fontSize: 12 }}>
      <div style={{ marginBottom: 4 }}>{selectedGene}</div>
      <div
        style={{
          width: 20,
          height: 200,
          background: colorScaleGradient(COLOR_SCALES[colorScaleName], "to bottom"),
          borderRadius: 2,
        }}
      />
      <div style={{ display: "flex", flexDirection: "column", justifyContent: "space-between", height: 200, marginLeft: 4, position: "relative", top: -200 }}>
        <span>{expressionRange.max.toFixed(2)}</span>
        <span>{((expressionRange.max + expressionRange.min) / 2).toFixed(2)}</span>
        <span>{expressionRange.min.toFixed(2)}</span>
      </div>
    </div>
  );
}
