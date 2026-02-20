export default function HoverTooltip({
  hoverInfo,
  colorColumn,
  selectedGene,
  tooltipData,
  hasColorData,
  hasGeneExpression,
}) {
  return (
    <div
      style={{
        position: "absolute",
        left: hoverInfo.x + 10,
        top: hoverInfo.y + 10,
        background: "rgba(0,0,0,0.8)",
        color: "white",
        padding: "4px 8px",
        borderRadius: 4,
        fontSize: 12,
        pointerEvents: "none",
      }}
    >
      {hoverInfo.object.position ? (
        <>
          <div>x: {hoverInfo.object.position[0].toFixed(4)}</div>
          <div>y: {hoverInfo.object.position[1].toFixed(4)}</div>
          {hasColorData && <div>{colorColumn}: {hoverInfo.object.category}</div>}
          {hasGeneExpression && <div>{selectedGene}: {hoverInfo.object.expression?.toFixed(4)}</div>}
          {Object.entries(tooltipData).map(([col, values]) => (
            <div key={col}>{col}: {values[hoverInfo.object.index]}</div>
          ))}
        </>
      ) : hoverInfo.object.hexCount != null ? (
        <>
          <div>Count: {hoverInfo.object.hexCount.toLocaleString()}</div>
          {hoverInfo.object.meanExpression != null && (
            <div>Mean {selectedGene}: {hoverInfo.object.meanExpression.toFixed(4)}</div>
          )}
          {hoverInfo.object.dominantCategory != null && (
            <div>{colorColumn}: {hoverInfo.object.dominantCategory} ({hoverInfo.object.dominantCount})</div>
          )}
        </>
      ) : (
        <div>Bin selected</div>
      )}
    </div>
  );
}
