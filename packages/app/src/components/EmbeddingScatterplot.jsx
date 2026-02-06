import { useState, useMemo, useRef, useLayoutEffect } from "react";
import { Typography, Space, Select, Tag, Card, Input } from "antd";
import DeckGL from "@deck.gl/react";
import { ScatterplotLayer } from "@deck.gl/layers";
import { OrthographicView } from "@deck.gl/core";
import useAppStore from "../store/useAppStore";
import { calculatePlotDimensions } from "../utils/calculatePlotDimensions";

const { Text } = Typography;

// Viridis-like color scale for continuous data
const VIRIDIS = [
  [68, 1, 84],
  [72, 40, 120],
  [62, 74, 137],
  [49, 104, 142],
  [38, 130, 142],
  [31, 158, 137],
  [53, 183, 121],
  [109, 205, 89],
  [180, 222, 44],
  [253, 231, 37],
];

function interpolateViridis(t) {
  const idx = t * (VIRIDIS.length - 1);
  const lower = Math.floor(idx);
  const upper = Math.min(lower + 1, VIRIDIS.length - 1);
  const frac = idx - lower;
  return [
    Math.round(VIRIDIS[lower][0] + (VIRIDIS[upper][0] - VIRIDIS[lower][0]) * frac),
    Math.round(VIRIDIS[lower][1] + (VIRIDIS[upper][1] - VIRIDIS[lower][1]) * frac),
    Math.round(VIRIDIS[lower][2] + (VIRIDIS[upper][2] - VIRIDIS[lower][2]) * frac),
  ];
}

const COLORS = [
  [31, 119, 180],   // #1f77b4
  [255, 127, 14],   // #ff7f0e
  [44, 160, 44],    // #2ca02c
  [214, 39, 40],    // #d62728
  [148, 103, 189],  // #9467bd
  [140, 86, 75],    // #8c564b
  [227, 119, 194],  // #e377c2
  [127, 127, 127],  // #7f7f7f
  [188, 189, 34],   // #bcbd22
  [23, 190, 207],   // #17becf
  [174, 199, 232],  // #aec7e8
  [255, 187, 120],  // #ffbb78
  [152, 223, 138],  // #98df8a
  [255, 152, 150],  // #ff9896
  [197, 176, 213],  // #c5b0d5
];

export default function EmbeddingScatterplot({
  data,
  shape,
  label,
  maxPoints = 50000,
}) {
  const {
    metadata,
    // Color by obs column
    colorColumn,
    colorData,
    colorLoading,
    setColorColumn,
    // Gene expression
    selectedGene,
    geneExpression,
    geneLoading,
    setSelectedGene,
    clearGeneSelection,
  } = useAppStore();

  const { obsColumns, geneNames } = metadata;

  const [hoverInfo, setHoverInfo] = useState(null);
  const [geneSearchText, setGeneSearchText] = useState("");
  const containerRef = useRef(null);
  const [containerSize, setContainerSize] = useState({ width: 600, height: 600 });

  // Compute expression range for normalization
  const expressionRange = useMemo(() => {
    if (!geneExpression) return null;
    let min = Infinity, max = -Infinity;
    for (let i = 0; i < geneExpression.length; i++) {
      const val = geneExpression[i];
      if (val < min) min = val;
      if (val > max) max = val;
    }
    return { min, max };
  }, [geneExpression]);

  const { points, categoryColorMap, bounds } = useMemo(() => {
    if (!data || !shape) return { points: [], categoryColorMap: {}, bounds: null };

    const cols = shape[1];
    const step = Math.max(1, Math.floor(shape[0] / maxPoints));

    // Build category color map (for obs columns)
    const categories = new Map();
    if (colorData) {
      for (let i = 0; i < shape[0]; i += step) {
        const cat = String(colorData[i]);
        if (!categories.has(cat)) {
          categories.set(cat, (categories.size % COLORS.length));
        }
      }
    }

    // Build points array and calculate bounds
    const pts = [];
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;

    for (let i = 0; i < shape[0]; i += step) {
      const x = data[i * cols];
      const y = -data[i * cols + 1]; // Flip Y axis

      minX = Math.min(minX, x);
      maxX = Math.max(maxX, x);
      minY = Math.min(minY, y);
      maxY = Math.max(maxY, y);

      const cat = colorData ? String(colorData[i]) : "All";
      const exprValue = geneExpression ? geneExpression[i] : null;

      pts.push({
        position: [x, y],
        category: cat,
        colorIndex: categories.get(cat) ?? 0,
        index: i,
        expression: exprValue,
      });
    }

    // Convert categories map to object for legend
    const colorMap = {};
    categories.forEach((colorIdx, cat) => {
      colorMap[cat] = COLORS[colorIdx];
    });

    return {
      points: pts,
      categoryColorMap: colorMap,
      bounds: { minX, maxX, minY, maxY },
    };
  }, [data, shape, colorData, geneExpression, maxPoints]);

  // Calculate container dimensions based on available space and data aspect ratio
  useLayoutEffect(() => {
    const updateSize = () => {
      if (!containerRef.current || !bounds) return;

      const parentWidth = containerRef.current.parentElement?.clientWidth || 800;
      // Reserve space for gene list (220px) and legend (~150px) and gaps
      const availableWidth = Math.min(parentWidth - 400, 800);

      const dimensions = calculatePlotDimensions({
        bounds,
        availableWidth,
        maxHeight: 600,
        minWidth: 400,
      });

      setContainerSize(dimensions);
    };

    updateSize();
    window.addEventListener("resize", updateSize);
    return () => window.removeEventListener("resize", updateSize);
  }, [bounds]);

  // Filtered genes for the gene list
  const filteredGenes = useMemo(() => {
    if (!geneNames) return [];
    if (!geneSearchText) return geneNames;
    const search = geneSearchText.toLowerCase();
    return geneNames.filter(name => name.toLowerCase().includes(search));
  }, [geneNames, geneSearchText]);

  const initialViewState = useMemo(() => {
    if (!bounds) return { target: [0, 0], zoom: 0 };

    const centerX = (bounds.minX + bounds.maxX) / 2;
    const centerY = (bounds.minY + bounds.maxY) / 2;
    const rangeX = bounds.maxX - bounds.minX;
    const rangeY = bounds.maxY - bounds.minY;
    const maxRange = Math.max(rangeX, rangeY);

    // Use the smaller container dimension to calculate zoom
    const viewSize = Math.min(containerSize.width, containerSize.height);
    const zoom = Math.log2(viewSize / maxRange) - 1;

    return {
      target: [centerX, centerY],
      zoom: Math.max(-5, Math.min(zoom, 10)),
    };
  }, [bounds, containerSize]);

  const layers = [
    new ScatterplotLayer({
      id: "scatterplot",
      data: points,
      getPosition: (d) => d.position,
      getFillColor: (d) => {
        if (geneExpression && expressionRange && expressionRange.max > expressionRange.min) {
          const t = (d.expression - expressionRange.min) / (expressionRange.max - expressionRange.min);
          return interpolateViridis(t);
        }
        if (colorData) {
          return COLORS[d.colorIndex];
        }
        return [24, 144, 255];
      },
      getRadius: 1,
      radiusUnits: "pixels",
      radiusMinPixels: 0.5,
      radiusMaxPixels: 1,
      opacity: 0.7,
      pickable: true,
      onHover: (info) => setHoverInfo(info.object ? info : null),
      updateTriggers: {
        getFillColor: [colorData, geneExpression, expressionRange],
      },
    }),
  ];

  const sortedCategories = useMemo(() => {
    if (!colorData || Object.keys(categoryColorMap).length === 0) return [];
    return Object.entries(categoryColorMap)
      .sort((a, b) => {
        const countA = points.filter(p => p.category === a[0]).length;
        const countB = points.filter(p => p.category === b[0]).length;
        return countB - countA;
      });
  }, [categoryColorMap, colorData, points]);

  return (
    <>
      <Space style={{ marginBottom: 16 }} wrap align="center">
        <Text>Color by obs:</Text>
        <Select
          allowClear
          placeholder="Select obs column"
          style={{ width: 180 }}
          value={colorColumn}
          onChange={setColorColumn}
          loading={colorLoading}
          options={obsColumns?.map((c) => ({ label: c, value: c })) || []}
        />
        {selectedGene && (
          <>
            <Text type="secondary">|</Text>
            <Text>Gene:</Text>
            <Tag
              closable
              onClose={clearGeneSelection}
              color="blue"
              style={{ fontSize: 13, padding: "2px 8px" }}
            >
              {selectedGene}
            </Tag>
            {geneLoading && <Text type="secondary">Loading...</Text>}
          </>
        )}
        {points.length < shape?.[0] && (
          <>
            <Text type="secondary">|</Text>
            <Text type="secondary">
              Showing {points.length.toLocaleString()} of {shape[0].toLocaleString()} points
            </Text>
          </>
        )}
      </Space>

      <div ref={containerRef} style={{ display: "flex", gap: 16 }}>
        <div
          style={{
            width: containerSize.width,
            height: containerSize.height,
            position: "relative",
            border: "1px solid #d9d9d9",
            borderRadius: 4,
          }}
        >
          <DeckGL
            views={new OrthographicView({ id: "ortho" })}
            initialViewState={initialViewState}
            controller={true}
            layers={layers}
          />
          {hoverInfo && (
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
              <div>x: {hoverInfo.object.position[0].toFixed(4)}</div>
              <div>y: {hoverInfo.object.position[1].toFixed(4)}</div>
              {colorData && <div>{colorColumn}: {hoverInfo.object.category}</div>}
              {geneExpression && <div>{selectedGene}: {hoverInfo.object.expression?.toFixed(4)}</div>}
            </div>
          )}

          {/* Axis labels */}
          <div
            style={{
              position: "absolute",
              bottom: 4,
              left: "50%",
              transform: "translateX(-50%)",
              fontSize: 12,
              color: "#666",
            }}
          >
            {label}_1
          </div>
          <div
            style={{
              position: "absolute",
              left: 4,
              top: "50%",
              transform: "translateY(-50%) rotate(-90deg)",
              fontSize: 12,
              color: "#666",
            }}
          >
            {label}_2
          </div>
        </div>

        {/* Legend */}
        {colorColumn && sortedCategories.length > 1 && (
          <div style={{ maxHeight: containerSize.height, overflow: "auto", fontSize: 12 }}>
            {sortedCategories.map(([cat, color]) => (
              <div key={cat} style={{ display: "flex", alignItems: "center", gap: 4, marginBottom: 2 }}>
                <div
                  style={{
                    width: 12,
                    height: 12,
                    borderRadius: 2,
                    backgroundColor: `rgb(${color.join(",")})`,
                  }}
                />
                <span>{cat}</span>
              </div>
            ))}
          </div>
        )}

        {/* Gene expression color scale */}
        {geneExpression && expressionRange && (
          <div style={{ fontSize: 12 }}>
            <div style={{ marginBottom: 4 }}>{selectedGene}</div>
            <div
              style={{
                width: 20,
                height: 200,
                background: `linear-gradient(to bottom, ${VIRIDIS.slice().reverse().map(c => `rgb(${c.join(",")})`).join(", ")})`,
                borderRadius: 2,
              }}
            />
            <div style={{ display: "flex", flexDirection: "column", justifyContent: "space-between", height: 200, marginLeft: 4, position: "relative", top: -200 }}>
              <span>{expressionRange.max.toFixed(2)}</span>
              <span>{((expressionRange.max + expressionRange.min) / 2).toFixed(2)}</span>
              <span>{expressionRange.min.toFixed(2)}</span>
            </div>
          </div>
        )}

        {/* Gene list */}
        {geneNames && geneNames.length > 0 && (
          <Card
            size="small"
            title={`Genes (${geneNames.length.toLocaleString()})`}
            style={{ width: 220, height: containerSize.height }}
            bodyStyle={{ padding: 0, height: "calc(100% - 38px)", display: "flex", flexDirection: "column" }}
          >
            <div style={{ padding: 8, borderBottom: "1px solid #f0f0f0" }}>
              <Input.Search
                placeholder="Search genes..."
                size="small"
                value={geneSearchText}
                onChange={(e) => setGeneSearchText(e.target.value)}
                allowClear
              />
            </div>
            <div style={{ flex: 1, overflow: "auto" }}>
              {filteredGenes.map((gene) => (
                <div
                  key={gene}
                  style={{
                    padding: "4px 12px",
                    cursor: "pointer",
                    backgroundColor: selectedGene === gene ? "#e6f4ff" : undefined,
                    fontSize: 12,
                    whiteSpace: "nowrap",
                    overflow: "hidden",
                    textOverflow: "ellipsis",
                  }}
                  onClick={() => setSelectedGene(gene)}
                >
                  {gene}
                </div>
              ))}
              {filteredGenes.length === 0 && geneSearchText && (
                <div style={{ padding: "8px 12px", color: "#999", fontSize: 12 }}>
                  No genes found
                </div>
              )}
            </div>
            {geneSearchText && (
              <div style={{ padding: "4px 12px", fontSize: 11, color: "#999", borderTop: "1px solid #f0f0f0" }}>
                {filteredGenes.length.toLocaleString()} matches
              </div>
            )}
          </Card>
        )}
      </div>
    </>
  );
}
