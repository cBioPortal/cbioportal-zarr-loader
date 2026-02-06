import { useState, useMemo, useRef, useLayoutEffect } from "react";
import { Typography, Space, Select, Tag, Button } from "antd";
import { ExpandOutlined, CompressOutlined } from "@ant-design/icons";
import DeckGL from "@deck.gl/react";
import { ScatterplotLayer } from "@deck.gl/layers";
import { OrthographicView } from "@deck.gl/core";
import useAppStore from "../store/useAppStore";
import { calculatePlotDimensions } from "../utils/calculatePlotDimensions";
import { CATEGORICAL_COLORS, interpolateViridis, viridisGradient } from "../utils/colors";

const { Text } = Typography;

export default function EmbeddingScatterplot({
  data,
  shape,
  label,
  maxPoints = Infinity,
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
    clearGeneSelection,
  } = useAppStore();

  const { obsColumns } = metadata;

  const [hoverInfo, setHoverInfo] = useState(null);
  const [expanded, setExpanded] = useState(false);
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
          categories.set(cat, (categories.size % CATEGORICAL_COLORS.length));
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
      colorMap[cat] = CATEGORICAL_COLORS[colorIdx];
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

      // Use larger dimensions when expanded
      const availableWidth = expanded
        ? Math.min(parentWidth - 100, 1200)
        : Math.min(parentWidth - 400, 800);
      const maxHeight = expanded ? 900 : 600;
      const minWidth = expanded ? 600 : 400;

      const dimensions = calculatePlotDimensions({
        bounds,
        availableWidth,
        maxHeight,
        minWidth,
      });

      setContainerSize(dimensions);
    };

    updateSize();
    window.addEventListener("resize", updateSize);
    return () => window.removeEventListener("resize", updateSize);
  }, [bounds, expanded]);

  const initialViewState = useMemo(() => {
    if (!bounds) return { target: [0, 0], zoom: 0 };

    const centerX = (bounds.minX + bounds.maxX) / 2;
    const centerY = (bounds.minY + bounds.maxY) / 2;
    const rangeX = bounds.maxX - bounds.minX;
    const rangeY = bounds.maxY - bounds.minY;
    const maxRange = Math.max(rangeX, rangeY);

    // Use the smaller container dimension to calculate zoom
    const viewSize = Math.min(containerSize.width, containerSize.height);
    const zoom = Math.log2(viewSize / maxRange) - 0.1;

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
          return CATEGORICAL_COLORS[d.colorIndex];
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
        <Text type="secondary">|</Text>
        <Button
          size="small"
          icon={expanded ? <CompressOutlined /> : <ExpandOutlined />}
          onClick={() => setExpanded(!expanded)}
        >
          {expanded ? "Collapse" : "Expand"}
        </Button>
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
            key={`${containerSize.width}-${containerSize.height}`}
            width={containerSize.width}
            height={containerSize.height}
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
                background: viridisGradient("to bottom"),
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

      </div>
    </>
  );
}
