import { useState, useMemo, useRef, useCallback, useLayoutEffect } from "react";
import { Typography, Space, Button, Select } from "antd";
import { ExpandOutlined, CompressOutlined, SelectOutlined, EditOutlined } from "@ant-design/icons";
import DeckGL from "@deck.gl/react";
import { ScatterplotLayer } from "@deck.gl/layers";
import { OrthographicView } from "@deck.gl/core";
import useAppStore from "../store/useAppStore";
import { calculatePlotDimensions } from "../utils/calculatePlotDimensions";
import { CATEGORICAL_COLORS, COLOR_SCALES, interpolateColorScale, colorScaleGradient } from "../utils/colors";

const { Text } = Typography;

function pointInPolygon(x, y, polygon) {
  let inside = false;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i][0], yi = polygon[i][1];
    const xj = polygon[j][0], yj = polygon[j][1];
    if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
      inside = !inside;
    }
  }
  return inside;
}

export default function EmbeddingScatterplot({
  data,
  shape,
  label,
  maxPoints = Infinity,
}) {
  const {
    // Color by obs column
    colorColumn,
    colorData,
    // Gene expression
    selectedGene,
    geneExpression,
    // Tooltip
    tooltipData,
    // Selection
    selectedPointIndices,
    setSelectedPoints,
    clearSelectedPoints,
  } = useAppStore();

  const [hoverInfo, setHoverInfo] = useState(null);
  const [expanded, setExpanded] = useState(false);
  const [colorScaleName, setColorScaleName] = useState("viridis");
  const [selectMode, setSelectMode] = useState("pan");
  const dragStartRef = useRef(null);
  const dragEndRef = useRef(null);
  const selectionRectRef = useRef(null);
  const lassoPointsRef = useRef([]);
  const lassoSvgRef = useRef(null);
  const containerRef = useRef(null);
  const deckRef = useRef(null);
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

  // Selection rectangle mouse handlers — use refs + direct DOM updates to avoid re-renders during drag
  const updateSelectionRect = useCallback(() => {
    const el = selectionRectRef.current;
    const start = dragStartRef.current;
    const end = dragEndRef.current;
    if (!el || !start || !end) {
      if (el) el.style.display = "none";
      return;
    }
    el.style.display = "block";
    el.style.left = `${Math.min(start.x, end.x)}px`;
    el.style.top = `${Math.min(start.y, end.y)}px`;
    el.style.width = `${Math.abs(end.x - start.x)}px`;
    el.style.height = `${Math.abs(end.y - start.y)}px`;
  }, []);

  const handleMouseDown = useCallback((e) => {
    if (selectMode === "pan") return;
    const rect = e.currentTarget.getBoundingClientRect();
    const pos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    if (selectMode === "rectangle") {
      dragStartRef.current = pos;
      dragEndRef.current = pos;
      updateSelectionRect();
    } else if (selectMode === "lasso") {
      lassoPointsRef.current = [pos];
      const svg = lassoSvgRef.current;
      if (svg) {
        svg.style.display = "block";
        svg.querySelector("polyline").setAttribute("points", `${pos.x},${pos.y}`);
      }
    }
  }, [selectMode, updateSelectionRect]);

  const handleMouseMove = useCallback((e) => {
    if (selectMode === "pan") return;
    const rect = e.currentTarget.getBoundingClientRect();
    const pos = { x: e.clientX - rect.left, y: e.clientY - rect.top };
    if (selectMode === "rectangle") {
      if (!dragStartRef.current) return;
      dragEndRef.current = pos;
      updateSelectionRect();
    } else if (selectMode === "lasso") {
      if (lassoPointsRef.current.length === 0) return;
      lassoPointsRef.current.push(pos);
      const svg = lassoSvgRef.current;
      if (svg) {
        const pointsStr = lassoPointsRef.current.map(p => `${p.x},${p.y}`).join(" ");
        svg.querySelector("polyline").setAttribute("points", pointsStr);
      }
    }
  }, [selectMode, updateSelectionRect]);

  const handleMouseUp = useCallback(() => {
    if (selectMode === "pan") return;

    if (selectMode === "rectangle") {
      if (!dragStartRef.current || !dragEndRef.current) return;
      const start = dragStartRef.current;
      const end = dragEndRef.current;
      const viewport = deckRef.current?.deck?.getViewports()?.[0];
      if (viewport) {
        const [wx1, wy1] = viewport.unproject([start.x, start.y]);
        const [wx2, wy2] = viewport.unproject([end.x, end.y]);
        const minWx = Math.min(wx1, wx2);
        const maxWx = Math.max(wx1, wx2);
        const minWy = Math.min(wy1, wy2);
        const maxWy = Math.max(wy1, wy2);

        const indices = [];
        for (const pt of points) {
          const [px, py] = pt.position;
          if (px >= minWx && px <= maxWx && py >= minWy && py <= maxWy) {
            indices.push(pt.index);
          }
        }
        setSelectedPoints(indices);
      }
      dragStartRef.current = null;
      dragEndRef.current = null;
      updateSelectionRect();
    } else if (selectMode === "lasso") {
      const lassoPoints = lassoPointsRef.current;
      if (lassoPoints.length < 3) {
        lassoPointsRef.current = [];
        const svg = lassoSvgRef.current;
        if (svg) svg.style.display = "none";
        return;
      }
      const viewport = deckRef.current?.deck?.getViewports()?.[0];
      if (viewport) {
        const worldPolygon = lassoPoints.map(p => viewport.unproject([p.x, p.y]));
        const indices = [];
        for (const pt of points) {
          const [px, py] = pt.position;
          if (pointInPolygon(px, py, worldPolygon)) {
            indices.push(pt.index);
          }
        }
        setSelectedPoints(indices);
      }
      lassoPointsRef.current = [];
      const svg = lassoSvgRef.current;
      if (svg) svg.style.display = "none";
    }
  }, [selectMode, points, setSelectedPoints, updateSelectionRect]);

  // Build selected indices set for fast lookup in getFillColor
  const selectedSet = useMemo(
    () => new Set(selectedPointIndices),
    [selectedPointIndices],
  );

  const layers = [
    new ScatterplotLayer({
      id: "scatterplot",
      data: points,
      getPosition: (d) => d.position,
      getFillColor: (d) => {
        const dimmed = selectedSet.size > 0 && !selectedSet.has(d.index);
        if (dimmed) return [180, 180, 180, 60];

        const scale = COLOR_SCALES[colorScaleName];
        if (geneExpression && expressionRange && expressionRange.max > expressionRange.min) {
          const t = (d.expression - expressionRange.min) / (expressionRange.max - expressionRange.min);
          return interpolateColorScale(t, scale);
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
        getFillColor: [colorData, geneExpression, expressionRange, colorScaleName, selectedPointIndices],
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
        {geneExpression && (
          <>
            <Text>Color scale:</Text>
            <Select
              size="small"
              value={colorScaleName}
              onChange={setColorScaleName}
              style={{ width: 100 }}
              options={Object.keys(COLOR_SCALES).map((name) => ({ label: name, value: name }))}
            />
          </>
        )}
        {points.length < shape?.[0] && (
          <Text type="secondary">
            Showing {points.length.toLocaleString()} of {shape[0].toLocaleString()} points
          </Text>
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
            cursor: selectMode !== "pan" ? "crosshair" : undefined,
          }}
          onMouseDown={handleMouseDown}
          onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp}
        >
          <DeckGL
            ref={deckRef}
            key={`${containerSize.width}-${containerSize.height}`}
            width={containerSize.width}
            height={containerSize.height}
            views={new OrthographicView({ id: "ortho" })}
            initialViewState={initialViewState}
            controller={{ dragPan: selectMode === "pan" }}
            layers={layers}
          />
          {/* Selection rectangle overlay — styled directly via ref to avoid re-renders */}
          <div
            ref={selectionRectRef}
            style={{
              display: "none",
              position: "absolute",
              backgroundColor: "rgba(24, 144, 255, 0.15)",
              border: "1px solid rgba(24, 144, 255, 0.6)",
              pointerEvents: "none",
              zIndex: 2,
            }}
          />
          {/* Lasso SVG overlay — updated directly via ref to avoid re-renders */}
          <svg
            ref={lassoSvgRef}
            style={{
              position: "absolute",
              top: 0,
              left: 0,
              width: "100%",
              height: "100%",
              pointerEvents: "none",
              zIndex: 2,
              display: "none",
            }}
          >
            <polyline
              fill="rgba(24, 144, 255, 0.15)"
              stroke="rgba(24, 144, 255, 0.6)"
              strokeWidth="1"
              points=""
            />
          </svg>
          <div style={{ position: "absolute", top: 8, left: 8, zIndex: 1, display: "flex", gap: 4 }}>
            <Button
              size="small"
              type={selectMode === "rectangle" ? "primary" : "default"}
              icon={<SelectOutlined />}
              onClick={() => {
                if (selectMode === "rectangle") {
                  setSelectMode("pan");
                  clearSelectedPoints();
                } else {
                  setSelectMode("rectangle");
                  clearSelectedPoints();
                }
              }}
              style={{ opacity: 0.85 }}
              title="Rectangle select"
            />
            <Button
              size="small"
              type={selectMode === "lasso" ? "primary" : "default"}
              icon={<EditOutlined />}
              onClick={() => {
                if (selectMode === "lasso") {
                  setSelectMode("pan");
                  clearSelectedPoints();
                } else {
                  setSelectMode("lasso");
                  clearSelectedPoints();
                }
              }}
              style={{ opacity: 0.85 }}
              title="Lasso select"
            />
          </div>
          <Button
            size="small"
            icon={expanded ? <CompressOutlined /> : <ExpandOutlined />}
            onClick={() => setExpanded(!expanded)}
            style={{
              position: "absolute",
              top: 8,
              right: 8,
              zIndex: 1,
              opacity: 0.85,
            }}
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
              {Object.entries(tooltipData).map(([col, values]) => (
                <div key={col}>{col}: {values[hoverInfo.object.index]}</div>
              ))}
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
          {selectedPointIndices.length > 0 && (
            <div
              style={{
                position: "absolute",
                bottom: 8,
                right: 8,
                background: "rgba(0,0,0,0.7)",
                color: "white",
                padding: "2px 8px",
                borderRadius: 4,
                fontSize: 12,
                zIndex: 1,
              }}
            >
              {selectedPointIndices.length.toLocaleString()} selected
            </div>
          )}
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
        )}

      </div>
    </>
  );
}
