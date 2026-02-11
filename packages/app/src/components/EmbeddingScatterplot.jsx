import { useState, useMemo, useRef, useCallback, useLayoutEffect } from "react";
import { Typography, Space, Button, Select } from "antd";
import { ExpandOutlined, CompressOutlined, SelectOutlined, EditOutlined, CloseCircleOutlined } from "@ant-design/icons";
import DeckGL from "@deck.gl/react";
import { ScatterplotLayer } from "@deck.gl/layers";
import { OrthographicView } from "@deck.gl/core";
import useAppStore from "../store/useAppStore";
import { calculatePlotDimensions } from "../utils/calculatePlotDimensions";
import { COLOR_SCALES } from "../utils/colors";
import {
  pointInPolygon,
  buildScatterplotPoints,
  buildSelectionSummary,
  computeRange,
  computeViewState,
  getPointFillColor,
  sortCategoriesByCount,
} from "../utils/scatterplotUtils";
import HoverTooltip from "./HoverTooltip";
import ExpressionLegend from "./ExpressionLegend";
import SelectionSummaryPanel from "./SelectionSummaryPanel";
import CollapsibleLegend from "./CollapsibleLegend";

const { Text } = Typography;

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

  const expressionRange = useMemo(
    () => computeRange(geneExpression),
    [geneExpression],
  );

  const { points, categoryColorMap, bounds } = useMemo(
    () => buildScatterplotPoints({ data, shape, maxPoints, colorData, geneExpression }),
    [data, shape, colorData, geneExpression, maxPoints],
  );

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

  const initialViewState = useMemo(
    () => computeViewState(bounds, containerSize),
    [bounds, containerSize],
  );

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
      getFillColor: (d) => getPointFillColor(d, {
        selectedSet,
        geneExpression,
        expressionRange,
        hasColorData: !!colorData,
        colorScale: COLOR_SCALES[colorScaleName],
      }),
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

  const sortedCategories = useMemo(
    () => colorData ? sortCategoriesByCount(categoryColorMap, points) : [],
    [categoryColorMap, colorData, points],
  );

  const selectionSummary = useMemo(
    () => buildSelectionSummary({
      selectedSet,
      points,
      hasColorData: !!colorData,
      hasGeneExpression: !!geneExpression,
      tooltipData,
    }),
    [selectedSet, points, colorData, geneExpression, tooltipData],
  );

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
            overflow: "hidden",
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
          <div style={{ position: "absolute", top: 8, left: 8, zIndex: 1, display: "flex", gap: 4 }} onMouseDown={(e) => e.stopPropagation()}>
            <Button
              size="small"
              type={selectMode === "rectangle" ? "primary" : "default"}
              icon={<SelectOutlined />}
              onClick={() => {
                if (selectMode === "rectangle") {
                  setSelectMode("pan");
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
            onMouseDown={(e) => e.stopPropagation()}
            style={{
              position: "absolute",
              top: 8,
              right: 8,
              zIndex: 1,
              opacity: 0.85,
            }}
          />
          {hoverInfo && (
            <HoverTooltip
              hoverInfo={hoverInfo}
              colorColumn={colorColumn}
              selectedGene={selectedGene}
              tooltipData={tooltipData}
              hasColorData={!!colorData}
              hasGeneExpression={!!geneExpression}
            />
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
                display: "flex",
                alignItems: "center",
                gap: 6,
              }}
            >
              {selectedPointIndices.length.toLocaleString()} selected
              <CloseCircleOutlined
                onClick={clearSelectedPoints}
                style={{ cursor: "pointer", fontSize: 14 }}
                title="Clear selection"
              />
            </div>
          )}
        </div>

        {/* Legend */}
        {colorColumn && sortedCategories.length > 1 && (
          <CollapsibleLegend categories={sortedCategories} maxHeight={containerSize.height} />
        )}

        {/* Gene expression color scale */}
        {geneExpression && expressionRange && (
          <ExpressionLegend
            selectedGene={selectedGene}
            expressionRange={expressionRange}
            colorScaleName={colorScaleName}
          />
        )}

        {/* Selection summary */}
        {selectionSummary && (
          <SelectionSummaryPanel
            selectionSummary={selectionSummary}
            selectedCount={selectedPointIndices.length}
            categoryColorMap={categoryColorMap}
            colorColumn={colorColumn}
            selectedGene={selectedGene}
            maxHeight={containerSize.height}
          />
        )}

      </div>
    </>
  );
}
