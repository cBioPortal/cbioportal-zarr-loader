import { useState, useMemo, useRef, useCallback, useLayoutEffect } from "react";
import { Typography, Button, Select, Popover } from "antd";
import { ExpandOutlined, CompressOutlined, SelectOutlined, EditOutlined, CloseCircleOutlined, SaveOutlined, SettingOutlined, HeatMapOutlined, DotChartOutlined } from "@ant-design/icons";
import DeckGL from "@deck.gl/react";
import { ScatterplotLayer } from "@deck.gl/layers";
import { HexagonLayer } from "@deck.gl/aggregation-layers";
import { OrthographicView } from "@deck.gl/core";
import useAppStore from "../store/useAppStore";
import { calculatePlotDimensions } from "../utils/calculatePlotDimensions";
import { COLOR_SCALES, CATEGORICAL_COLORS } from "../utils/colors";
import {
  pointInPolygon,
  simplifyPolygon,
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
  onSaveSelection,
  showHexbinToggle = false,
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
    tooltipColumns,
    toggleTooltipColumn,
    tooltipColumnLoading,
    metadata,
    // Selection
    selectedPointIndices,
    setSelectedPoints,
    clearSelectedPoints,
    selectionGeometry,
    setSelectionGeometry,
    // Color scale
    colorScaleName,
    setColorScaleName,
  } = useAppStore();

  const [hoverInfo, setHoverInfo] = useState(null);
  const [expanded, setExpanded] = useState(false);
  const [layerMode, setLayerMode] = useState(showHexbinToggle ? "hexbin" : "scatter");
  const [selectMode, setSelectMode] = useState("pan");
  const [hoveredCategory, setHoveredCategory] = useState(null);
  const [hoveredExpression, setHoveredExpression] = useState(null);
  const [hoveredTooltipFilter, setHoveredTooltipFilter] = useState(null);
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
      const last = lassoPointsRef.current[lassoPointsRef.current.length - 1];
      if ((pos.x - last.x) ** 2 + (pos.y - last.y) ** 2 < 25) return; // 5px min distance
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
        setSelectionGeometry({ type: "rectangle", bounds: [minWx, minWy, maxWx, maxWy] });
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
        setSelectionGeometry({ type: "lasso", polygon: simplifyPolygon(worldPolygon) });
        setSelectedPoints(indices);
      }
      lassoPointsRef.current = [];
      const svg = lassoSvgRef.current;
      if (svg) svg.style.display = "none";
    }
  }, [selectMode, points, setSelectedPoints, setSelectionGeometry, updateSelectionRect]);

  // Build selected indices set for fast lookup in getFillColor
  const selectedSet = useMemo(
    () => new Set(selectedPointIndices),
    [selectedPointIndices],
  );

  // Determine hexbin coloring mode and config
  const hexColorMode = geneExpression ? "expression" : colorData ? "category" : "density";

  const hexColorConfig = useMemo(() => {
    if (hexColorMode === "expression") {
      return {
        getColorWeight: (d) => d.expression ?? 0,
        colorAggregation: "MEAN",
        colorRange: COLOR_SCALES[colorScaleName],
        colorDomain: expressionRange ? [expressionRange.min, expressionRange.max] : undefined,
        colorScaleType: "linear",
      };
    }
    if (hexColorMode === "category") {
      // Build list of unique categories in order
      const uniqueCats = [...new Set(points.map((p) => p.category))].sort();
      const catToIndex = Object.fromEntries(uniqueCats.map((c, i) => [c, i]));
      const catColors = uniqueCats.map((_, i) => CATEGORICAL_COLORS[i % CATEGORICAL_COLORS.length]);
      return {
        getColorValue: (pts) => {
          // Find dominant category in this bin
          const counts = {};
          for (const p of pts) {
            counts[p.category] = (counts[p.category] || 0) + 1;
          }
          let maxCount = 0, dominant = pts[0].category;
          for (const [cat, cnt] of Object.entries(counts)) {
            if (cnt > maxCount) { maxCount = cnt; dominant = cat; }
          }
          return catToIndex[dominant] ?? 0;
        },
        colorRange: catColors,
        colorDomain: [0, catColors.length - 1],
        colorScaleType: "ordinal",
        _uniqueCats: uniqueCats,
      };
    }
    // Density fallback
    return {
      colorRange: [
        [237, 248, 233],
        [199, 233, 192],
        [161, 217, 155],
        [116, 196, 118],
        [49, 163, 84],
        [0, 109, 44],
      ],
    };
  }, [hexColorMode, colorScaleName, expressionRange, points]);

  const hexHover = useCallback((info) => {
    if (!info.object) { setHoverInfo(null); return; }
    const pts = info.object.points;
    const count = info.object.count ?? pts?.length ?? info.object.colorValue;
    const hex = { hexCount: count };

    if (pts && pts.length > 0) {
      // Each entry in pts may be the original point directly or wrapped in { source }
      const unwrap = (p) => p.source ?? p;
      if (hexColorMode === "expression") {
        const sum = pts.reduce((s, p) => s + (unwrap(p).expression ?? 0), 0);
        hex.meanExpression = sum / pts.length;
      }
      if (hexColorMode === "category") {
        const counts = {};
        for (const p of pts) {
          const cat = unwrap(p).category;
          counts[cat] = (counts[cat] || 0) + 1;
        }
        let maxCount = 0;
        for (const [cat, cnt] of Object.entries(counts)) {
          if (cnt > maxCount) { maxCount = cnt; hex.dominantCategory = cat; hex.dominantCount = cnt; }
        }
      }
      // Aggregate tooltip column breakdowns
      if (tooltipData && Object.keys(tooltipData).length > 0) {
        const breakdowns = {};
        for (const [col, values] of Object.entries(tooltipData)) {
          const counts = {};
          for (const p of pts) {
            const val = String(values[unwrap(p).index] ?? "");
            counts[val] = (counts[val] || 0) + 1;
          }
          // Sort by count descending, take top 5
          breakdowns[col] = Object.entries(counts)
            .sort((a, b) => b[1] - a[1])
            .slice(0, 5);
        }
        hex.tooltipBreakdowns = breakdowns;
      }
    }
    setHoverInfo({ x: info.x, y: info.y, object: hex });
  }, [hexColorMode, tooltipData]);

  const layers = layerMode === "hexbin"
    ? [
        new HexagonLayer({
          id: "hexbin",
          data: points,
          getPosition: (d) => d.position,
          gpuAggregation: false,
          radius: 0.3,
          elevationScale: 0,
          extruded: false,
          coverage: 0.9,
          opacity: 0.8,
          pickable: true,
          onHover: hexHover,
          ...hexColorConfig,
          updateTriggers: {
            getColorWeight: [geneExpression, expressionRange],
            getColorValue: [colorData],
          },
        }),
      ]
    : [
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
          getRadius: (d) => {
            if (hoveredCategory != null && d.category === hoveredCategory) return 3;
            if (hoveredExpression != null && d.expression != null && expressionRange) {
              const tolerance = (expressionRange.max - expressionRange.min) * 0.05;
              if (Math.abs(d.expression - hoveredExpression) <= tolerance) return 3;
            }
            if (hoveredTooltipFilter != null) {
              const colValues = tooltipData[hoveredTooltipFilter.col];
              if (colValues && String(colValues[d.index]) === hoveredTooltipFilter.value) return 3;
            }
            return 1;
          },
          radiusUnits: "pixels",
          radiusMinPixels: 0.5,
          radiusMaxPixels: (hoveredCategory != null || hoveredExpression != null || hoveredTooltipFilter != null) ? 3 : 1,
          opacity: 0.7,
          pickable: true,
          onHover: (info) => setHoverInfo(info.object ? info : null),
          updateTriggers: {
            getFillColor: [colorData, geneExpression, expressionRange, colorScaleName, selectedPointIndices],
            getRadius: [hoveredCategory, hoveredExpression, hoveredTooltipFilter],
          },
        }),
      ];

  const sortedCategories = useMemo(
    () => colorData ? sortCategoriesByCount(categoryColorMap, points) : [],
    [categoryColorMap, colorData, points],
  );

  const handleTooltipChange = useCallback((newValues) => {
    const oldSet = new Set(tooltipColumns);
    const newSet = new Set(newValues);
    // Added columns
    for (const col of newValues) {
      if (!oldSet.has(col)) toggleTooltipColumn(col);
    }
    // Removed columns
    for (const col of tooltipColumns) {
      if (!newSet.has(col)) toggleTooltipColumn(col);
    }
  }, [tooltipColumns, toggleTooltipColumn]);

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
      {points.length < shape?.[0] && (
        <Text type="secondary" style={{ marginBottom: 16, display: "block" }}>
          Showing {points.length.toLocaleString()} of {shape[0].toLocaleString()} points
        </Text>
      )}

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
            {showHexbinToggle && (
              <Button
                size="small"
                type={layerMode === "hexbin" ? "primary" : "default"}
                icon={layerMode === "hexbin" ? <HeatMapOutlined /> : <DotChartOutlined />}
                onClick={() => setLayerMode(layerMode === "hexbin" ? "scatter" : "hexbin")}
                style={{ opacity: 0.85 }}
                title={layerMode === "hexbin" ? "Switch to scatter" : "Switch to hexbin"}
              />
            )}
          </div>
          <div style={{ position: "absolute", top: 8, right: 8, zIndex: 1, display: "flex", gap: 4 }} onMouseDown={(e) => e.stopPropagation()}>
            {geneExpression && (
              <Popover
                trigger="click"
                placement="bottomRight"
                content={
                  <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
                    <Text>Color scale:</Text>
                    <Select
                      size="small"
                      value={colorScaleName}
                      onChange={setColorScaleName}
                      style={{ width: 100 }}
                      options={Object.keys(COLOR_SCALES).map((name) => ({ label: name, value: name }))}
                    />
                  </div>
                }
              >
                <Button
                  size="small"
                  icon={<SettingOutlined />}
                  style={{ opacity: 0.85 }}
                  title="Plot settings"
                />
              </Popover>
            )}
            <Button
              size="small"
              icon={expanded ? <CompressOutlined /> : <ExpandOutlined />}
              onClick={() => setExpanded(!expanded)}
              style={{ opacity: 0.85 }}
            />
          </div>
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
              onMouseDown={(e) => e.stopPropagation()}
            >
              {selectedPointIndices.length.toLocaleString()} selected
              {selectionGeometry && onSaveSelection && (
                <SaveOutlined
                  onClick={onSaveSelection}
                  style={{ cursor: "pointer", fontSize: 14 }}
                  title="Save selection to config"
                />
              )}
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
          <CollapsibleLegend categories={sortedCategories} maxHeight={containerSize.height} onHoverCategory={setHoveredCategory} />
        )}

        {/* Gene expression color scale */}
        {geneExpression && expressionRange && (
          <ExpressionLegend
            selectedGene={selectedGene}
            expressionRange={expressionRange}
            colorScaleName={colorScaleName}
            onHoverExpression={setHoveredExpression}
          />
        )}

        {/* Selection summary */}
        {selectionSummary && (
          <SelectionSummaryPanel
            selectionSummary={selectionSummary}
            selectedCount={selectedPointIndices.length || points.length}
            categoryColorMap={categoryColorMap}
            colorColumn={colorColumn}
            selectedGene={selectedGene}
            maxHeight={containerSize.height}
            onHoverCategory={setHoveredCategory}
            onHoverTooltipValue={setHoveredTooltipFilter}
            obsColumns={metadata?.obsColumns}
            tooltipColumns={tooltipColumns}
            onTooltipChange={handleTooltipChange}
            tooltipColumnLoading={tooltipColumnLoading}
          />
        )}

      </div>
    </>
  );
}
