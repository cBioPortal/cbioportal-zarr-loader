import { useState, useMemo, useRef, useEffect, useCallback } from "react";
import { Button } from "antd";
import {
  SelectOutlined,
  EditOutlined,
  CloseCircleOutlined,
  ExpandOutlined,
  CompressOutlined,
} from "@ant-design/icons";
import DeckGL from "@deck.gl/react";
import { ScatterplotLayer } from "@deck.gl/layers";
import { OrthographicView } from "@deck.gl/core";
import { _StatsWidget as StatsWidget } from "@deck.gl/widgets";
import "@deck.gl/widgets/stylesheet.css";
import { Stats } from "@probe.gl/stats";
import { computeBoundsFromBuffer, computeViewState } from "../../utils/scatterplotUtils";
import SelectionOverlay from "../ui/SelectionOverlay";
import useSelectionInteraction from "../../hooks/useSelectionInteraction";

export type SelectionGeometry =
  | { type: "rectangle"; bounds: [number, number, number, number] }
  | { type: "lasso"; polygon: [number, number][] };

interface EmbeddingScatterplotGLProps {
  data: Float32Array;
  shape: [number, number];
  label: string;
  selectedPointIndices: number[];
  setSelectedPoints: (indices: number[]) => void;
  clearSelectedPoints: () => void;
  selectionGeometry: SelectionGeometry | null;
  setSelectionGeometry: (geo: SelectionGeometry | null) => void;
  debugMode?: boolean;
}

export default function EmbeddingScatterplotGL({
  data,
  shape,
  label,
  selectedPointIndices,
  setSelectedPoints,
  clearSelectedPoints,
  selectionGeometry,
  setSelectionGeometry,
  debugMode = false,
}: EmbeddingScatterplotGLProps) {
  const [expanded, setExpanded] = useState(true);
  const containerRef = useRef<HTMLDivElement>(null);
  const deckRef = useRef(null);
  const [containerSize, setContainerSize] = useState({ width: 600, height: 400 });

  // Debug stats widget — mirrors deck.gl _onMetrics into a stable Stats object
  const customStatsRef = useRef<InstanceType<typeof Stats> | null>(null);
  if (debugMode && !customStatsRef.current) {
    customStatsRef.current = new Stats({ id: "deck.gl-metrics" });
  }

  const debugWidgets = useMemo(
    () =>
      debugMode && customStatsRef.current
        ? [
            new StatsWidget({
              id: "deck-stats",
              type: "custom",
              stats: customStatsRef.current,
              title: "Deck Stats",
              framesPerUpdate: 1,
              formatters: {
                fps: "fps",
                "CPU Time": "averageTime",
                "GPU Time": "averageTime",
                "Redraw Count": "count",
                "Pick Count": "count",
                "setProps Time": "totalTime",
                "Pick Time": "totalTime",
              },
            }),
          ]
        : [],
    [debugMode],
  );

  const onDebugMetrics = useCallback((metrics: Record<string, number>) => {
    console.table(metrics);
    const s = customStatsRef.current;
    if (!s) return;
    const fpsStat = s.get("fps", "count");
    fpsStat.reset();
    fpsStat.addCount(Math.round(metrics.fps));
    const cpuStat = s.get("CPU Time", "totalTime");
    cpuStat.reset();
    cpuStat.addTime(metrics.cpuTime ?? 0);
    const gpuStat = s.get("GPU Time", "totalTime");
    gpuStat.reset();
    gpuStat.addTime(metrics.gpuTime ?? 0);
    const redrawStat = s.get("Redraw Count", "count");
    redrawStat.reset();
    redrawStat.addCount(metrics.framesRedrawn ?? 0);
    const pickCountStat = s.get("Pick Count", "count");
    pickCountStat.reset();
    pickCountStat.addCount(metrics.pickCount ?? 0);
    const setPropsStat = s.get("setProps Time", "totalTime");
    setPropsStat.reset();
    setPropsStat.addTime(metrics.setPropsTime ?? 0);
    const pickTimeStat = s.get("Pick Time", "totalTime");
    pickTimeStat.reset();
    pickTimeStat.addTime(metrics.pickTime ?? 0);
    for (const w of debugWidgets) {
      w.updateHTML();
    }
  }, [debugWidgets]);

  // Fill available space via ResizeObserver
  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const observer = new ResizeObserver((entries) => {
      const { width, height } = entries[0].contentRect;
      if (width > 0 && height > 0) {
        setContainerSize({ width: Math.round(width), height: Math.round(height) });
      }
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, []);

  const bounds = useMemo(
    () => computeBoundsFromBuffer(data, shape),
    [data, shape],
  );

  const initialViewState = useMemo(
    () => ({
      ...computeViewState(bounds, containerSize),
      minZoom: -5,
      maxZoom: 20,
    }),
    [bounds, containerSize],
  );

  const numPoints = shape[0];

  const {
    selectMode,
    setSelectMode,
    selectionRectRef,
    lassoSvgRef,
    handleMouseDown,
    handleMouseMove,
    handleMouseUp,
  } = useSelectionInteraction({
    deckRef,
    positionBuffer: data,
    stride: shape[1],
    numPoints,
    setSelectedPoints,
    setSelectionGeometry,
    clearSelectedPoints,
  });

  // Stable data descriptor — same reference unless data actually changes.
  // Prevents deck.gl from re-uploading GPU buffers on unrelated re-renders.
  const layerData = useMemo(
    () => ({
      length: numPoints,
      attributes: {
        getPosition: { value: data, size: 2 },
      },
    }),
    [data, numPoints],
  );

  const layers = useMemo(
    () => [
      new ScatterplotLayer({
        id: "scatterplot",
        data: layerData,
        dataComparator: (newData, oldData) => newData === oldData,
        getFillColor: [180, 180, 180, 200],
        getRadius: 1,
        radiusUnits: "pixels",
        radiusMinPixels: 0.5,
        radiusMaxPixels: 2,
        opacity: 0.7,
        pickable: false,
      }),
    ],
    [layerData],
  );

  return (
    <div
      ref={containerRef}
      style={{
        width: "100%",
        height: expanded ? "90vh" : "100%",
        minHeight: 300,
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
        width={containerSize.width}
        height={containerSize.height}
        views={new OrthographicView({ id: "ortho" })}
        initialViewState={initialViewState}
        controller={{ dragPan: selectMode === "pan" }}
        layers={layers}
        widgets={debugWidgets}
        {...(debugMode ? { _onMetrics: onDebugMetrics } : {})}
      />
      <SelectionOverlay selectionRectRef={selectionRectRef} lassoSvgRef={lassoSvgRef} />

      {/* Toolbar: select modes */}
      <div
        style={{ position: "absolute", top: 8, left: 8, zIndex: 1, display: "flex", gap: 4 }}
        onMouseDown={(e) => e.stopPropagation()}
      >
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

      {/* Expand/collapse */}
      <div
        style={{ position: "absolute", top: 8, right: 8, zIndex: 1 }}
        onMouseDown={(e) => e.stopPropagation()}
      >
        <Button
          size="small"
          icon={expanded ? <CompressOutlined /> : <ExpandOutlined />}
          onClick={() => setExpanded(!expanded)}
          style={{ opacity: 0.85 }}
        />
      </div>

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

      {/* Selection count badge */}
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
          <CloseCircleOutlined
            onClick={clearSelectedPoints}
            style={{ cursor: "pointer", fontSize: 14 }}
            title="Clear selection"
          />
        </div>
      )}
    </div>
  );
}
