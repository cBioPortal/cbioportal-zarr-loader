import { useState, useMemo, useRef, useEffect } from "react";
import { Button } from "antd";
import {
  SelectOutlined,
  EditOutlined,
  CloseCircleOutlined,
} from "@ant-design/icons";
import DeckGL from "@deck.gl/react";
import { ScatterplotLayer } from "@deck.gl/layers";
import { OrthographicView } from "@deck.gl/core";
import { _StatsWidget as StatsWidget } from "@deck.gl/widgets";
import "@deck.gl/widgets/stylesheet.css";
// bounds/viewState utils disabled for perf baseline
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
  loading?: boolean;
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
  loading = false,
  clearSelectedPoints,
  selectionGeometry,
  setSelectionGeometry,
  debugMode = false,
}: EmbeddingScatterplotGLProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const deckRef = useRef(null);
  const [containerSize, setContainerSize] = useState({ width: 600, height: 400 });

  const debugWidgets = useMemo(
    () =>
      debugMode
        ? [new StatsWidget({ id: "deck-stats", type: "deck", framesPerUpdate: 1, placement: 'bottom-left' })]
        : [],
    [debugMode],
  );

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

  const initialViewState = useMemo(() => {
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (let i = 0; i < data.length; i += 2) {
      const x = data[i], y = data[i + 1];
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
    const rangeX = maxX - minX || 1;
    const rangeY = maxY - minY || 1;
    const zoom = Math.log2(Math.min(containerSize.width / rangeX, containerSize.height / rangeY)) + 0.1;
    return {
      target: [(minX + maxX) / 2, (minY + maxY) / 2] as [number, number],
      zoom,
      minZoom: zoom - 4,
      maxZoom: zoom + 6,
    };
  }, [data, containerSize]);

  // No pan clamping — bounds computation disabled for perf testing.

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
    stride: 2,
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

  const views = useMemo(() => new OrthographicView({ id: "ortho" }), []);

  const controller = useMemo(
    () => ({ dragPan: selectMode === "pan", scrollZoom: { speed: 0.01, smooth: false }, inertia: false }),
    [selectMode],
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
        height: "100%",
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
        views={views}
        initialViewState={initialViewState}
        controller={controller}
        useDevicePixels={false}
        layers={layers}
        widgets={debugWidgets}
      />
      <SelectionOverlay selectionRectRef={selectionRectRef} lassoSvgRef={lassoSvgRef} />

      {loading && (
        <div
          style={{
            position: "absolute",
            inset: 0,
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            background: "rgba(255,255,255,0.8)",
            zIndex: 2,
            fontSize: 14,
            color: "#999",
          }}
        >
          Loading {label}...
        </div>
      )}

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
