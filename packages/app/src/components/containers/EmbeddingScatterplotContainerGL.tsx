import { useMemo } from "react";
import { useShallow } from "zustand/react/shallow";
import useAppStore from "../../store/useAppStore";
import { ensureFloat32 } from "../../utils/ensureFloat32";
import EmbeddingScatterplotGL from "../charts/EmbeddingScatterplotGL";

interface Props {
  data: Float32Array | Float64Array;
  shape: [number, number];
  label: string;
  debugMode?: boolean;
}

export default function EmbeddingScatterplotContainerGL({
  data,
  shape,
  label,
  debugMode = false,
}: Props) {
  const {
    obsmLoading,
    selectedPointIndices,
    setSelectedPoints,
    clearSelectedPoints,
    selectionGeometry,
    setSelectionGeometry,
  } = useAppStore(
    useShallow((s) => ({
      obsmLoading: s.obsmLoading,
      selectedPointIndices: s.selectedPointIndices,
      setSelectedPoints: s.setSelectedPoints,
      clearSelectedPoints: s.clearSelectedPoints,
      selectionGeometry: s.selectionGeometry,
      setSelectionGeometry: s.setSelectionGeometry,
    }))
  );

  const f32Data = useMemo(() => ensureFloat32(data), [data]);

  return (
    <EmbeddingScatterplotGL
      data={f32Data}
      shape={shape}
      label={label}
      loading={obsmLoading}
      selectedPointIndices={selectedPointIndices}
      setSelectedPoints={setSelectedPoints}
      clearSelectedPoints={clearSelectedPoints}
      selectionGeometry={selectionGeometry}
      setSelectionGeometry={setSelectionGeometry}
      debugMode={debugMode}
    />
  );
}
