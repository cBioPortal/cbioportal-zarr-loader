import { useMemo } from "react";
import useAppStore from "../../store/useAppStore";
import { ensureFloat32 } from "../../utils/ensureFloat32";
import EmbeddingScatterplotGL from "../charts/EmbeddingScatterplotGL";
import type { SelectionGeometry } from "../charts/EmbeddingScatterplotGL";

interface Props {
  data: Float32Array | Float64Array;
  shape: [number, number];
  label: string;
  debugMode?: boolean;
}

interface StoreSlice {
  obsmLoading: boolean;
  selectedPointIndices: number[];
  setSelectedPoints: (indices: number[]) => void;
  clearSelectedPoints: () => void;
  selectionGeometry: SelectionGeometry | null;
  setSelectionGeometry: (geo: SelectionGeometry | null) => void;
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
  } = useAppStore() as StoreSlice;

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
