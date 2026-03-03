import { useState, useMemo } from "react";
import { Select } from "antd";
import useAppStore from "../../store/useAppStore";

const DEBOUNCE_MS = 250;

const sectionLabelStyle = {
  fontSize: 12,
  fontWeight: 600,
  color: "#666",
  textTransform: "uppercase",
};

export default function SidebarColorBy() {
  const {
    metadata,
    colorColumn,
    selectedGene,
    colorLoading,
    geneLoading,
    setColorColumn,
    setSelectedGene,
    clearGeneSelection,
  } = useAppStore();

  const { obsColumns, geneNames } = metadata;

  const [mode, setMode] = useState("columns");

  const isColumns = mode === "columns";
  const items = isColumns ? obsColumns : geneNames;
  const selected = isColumns ? colorColumn : selectedGene;
  const loading = isColumns ? colorLoading : geneLoading;

  const options = useMemo(
    () => (items || []).map((name) => ({ label: name, value: name })),
    [items],
  );

  const debouncedSelect = useMemo(() => {
    let timer;
    return (value) => {
      clearTimeout(timer);
      timer = setTimeout(() => {
        if (isColumns) {
          setColorColumn(value);
        } else {
          setSelectedGene(value);
        }
      }, DEBOUNCE_MS);
    };
  }, [isColumns, setColorColumn, setSelectedGene]);

  const handleModeChange = (newMode) => {
    setMode(newMode);
    if (newMode === "columns") {
      clearGeneSelection();
    } else {
      setColorColumn(null);
    }
  };

  return (
    <div
      style={{
        display: "flex",
        flexDirection: "column",
        padding: "12px 16px",
        borderBottom: "1px solid #f0f0f0",
        gap: 8,
      }}
    >
      {/* Header row */}
      <div
        style={{
          display: "flex",
          alignItems: "center",
          justifyContent: "space-between",
        }}
      >
        <div style={sectionLabelStyle}>Color by</div>
        <Select
          size="small"
          value={mode}
          onChange={handleModeChange}
          style={{ width: 100 }}
          options={[
            { label: "Columns", value: "columns" },
            { label: "Genes", value: "genes" },
          ]}
        />
      </div>

      {/* Searchable dropdown */}
      <Select
        size="small"
        style={{ width: "100%" }}
        showSearch
        value={selected}
        onChange={debouncedSelect}
        options={options}
        loading={loading}
        disabled={loading}
        placeholder={isColumns ? "Select column..." : "Select gene..."}
        optionFilterProp="label"
        virtual
      />
    </div>
  );
}
