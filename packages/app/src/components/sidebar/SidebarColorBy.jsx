import { useState, useMemo } from "react";
import { Select } from "antd";
import useAppStore from "../../store/useAppStore";

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
    setColorColumn,
    setSelectedGene,
    clearGeneSelection,
  } = useAppStore();

  const { obsColumns, geneNames } = metadata;

  const [mode, setMode] = useState("columns");

  const isColumns = mode === "columns";
  const items = isColumns ? obsColumns : geneNames;
  const selected = isColumns ? colorColumn : selectedGene;

  const options = useMemo(
    () => (items || []).map((name) => ({ label: name, value: name })),
    [items],
  );

  const handleModeChange = (newMode) => {
    setMode(newMode);
    if (newMode === "columns") {
      clearGeneSelection();
    } else {
      setColorColumn(null);
    }
  };

  const handleSelect = (value) => {
    if (isColumns) {
      setColorColumn(value);
    } else {
      setSelectedGene(value);
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
        onChange={handleSelect}
        options={options}
        loading={colorLoading}
        placeholder={isColumns ? "Select column..." : "Select gene..."}
        optionFilterProp="label"
        virtual
      />
    </div>
  );
}
