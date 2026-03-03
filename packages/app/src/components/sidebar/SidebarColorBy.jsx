import { useState, useMemo } from "react";
import { Select } from "antd";
import useAppStore from "../../store/useAppStore";

const DEBOUNCE_MS = 250;
const MAX_GENE_OPTIONS = 100;

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
  const [searchText, setSearchText] = useState("");

  const isColumns = mode === "columns";
  const selected = isColumns ? colorColumn : selectedGene;
  const loading = isColumns ? colorLoading : geneLoading;

  // Columns: show all (typically small). Genes: only show filtered matches
  // to avoid rc-virtual-list scrollTo failure with 20k+ items.
  const options = useMemo(() => {
    if (isColumns) {
      return (obsColumns || []).map((name) => ({ label: name, value: name }));
    }
    if (!geneNames || !searchText) return [];
    const lower = searchText.toLowerCase();
    const matches = [];
    for (const name of geneNames) {
      if (name.toLowerCase().includes(lower)) {
        matches.push({ label: name, value: name });
        if (matches.length >= MAX_GENE_OPTIONS) break;
      }
    }
    return matches;
  }, [isColumns, obsColumns, geneNames, searchText]);

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
    setSearchText("");
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
        onSearch={setSearchText}
        options={options}
        loading={loading}
        disabled={loading}
        placeholder={isColumns ? "Select column..." : "Type to search genes..."}
        filterOption={isColumns}
        optionFilterProp="label"
        notFoundContent={
          !isColumns && !searchText ? "Type to search..." : undefined
        }
        virtual
      />
    </div>
  );
}
