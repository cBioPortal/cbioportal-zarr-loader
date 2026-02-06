import { useState, useMemo } from "react";
import { Card, Input } from "antd";
import useAppStore from "../store/useAppStore";

/**
 * Gene selection list with search functionality.
 * Reads gene names and selection state from the store.
 */
export default function GeneList({ height = 300, width = 220, style = {} }) {
  const {
    metadata,
    selectedGene,
    setSelectedGene,
  } = useAppStore();

  const { geneNames } = metadata || {};
  const [searchText, setSearchText] = useState("");

  const filteredGenes = useMemo(() => {
    if (!geneNames) return [];
    if (!searchText) return geneNames;
    const search = searchText.toLowerCase();
    return geneNames.filter(name => name.toLowerCase().includes(search));
  }, [geneNames, searchText]);

  if (!geneNames || geneNames.length === 0) return null;

  return (
    <Card
      size="small"
      title={`Genes (${geneNames.length.toLocaleString()})`}
      style={{ width, height, ...style }}
      styles={{ body: { padding: 0, height: "calc(100% - 38px)", display: "flex", flexDirection: "column", overflow: "hidden" } }}
    >
      <div style={{ padding: 8, borderBottom: "1px solid #f0f0f0" }}>
        <Input.Search
          placeholder="Search genes..."
          size="small"
          value={searchText}
          onChange={(e) => setSearchText(e.target.value)}
          allowClear
        />
      </div>
      <div style={{ flex: 1, overflow: "auto" }}>
        {filteredGenes.map((gene) => (
          <div
            key={gene}
            style={{
              padding: "4px 12px",
              cursor: "pointer",
              backgroundColor: selectedGene === gene ? "#e6f4ff" : undefined,
              fontSize: 12,
              whiteSpace: "nowrap",
              overflow: "hidden",
              textOverflow: "ellipsis",
            }}
            onClick={() => setSelectedGene(gene)}
          >
            {gene}
          </div>
        ))}
        {filteredGenes.length === 0 && searchText && (
          <div style={{ padding: "8px 12px", color: "#999", fontSize: 12 }}>
            No genes found
          </div>
        )}
      </div>
      {searchText && (
        <div style={{ padding: "4px 12px", fontSize: 11, color: "#999", borderTop: "1px solid #f0f0f0" }}>
          {filteredGenes.length.toLocaleString()} matches
        </div>
      )}
    </Card>
  );
}
