import { useState } from "react";
import { Card, Button, Segmented } from "antd";
import useAppStore from "../store/useAppStore";
import SearchableList from "./SearchableList";

export default function ColorByPanel({ height = 300, width = 220, style = {} }) {
  const {
    metadata,
    colorColumn,
    colorLoading,
    setColorColumn,
    selectedGene,
    setSelectedGene,
    clearGeneSelection,
  } = useAppStore();

  const { obsColumns, geneNames } = metadata || {};
  const [activeTab, setActiveTab] = useState("columns");

  const isColumns = activeTab === "columns";
  const selected = isColumns ? colorColumn : selectedGene;
  const onClear = isColumns ? () => setColorColumn(null) : clearGeneSelection;
  const selectedCount = selected ? 1 : 0;

  return (
    <Card
      size="small"
      title={
        <Segmented
          size="small"
          value={activeTab}
          onChange={setActiveTab}
          options={[
            { label: `Columns (${obsColumns?.length ?? 0})`, value: "columns" },
            { label: `Genes (${geneNames?.length ?? 0})`, value: "genes" },
          ]}
        />
      }
      extra={onClear && selectedCount > 0 ? (
        <Button type="link" size="small" onClick={onClear} style={{ padding: 0 }}>
          Clear
        </Button>
      ) : null}
      style={{ width, height, ...style }}
      styles={{
        body: {
          padding: 0,
          height: "calc(100% - 38px)",
          display: "flex",
          flexDirection: "column",
          overflow: "hidden",
        },
      }}
    >
      {isColumns ? (
        <SearchableList
          bare
          items={obsColumns}
          selected={colorColumn}
          onSelect={setColorColumn}
          loading={colorLoading ? colorColumn : null}
          placeholder="Search columns..."
        />
      ) : (
        <SearchableList
          bare
          items={geneNames}
          selected={selectedGene}
          onSelect={setSelectedGene}
          placeholder="Search genes..."
        />
      )}
    </Card>
  );
}
