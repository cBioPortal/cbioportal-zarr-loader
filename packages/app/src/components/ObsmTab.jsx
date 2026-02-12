import { useEffect, useState } from "react";
import {
  Card,
  Typography,
  Spin,
  Alert,
  Space,
  Button,
  Input,
  Select,
  message,
} from "antd";
import { ReloadOutlined } from "@ant-design/icons";
import EmbeddingScatterplot from "./EmbeddingScatterplot";
import SearchableList from "./SearchableList";
import ColorColumnList from "./ColorColumnList";
import GeneList from "./GeneList";
import TooltipColumnList from "./TooltipColumnList";
import TabLayout from "./TabLayout";
import { z } from "zod";
import useAppStore from "../store/useAppStore";

const { Text } = Typography;

const selectionSchema = z.object({
  target: z.string(),
  values: z.array(z.union([z.string(), z.number()])),
  active_tooltips: z.array(z.string()).optional(),
});

const filterSchema = z.object({
  selection: selectionSchema,
  saved_selections: z.array(selectionSchema).optional(),
});

export default function ObsmTab() {
  const {
    metadata,
    selectedObsm,
    obsmData,
    obsmLoading,
    obsmTime,
    fetchObsm,
    tooltipData,
    tooltipColumns,
    toggleTooltipColumn,
    clearTooltipColumns,
    setSelectedPoints,
  } = useAppStore();

  const [filterJson, setFilterJson] = useState(JSON.stringify({
    selection: {
      target: "donor_id",
      values: ["SPECTRUM-OV-070"],
      active_tooltips: ["cell_type", "author_sample_id"]
    },
    saved_selections: [
      {
        target: "donor_id",
        values: ["SPECTRUM-OV-090", "SPECTRUM-OV-022"],
        active_tooltips: ["cell_type", "author_sample_id", "Phase"]
      }
    ]
  }, null, 2));
  const [appliedSelections, setAppliedSelections] = useState([]);
  const [activeSelectionIndex, setActiveSelectionIndex] = useState(undefined);

  const { obsmKeys } = metadata;
  const isEmbedding = selectedObsm && /umap|tsne|pca/i.test(selectedObsm) && obsmData?.shape?.[1] >= 2;

  // Auto-fetch UMAP embedding on mount
  useEffect(() => {
    if (!selectedObsm && obsmKeys.length > 0) {
      const umapKey = obsmKeys.find(k => /umap/i.test(k));
      if (umapKey) {
        fetchObsm(umapKey);
      }
    }
  }, [obsmKeys, selectedObsm, fetchObsm]);

  const applySelection = async ({ target, values, active_tooltips = [] }) => {
    // Clear existing tooltips and load only the ones specified
    clearTooltipColumns();
    const columnsToLoad = [target, ...active_tooltips];
    await Promise.all(columnsToLoad.map(col => toggleTooltipColumn(col)));

    const columnData = useAppStore.getState().tooltipData[target];
    if (!columnData) {
      message.error(`Failed to load column "${target}".`);
      return;
    }

    const valueSet = new Set(values.map(String));
    const matchingIndices = [];
    for (let i = 0; i < columnData.length; i++) {
      if (valueSet.has(String(columnData[i]))) {
        matchingIndices.push(i);
      }
    }

    setSelectedPoints(matchingIndices);
    message.success(`Selected ${matchingIndices.length} points`);
  };

  const handleFilterApply = async () => {
    let raw;
    try {
      raw = JSON.parse(filterJson);
    } catch {
      message.error("Invalid JSON");
      return;
    }

    const result = filterSchema.safeParse(raw);
    if (!result.success) {
      message.error(result.error.issues.map(i => i.message).join("; "));
      return;
    }

    await applySelection(result.data.selection);

    // Save the applied selection and any saved_selections to the dropdown history
    const selectionKey = JSON.stringify(result.data.selection);
    const newSelections = [result.data.selection, ...(result.data.saved_selections ?? [])];
    setAppliedSelections(prev => {
      const existing = new Set(prev.map(s => JSON.stringify(s)));
      const toAdd = newSelections.filter(s => !existing.has(JSON.stringify(s)));
      const next = toAdd.length > 0 ? [...prev, ...toAdd] : prev;
      setActiveSelectionIndex(next.findIndex(s => JSON.stringify(s) === selectionKey));
      return next;
    });
  };

  const handleSelectionPick = async (index) => {
    setActiveSelectionIndex(index);
    await applySelection(appliedSelections[index]);
  };

  return (
    <TabLayout
      sidebar={
        <>
          <SearchableList
            title="Keys"
            items={obsmKeys}
            selected={selectedObsm}
            onSelect={fetchObsm}
            loading={obsmLoading ? selectedObsm : null}
            placeholder="Search keys..."
            height={200}
          />
          <ColorColumnList height={250} style={{ marginTop: 16 }} />
          <GeneList height={300} style={{ marginTop: 16 }} />
          <TooltipColumnList height={250} style={{ marginTop: 16 }} />
        </>
      }
    >
      {selectedObsm ? (
        <>
          <Card
            title={`obsm: ${selectedObsm}`}
            size="small"
            extra={
              <Button
                size="small"
                icon={<ReloadOutlined />}
                onClick={() => fetchObsm(selectedObsm)}
                loading={obsmLoading}
              />
            }
          >
            {obsmLoading ? (
              <Spin />
            ) : obsmData?.error ? (
              <Alert type="error" message={obsmData.error} />
            ) : (
              <>
                <Space style={{ marginBottom: 16 }} wrap>
                  <Text>Fetched in {obsmTime?.toFixed(1)} ms</Text>
                  <Text type="secondary">|</Text>
                  <Text>Shape: {obsmData?.shape?.join(" × ")}</Text>
                </Space>

                {isEmbedding && (
                  <EmbeddingScatterplot
                    data={obsmData.data}
                    shape={obsmData.shape}
                    label={selectedObsm}
                  />
                )}
              </>
            )}
          </Card>
          <Card title="Selection Filter" size="small" style={{ marginTop: 16 }}>
            <Select
              placeholder="No applied selections"
              style={{ width: "100%", marginBottom: 8 }}
              onChange={handleSelectionPick}
              value={activeSelectionIndex}
              options={appliedSelections.map((s, i) => ({
                value: i,
                label: `${s.target}: ${s.values.join(", ")}`,
              }))}
            />
            <Input.TextArea
              autoSize={{ minRows: 2 }}
              value={filterJson}
              onChange={e => setFilterJson(e.target.value)}
              placeholder='{"selection": {"target": "donor_id", "values": ["SPECTRUM-OV-070"]}}'
            />
            <Space style={{ marginTop: 8 }}>
              <Button
                size="small"
                onClick={() => {
                  try {
                    setFilterJson(JSON.stringify(JSON.parse(filterJson), null, 2));
                  } catch {
                    message.error("Invalid JSON — cannot format");
                  }
                }}
              >
                Format
              </Button>
              <Button
                type="primary"
                size="small"
                onClick={handleFilterApply}
              >
                Apply Filter
              </Button>
            </Space>
          </Card>
        </>
      ) : (
        <Card size="small">
          <Text type="secondary">Select a key to view its data</Text>
        </Card>
      )}
    </TabLayout>
  );
}
