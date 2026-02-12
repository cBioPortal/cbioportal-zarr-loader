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

const ColorBySchema = z.object({
  type: z.enum(["category", "gene"]),
  value: z.string(),
  color_scale: z.enum(["viridis", "magma"]).optional(),
});

const SelectionSchema = z.object({
  target: z.string(),
  values: z.array(z.union([z.string(), z.number()])),
});

const ViewSchema = z.object({
  name: z.string().optional(),
  embedding_key: z.string().optional(),
  selection: SelectionSchema,
  active_tooltips: z.array(z.string()).optional(),
  color_by: ColorBySchema.optional(),
});

const DefaultsSchema = z.object({
  embedding_key: z.string().optional(),
  active_tooltips: z.array(z.string()).optional(),
  color_by: ColorBySchema.optional(),
});

const FilterSchema = z.object({
  defaults: DefaultsSchema.optional(),
  initial_view: z.union([z.string(), z.number().int().nonnegative()]),
  saved_views: z.array(ViewSchema),
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
    setColorColumn,
    setSelectedGene,
    setColorScaleName,
  } = useAppStore();

  const [filterJson, setFilterJson] = useState(JSON.stringify({
    defaults: { embedding_key: "X_umap50", active_tooltips: ["cell_type", "author_sample_id"], color_by: { type: "category", value: "cell_type" } },
    initial_view: "OV-070 by cell_type",
    saved_views: [
      {
        name: "OV-070 by cell_type",
        selection: { target: "donor_id", values: ["SPECTRUM-OV-070"] },
        active_tooltips: ["cell_type", "author_sample_id"],
        color_by: { type: "category", value: "cell_type" }
      },
      {
        name: "OV-090 & OV-022 by cell_type",
        selection: { target: "donor_id", values: ["SPECTRUM-OV-090", "SPECTRUM-OV-022"] },
        active_tooltips: ["cell_type", "author_sample_id", "Phase"],
        color_by: { type: "category", value: "cell_type" }
      },
      {
        selection: { target: "donor_id", values: ["SPECTRUM-OV-041"] },
        active_tooltips: ["cell_type", "author_sample_id", "Phase"],
        color_by: { type: "gene", value: "dapl1", color_scale: "magma" }
      }
    ]
  }, null, 2));
  const [appliedSelections, setAppliedSelections] = useState([]);
  const [activeSelectionIndex, setActiveSelectionIndex] = useState(undefined);
  const [defaults, setDefaults] = useState({});

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

  const applyView = async ({ embedding_key, selection, active_tooltips, color_by }) => {
    // Load embedding: use view-specific key, fall back to default
    const embeddingKey = embedding_key || defaults.embedding_key;
    if (embeddingKey) {
      await fetchObsm(embeddingKey);
    }

    const { target, values } = selection;

    // Clear existing tooltips and load only the ones specified
    clearTooltipColumns();
    const tooltips = active_tooltips || defaults.active_tooltips || [];
    const columnsToLoad = [target, ...tooltips];
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

    // Apply color_by: use view-specific, fall back to default
    const resolvedColorBy = color_by || defaults.color_by;
    if (resolvedColorBy) {
      try {
        if (resolvedColorBy.type === "category") {
          await setColorColumn(resolvedColorBy.value);
        } else if (resolvedColorBy.type === "gene") {
          const { geneNames } = useAppStore.getState().metadata;
          const match = geneNames.find(g => g.toLowerCase() === resolvedColorBy.value.toLowerCase());
          if (!match) {
            throw new Error(`Gene "${resolvedColorBy.value}" not found`);
          }
          await setSelectedGene(match);
        }
        if (resolvedColorBy.color_scale) {
          setColorScaleName(resolvedColorBy.color_scale);
        }
      } catch (err) {
        message.error(`color_by: ${err.message}`);
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

    const result = FilterSchema.safeParse(raw);
    if (!result.success) {
      message.error(result.error.issues.map(i => i.message).join("; "));
      return;
    }

    const { defaults: parsedDefaults = {}, initial_view, saved_views } = result.data;
    setDefaults(parsedDefaults);
    let initialMatch;
    if (typeof initial_view === "number") {
      if (initial_view >= saved_views.length) {
        message.error(`Index ${initial_view} out of range (${saved_views.length} saved views, use 0-${saved_views.length - 1})`);
        return;
      }
      initialMatch = saved_views[initial_view];
    } else {
      initialMatch = saved_views.find(s => s.name === initial_view);
      if (!initialMatch) {
        message.error(`View "${initial_view}" not found in saved_views`);
        return;
      }
    }

    await applyView(initialMatch);

    // Populate the dropdown with all saved_views
    setAppliedSelections(saved_views);
    setActiveSelectionIndex(saved_views.indexOf(initialMatch));
  };

  const handleSelectionPick = async (index) => {
    setActiveSelectionIndex(index);
    await applyView(appliedSelections[index]);
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
                label: s.name || `${s.selection.target}: ${s.selection.values.join(", ")}`,
              }))}
            />
            <Input.TextArea
              autoSize={{ minRows: 2 }}
              value={filterJson}
              onChange={e => setFilterJson(e.target.value)}
              placeholder='{"initial_view": "my view", "saved_views": [{"name": "my view", "target": "donor_id", "values": ["..."]}]}'
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
