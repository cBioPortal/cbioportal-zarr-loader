import { useEffect, useState } from "react";
import {
  Card,
  Drawer,
  Typography,
  Spin,
  Alert,
  Space,
  Button,
  Input,
  Select,
  message,
} from "antd";
import { ReloadOutlined, EditOutlined } from "@ant-design/icons";
import EmbeddingScatterplot from "./EmbeddingScatterplot";
import SearchableList from "./SearchableList";
import ColorColumnList from "./ColorColumnList";
import GeneList from "./GeneList";
import TooltipColumnList from "./TooltipColumnList";
import TabLayout from "./TabLayout";
import useAppStore from "../store/useAppStore";
import {
  FilterSchema,
  findMatchingIndices,
  resolveInitialView,
  resolveViewWithDefaults,
} from "../utils/filterUtils";

const { Text } = Typography;

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
  const [drawerOpen, setDrawerOpen] = useState(false);

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

  const applyView = async (view) => {
    const resolved = resolveViewWithDefaults(view, defaults);

    // Load embedding
    if (resolved.embeddingKey) {
      await fetchObsm(resolved.embeddingKey);
    }

    const { target, values } = resolved.selection;

    // Clear existing tooltips and load only the ones specified
    clearTooltipColumns();
    const columnsToLoad = [target, ...resolved.activeTooltips];
    await Promise.all(columnsToLoad.map(col => toggleTooltipColumn(col)));

    const columnData = useAppStore.getState().tooltipData[target];
    if (!columnData) {
      message.error(`Failed to load column "${target}".`);
      return;
    }

    const matchingIndices = findMatchingIndices(columnData, values);

    // Apply color_by
    if (resolved.colorBy) {
      try {
        if (resolved.colorBy.type === "category") {
          await setColorColumn(resolved.colorBy.value);
        } else if (resolved.colorBy.type === "gene") {
          const { geneNames } = useAppStore.getState().metadata;
          const match = geneNames.find(g => g.toLowerCase() === resolved.colorBy.value.toLowerCase());
          if (!match) {
            throw new Error(`Gene "${resolved.colorBy.value}" not found`);
          }
          await setSelectedGene(match);
        }
        if (resolved.colorBy.color_scale) {
          setColorScaleName(resolved.colorBy.color_scale);
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
    const initialMatch = resolveInitialView(initial_view, saved_views);
    if (!initialMatch) {
      const msg = typeof initial_view === "number"
        ? `Index ${initial_view} out of range (${saved_views.length} saved views, use 0-${saved_views.length - 1})`
        : `View "${initial_view}" not found in saved_views`;
      message.error(msg);
      return;
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
              <Space>
                <Select
                  placeholder="No applied views"
                  size="small"
                  style={{ width: 240 }}
                  onChange={handleSelectionPick}
                  value={activeSelectionIndex}
                  options={appliedSelections.map((s, i) => ({
                    value: i,
                    label: `${i}: ${s.name || `${s.selection.target}: ${s.selection.values.join(", ")}`}`,
                  }))}
                />
                <Button
                  size="small"
                  icon={<EditOutlined />}
                  onClick={() => setDrawerOpen(true)}
                  title="Edit JSON filter"
                />
                <Button
                  size="small"
                  icon={<ReloadOutlined />}
                  onClick={() => fetchObsm(selectedObsm)}
                  loading={obsmLoading}
                />
              </Space>
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
          <Drawer
            title="Edit JSON Filter"
            placement="right"
            width={480}
            open={drawerOpen}
            onClose={() => setDrawerOpen(false)}
          >
            <Input.TextArea
              autoSize={{ minRows: 2 }}
              value={filterJson}
              onChange={e => setFilterJson(e.target.value)}
              placeholder='{"initial_view": "my view", "saved_views": [{"name": "my view", ...}]}'
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
          </Drawer>
        </>
      ) : (
        <Card size="small">
          <Text type="secondary">Select a key to view its data</Text>
        </Card>
      )}
    </TabLayout>
  );
}
