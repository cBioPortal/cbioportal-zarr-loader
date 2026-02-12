import { useEffect, useState } from "react";
import {
  Card,
  Typography,
  Spin,
  Alert,
  Space,
  Button,
  Input,
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

const filterSchema = z.object({
  obs_category: z.string(),
  values: z.array(z.union([z.string(), z.number()])),
  tooltips: z.array(z.string()).optional(),
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
    setSelectedPoints,
  } = useAppStore();

  const [filterJson, setFilterJson] = useState(JSON.stringify({ obs_category: "donor_id", values: ["SPECTRUM-OV-070"], tooltips: ["cell_type", "author_sample_id"] }, null, 2));

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

    const { obs_category, values, tooltips = [] } = result.data;

    // Load the obs_category and any extra tooltip columns if not already loaded
    const columnsToLoad = [obs_category, ...tooltips].filter(
      col => !tooltipColumns.includes(col)
    );
    await Promise.all(columnsToLoad.map(col => toggleTooltipColumn(col)));

    const columnData = useAppStore.getState().tooltipData[obs_category];
    if (!columnData) {
      message.error(`Failed to load column "${obs_category}".`);
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
            <Input.TextArea
              rows={4}
              value={filterJson}
              onChange={e => setFilterJson(e.target.value)}
              placeholder='{"obs_category": "donor_id", "values": ["SPECTRUM-OV-070"]}'
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
