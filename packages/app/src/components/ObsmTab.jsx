import { useEffect } from "react";
import {
  Card,
  Typography,
  Spin,
  Alert,
  Space,
  Button,
} from "antd";
import { ReloadOutlined } from "@ant-design/icons";
import EmbeddingScatterplot from "./EmbeddingScatterplot";
import SearchableList from "./SearchableList";
import ColorColumnList from "./ColorColumnList";
import GeneList from "./GeneList";
import TooltipColumnList from "./TooltipColumnList";
import TabLayout from "./TabLayout";
import useAppStore from "../store/useAppStore";

const { Text } = Typography;

export default function ObsmTab() {
  const {
    metadata,
    selectedObsm,
    obsmData,
    obsmLoading,
    obsmTime,
    fetchObsm,
  } = useAppStore();

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
        <Card
          title={`obsm: ${selectedObsm}`}
          size="small"
          extra={
            <Button
              size="small"
              icon={<ReloadOutlined />}
              onClick={() => fetchObsm(selectedObsm)}
              loading={obsmLoading}
            >
              Reload
            </Button>
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
      ) : (
        <Card size="small">
          <Text type="secondary">Select a key to view its data</Text>
        </Card>
      )}
    </TabLayout>
  );
}
