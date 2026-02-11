import {
  Card,
  Typography,
  Spin,
  Alert,
  Space,
  Button,
} from "antd";
import { ReloadOutlined } from "@ant-design/icons";
import EmbeddingScatterplot from "../components/visualizations/EmbeddingScatterplot";
import ColorColumnList from "../components/selectors/ColorColumnList";
import GeneList from "../components/selectors/GeneList";
import TooltipColumnList from "../components/selectors/TooltipColumnList";
import TabLayout from "../components/ui/TabLayout";

const { Text } = Typography;

export default function EmbeddingExplorerView({
  selectedObsm,
  obsmData,
  obsmLoading,
  obsmTime,
  onReload,
  sidebar,
}) {
  const isEmbedding = selectedObsm && /umap|tsne|pca/i.test(selectedObsm) && obsmData?.shape?.[1] >= 2;

  return (
    <TabLayout
      sidebar={
        <>
          {sidebar}
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
              onClick={onReload}
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
      ) : (
        <Card size="small">
          <Text type="secondary">Select a key to view its data</Text>
        </Card>
      )}
    </TabLayout>
  );
}
