import { Typography, Space } from "antd";
import EmbeddingScatterplot from "./EmbeddingScatterplot";
import useAppStore from "../../store/useAppStore";

const { Text } = Typography;

export default function StreamingEmbeddingScatterplot({ label }) {
  const { obsmStreamingData, obsmStreamingProgress, obsmStreamingTime } =
    useAppStore();

  if (!obsmStreamingData) return null;

  return (
    <>
      <Space style={{ marginBottom: 8 }} wrap>
        <Text>Shape: {obsmStreamingData.shape.join(" × ")}</Text>
        {obsmStreamingProgress !== null ? (
          <Text type="secondary">
            Loading: {Math.round(obsmStreamingProgress * 100)}%
          </Text>
        ) : (
          <Text>Fetched in {obsmStreamingTime?.toFixed(1)} ms</Text>
        )}
      </Space>
      <EmbeddingScatterplot
        data={obsmStreamingData.data}
        shape={obsmStreamingData.shape}
        label={label}
      />
    </>
  );
}
