import { Card, Typography } from "antd";
import useAppStore from "../store/useAppStore";

const { Text } = Typography;

export default function LayersPage() {
  const { metadata } = useAppStore();
  const { layerKeys } = metadata;

  return (
    <Card title="Layers" size="small">
      {layerKeys.length ? (
        <div>
          {layerKeys.map((k) => (
            <div key={k} style={{ padding: "4px 0" }}>
              {k}
            </div>
          ))}
        </div>
      ) : (
        <Text type="secondary">(none)</Text>
      )}
    </Card>
  );
}
