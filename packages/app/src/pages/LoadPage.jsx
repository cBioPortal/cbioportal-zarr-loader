import { useState } from "react";
import { useSearchParams, useNavigate } from "react-router";
import { Card, Input, Button, Typography, Space, Alert } from "antd";

const { Title, Text } = Typography;

const isEmbedded = window.self !== window.top || new URLSearchParams(window.location.search).has("embedded");

const DEFAULT_URL = "https://cbioportal-public-imaging.assets.cbioportal.org/msk_spectrum_tme_2022/zarr/spectrum_all_cells.zarr";

export default function LoadPage() {
  const [searchParams] = useSearchParams();
  const navigate = useNavigate();
  const [url, setUrl] = useState(searchParams.get("url") || DEFAULT_URL);

  const handleSubmit = () => {
    const trimmed = url.trim();
    if (!trimmed) return;
    const params = new URLSearchParams(searchParams);
    params.set("url", trimmed);
    navigate(`/?${params.toString()}`);
  };

  if (isEmbedded) {
    return (
      <div style={{ display: "flex", justifyContent: "center", alignItems: "center", minHeight: "60vh", padding: 24 }}>
        <Alert
          type="info"
          message="Load Dataset is not available in embedded mode"
          description="The dataset URL is controlled by the host application."
          showIcon
        />
      </div>
    );
  }

  return (
    <div style={{ display: "flex", justifyContent: "center", alignItems: "center", minHeight: "60vh", padding: 24 }}>
      <Card style={{ maxWidth: 640, width: "100%" }}>
        <Space direction="vertical" size="large" style={{ width: "100%" }}>
          <div>
            <Title level={3} style={{ marginBottom: 4 }}>Load Dataset</Title>
            <Text type="secondary">Enter the URL of a Zarr store to explore.</Text>
          </div>
          <Input
            placeholder="https://example.com/data.zarr"
            value={url}
            onChange={(e) => setUrl(e.target.value)}
            onPressEnter={handleSubmit}
            size="large"
          />
          <Button type="primary" size="large" onClick={handleSubmit} block>
            Load
          </Button>
        </Space>
      </Card>
    </div>
  );
}
