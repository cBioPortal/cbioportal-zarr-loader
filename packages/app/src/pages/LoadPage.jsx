import { useState, useEffect } from "react";
import { useSearchParams, useNavigate } from "react-router";
import { Card, Input, Button, Typography, Space, Alert, List } from "antd";
import { DeleteOutlined, ClearOutlined, PushpinFilled } from "@ant-design/icons";
import { getRecentUrls, removeRecentUrl, clearRecentUrls } from "../utils/recentUrls";
import { DEFAULT_URL } from "../constants";

const { Title, Text } = Typography;

const isEmbedded = window.self !== window.top || new URLSearchParams(window.location.search).has("embedded");

export default function LoadPage() {
  const [searchParams] = useSearchParams();
  const navigate = useNavigate();
  const [url, setUrl] = useState(searchParams.get("url") || "");
  const [recentUrls, setRecentUrls] = useState([]);

  useEffect(() => {
    setRecentUrls(getRecentUrls());
  }, []);

  const handleSubmit = () => {
    const trimmed = url.trim();
    if (!trimmed) return;
    const params = new URLSearchParams(window.location.search);
    params.set("url", trimmed);
    navigate(`/?${decodeURIComponent(params.toString())}`);
  };

  const handleRemove = (urlToRemove) => {
    removeRecentUrl(urlToRemove);
    setRecentUrls(getRecentUrls());
  };

  const handleClearAll = () => {
    clearRecentUrls();
    setRecentUrls([]);
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

  const datasetList = [
    ...(DEFAULT_URL ? [{ url: DEFAULT_URL, pinned: true }] : []),
    ...recentUrls.filter((entry) => entry.url !== DEFAULT_URL),
  ];

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
          {datasetList.length > 0 && (
            <div>
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 8 }}>
                <Text strong>Datasets</Text>
                {recentUrls.length > 0 && (
                  <Button
                    type="link"
                    size="small"
                    icon={<ClearOutlined />}
                    onClick={handleClearAll}
                  >
                    Clear recents
                  </Button>
                )}
              </div>
              <List
                size="small"
                bordered
                dataSource={datasetList}
                renderItem={(entry) => (
                  <List.Item
                    style={{ cursor: "pointer" }}
                    onClick={() => setUrl(entry.url)}
                    actions={
                      entry.pinned
                        ? []
                        : [
                            <Button
                              key="delete"
                              type="text"
                              size="small"
                              icon={<DeleteOutlined />}
                              onClick={(e) => {
                                e.stopPropagation();
                                handleRemove(entry.url);
                              }}
                            />,
                          ]
                    }
                  >
                    <div style={{ minWidth: 0, flex: 1 }}>
                      {entry.pinned && (
                        <Text type="secondary" style={{ fontSize: 12, display: "block" }}>
                          <PushpinFilled style={{ marginRight: 4 }} />
                          Default
                        </Text>
                      )}
                      <Text
                        ellipsis={{ tooltip: entry.url }}
                        style={{ maxWidth: "100%" }}
                      >
                        {entry.url}
                      </Text>
                    </div>
                  </List.Item>
                )}
              />
            </div>
          )}
        </Space>
      </Card>
    </div>
  );
}
