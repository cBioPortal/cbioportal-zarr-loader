import { useMemo, useState } from "react";
import { Select, Input } from "antd";
import { useNavigate, useSearchParams } from "react-router";
import useAppStore from "../../store/useAppStore";
import { getRecentUrls } from "../../utils/recentUrls";

const sectionLabelStyle = {
  fontSize: 12,
  fontWeight: 600,
  color: "#666",
  textTransform: "uppercase",
  marginBottom: 8,
};

function filenameFromUrl(url) {
  try {
    const pathname = new URL(url).pathname;
    const segments = pathname.split("/").filter(Boolean);
    return segments[segments.length - 1] || url;
  } catch {
    return url;
  }
}

export default function SidebarDatasetPicker() {
  const { url, isEmbedded } = useAppStore();
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const [pasteValue, setPasteValue] = useState("");

  const recentUrls = useMemo(() => getRecentUrls(), [url]);

  const options = useMemo(
    () =>
      recentUrls.map((entry) => ({
        label: filenameFromUrl(entry.url),
        value: entry.url,
        title: entry.url,
      })),
    [recentUrls],
  );

  if (isEmbedded) return null;

  const navigateWithParams = (newUrl) => {
    const params = new URLSearchParams(searchParams);
    params.set("url", newUrl);
    navigate(`/?${params.toString()}`);
  };

  const handleSelect = (selectedUrl) => {
    if (selectedUrl !== url) {
      navigateWithParams(selectedUrl);
    }
  };

  const handlePaste = () => {
    const trimmed = pasteValue.trim();
    if (trimmed) {
      setPasteValue("");
      navigateWithParams(trimmed);
    }
  };

  return (
    <div style={{ padding: "12px 16px", borderBottom: "1px solid #f0f0f0" }}>
      <div style={sectionLabelStyle}>Dataset</div>
      <Select
        size="small"
        style={{ width: "100%" }}
        value={url}
        onChange={handleSelect}
        options={options}
        optionLabelProp="label"
        placeholder="No recent datasets"
      />
      <Input
        size="small"
        style={{ marginTop: 8 }}
        placeholder="Paste URL..."
        value={pasteValue}
        onChange={(e) => setPasteValue(e.target.value)}
        onPressEnter={handlePaste}
      />
    </div>
  );
}
