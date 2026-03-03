import { useMemo } from "react";
import { Select } from "antd";
import useAppStore from "../../store/useAppStore";

const sectionLabelStyle = {
  fontSize: 12,
  fontWeight: 600,
  color: "#666",
  textTransform: "uppercase",
  marginBottom: 8,
};

export default function SidebarEmbeddingPicker() {
  const { metadata, selectedObsm, obsmLoading, fetchObsm } = useAppStore();
  const { obsmKeys } = metadata;

  const debouncedFetchObsm = useMemo(() => {
    let timer;
    return (key) => {
      clearTimeout(timer);
      useAppStore.setState({ selectedObsm: key });
      timer = setTimeout(() => fetchObsm(key), 250);
    };
  }, [fetchObsm]);

  const options = useMemo(
    () => obsmKeys.map((key) => ({ label: key, value: key })),
    [obsmKeys],
  );

  return (
    <div style={{ padding: "12px 16px", borderBottom: "1px solid #f0f0f0" }}>
      <div style={sectionLabelStyle}>Embedding</div>
      <Select
        size="small"
        style={{ width: "100%" }}
        value={selectedObsm}
        onChange={debouncedFetchObsm}
        options={options}
        loading={obsmLoading}
        showSearch
        placeholder="Select embedding..."
      />
    </div>
  );
}
