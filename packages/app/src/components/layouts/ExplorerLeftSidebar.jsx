import { useMemo } from "react";
import { Button } from "antd";
import { GithubOutlined, UploadOutlined } from "@ant-design/icons";
import { Link } from "react-router";
import useAppStore from "../../store/useAppStore";
import useLinkWithParams from "../../hooks/useLinkWithParams";
import SearchableList from "../ui/SearchableList";
import ColorByPanel from "../ui/ColorByPanel";

export default function ExplorerLeftSidebar() {
  const { metadata, featureFlags, selectedObsm, obsmLoading, fetchObsm } = useAppStore();
  const linkTo = useLinkWithParams();
  const { obsmKeys } = metadata;

  const debouncedFetchObsm = useMemo(() => {
    let timer;
    return (key) => {
      clearTimeout(timer);
      useAppStore.setState({ selectedObsm: key });
      timer = setTimeout(() => fetchObsm(key), 250);
    };
  }, [fetchObsm]);

  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <div
        style={{
          padding: "12px 16px",
          borderBottom: "1px solid #f0f0f0",
          display: "flex",
          alignItems: "center",
          justifyContent: "space-between",
        }}
      >
        <Link
          to={linkTo("/")}
          style={{
            fontSize: 14,
            fontWeight: 600,
            color: "inherit",
            textDecoration: "none",
          }}
        >
          ZExplorer
        </Link>
        <nav style={{ display: "flex", gap: 8, alignItems: "center" }}>
          {featureFlags.loadDataset && (
            <Link to={linkTo("/load")}>
              <Button type="text" size="small" icon={<UploadOutlined />} />
            </Link>
          )}
          <a
            href="https://github.com/cbioportal/cbioportal-zarr-loader"
            target="_blank"
            rel="noopener noreferrer"
          >
            <GithubOutlined style={{ fontSize: 16, color: "#666" }} />
          </a>
        </nav>
      </div>
      <div style={{ flex: 1, padding: 16, overflow: "auto" }}>
        <SearchableList
          title="Keys"
          items={obsmKeys}
          selected={selectedObsm}
          onSelect={debouncedFetchObsm}
          loading={obsmLoading ? selectedObsm : null}
          placeholder="Search keys..."
          height={200}
        />
        <ColorByPanel height={300} style={{ marginTop: 16 }} />
      </div>
    </div>
  );
}
