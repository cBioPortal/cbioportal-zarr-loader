import { Button } from "antd";
import { DashboardOutlined, GithubOutlined, UploadOutlined } from "@ant-design/icons";
import { Link } from "react-router";
import useAppStore from "../../store/useAppStore";
import useLinkWithParams from "../../hooks/useLinkWithParams";
import SidebarDatasetPicker from "../sidebar/SidebarDatasetPicker";
import SidebarEmbeddingPicker from "../sidebar/SidebarEmbeddingPicker";
import SidebarColorBy from "../sidebar/SidebarColorBy";

export default function ExplorerLeftSidebar() {
  const { featureFlags } = useAppStore();
  const linkTo = useLinkWithParams();

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
          {featureFlags.profile && (
            <Link to={linkTo("/profile")}>
              <Button type="text" size="small" icon={<DashboardOutlined />} />
            </Link>
          )}
          <a
            href="https://github.com/cbioportal/cbioportal-cell-explorer"
            target="_blank"
            rel="noopener noreferrer"
          >
            <GithubOutlined style={{ fontSize: 16, color: "#666" }} />
          </a>
        </nav>
      </div>
      <SidebarDatasetPicker />
      <SidebarEmbeddingPicker />
      <SidebarColorBy />
    </div>
  );
}
