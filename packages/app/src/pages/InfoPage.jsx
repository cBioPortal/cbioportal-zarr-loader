import { Card, Descriptions, Typography } from "antd";
import useAppStore from "../store/useAppStore";

const { Text } = Typography;

const URL = "https://cbioportal-public-imaging.assets.cbioportal.org/msk_spectrum_tme_2022/zarr/spectrum_all_cells.zarr";

export default function InfoPage() {
  const { adata, metadata } = useAppStore();
  const { chunks } = metadata;

  return (
    <Card title="Dataset" size="small">
      <Descriptions column={1} size="small">
        <Descriptions.Item label="Shape">
          {adata.nObs.toLocaleString()} obs × {adata.nVar.toLocaleString()} var
        </Descriptions.Item>
        <Descriptions.Item label="Chunk size">
          {chunks ? chunks.join(" × ") : "N/A"}
        </Descriptions.Item>
        <Descriptions.Item label="Encoding">
          {adata.attrs["encoding-type"]} v{adata.attrs["encoding-version"]}
        </Descriptions.Item>
        <Descriptions.Item label="URL">
          <Text copyable style={{ fontSize: 12 }}>
            {URL}
          </Text>
        </Descriptions.Item>
      </Descriptions>
    </Card>
  );
}
