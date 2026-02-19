import { Card, Typography } from "antd";
import SearchableList from "./SearchableList";
import TabLayout from "./TabLayout";
import useAppStore from "../store/useAppStore";

const { Text } = Typography;

export default function DotplotTab() {
  const {
    metadata,
    dotplotGenes,
    dotplotGeneLoading,
    toggleDotplotGene,
    clearDotplotGenes,
  } = useAppStore();

  const { geneNames } = metadata;

  return (
    <TabLayout
      sidebar={
        <SearchableList
          title="Genes"
          items={geneNames}
          selected={dotplotGenes}
          onSelect={toggleDotplotGene}
          onClear={clearDotplotGenes}
          loading={dotplotGeneLoading}
          multiSelect
          placeholder="Search genes..."
          height={350}
        />
      }
    >
      <Card title="Dotplot" size="small">
        <Text type="secondary">Coming soon</Text>
      </Card>
    </TabLayout>
  );
}
