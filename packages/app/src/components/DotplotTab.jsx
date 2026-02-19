import { Card, Tag, Typography, message } from "antd";
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

  const handleGeneSelect = async (geneName) => {
    const result = await toggleDotplotGene(geneName);
    if (result?.noExpression) {
      message.warning(`No expression data found for ${geneName}`);
    } else if (result?.error) {
      message.error(`Failed to fetch expression for ${geneName}`);
    }
  };

  return (
    <TabLayout
      sidebar={
        <SearchableList
          title="Genes"
          items={geneNames}
          selected={dotplotGenes}
          onSelect={handleGeneSelect}
          onClear={clearDotplotGenes}
          loading={dotplotGeneLoading}
          multiSelect
          placeholder="Search genes..."
          height={350}
        />
      }
    >
      <Card
        title={
          <span style={{ display: "flex", alignItems: "center", gap: 4, flexWrap: "wrap" }}>
            Dotplot
            {dotplotGenes.map((gene) => (
              <Tag key={gene} closable onClose={() => toggleDotplotGene(gene)} style={{ marginInlineEnd: 0 }}>
                {gene}
              </Tag>
            ))}
          </span>
        }
        size="small"
      >
        <Text type="secondary">Coming soon</Text>
      </Card>
    </TabLayout>
  );
}
