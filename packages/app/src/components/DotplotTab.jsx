import { useMemo } from "react";
import { Card, Tag, Typography, Spin, message } from "antd";
import SearchableList from "./SearchableList";
import TabLayout from "./TabLayout";
import Dotplot from "./charts/Dotplot";
import useAppStore from "../store/useAppStore";

const { Text } = Typography;

export default function DotplotTab() {
  const {
    metadata,
    dotplotGenes,
    dotplotGeneExpressions,
    dotplotGeneLoading,
    dotplotObsColumn,
    dotplotObsData,
    dotplotObsLoading,
    toggleDotplotGene,
    clearDotplotGenes,
    setDotplotObsColumn,
    clearDotplotObsColumn,
  } = useAppStore();

  const { geneNames, obsColumns } = metadata;

  const handleGeneSelect = async (geneName) => {
    const result = await toggleDotplotGene(geneName);
    if (result?.noExpression) {
      message.warning(`No expression data found for ${geneName}`);
    } else if (result?.error) {
      message.error(`Failed to fetch expression for ${geneName}`);
    }
  };

  // Compute unique groups from obs data
  const groups = useMemo(() => {
    if (!dotplotObsData) return [];
    return [...new Set(dotplotObsData)].sort();
  }, [dotplotObsData]);

  // Compute dotplot stats: for each gene × group, calculate fraction expressing and mean expression
  const dotplotData = useMemo(() => {
    if (dotplotGenes.length === 0 || !dotplotObsData || groups.length === 0) return null;

    // Build group → indices mapping once
    const groupIndices = {};
    for (const g of groups) groupIndices[g] = [];
    for (let i = 0; i < dotplotObsData.length; i++) {
      const g = dotplotObsData[i];
      if (groupIndices[g]) groupIndices[g].push(i);
    }

    const stats = [];
    for (const gene of dotplotGenes) {
      const expr = dotplotGeneExpressions[gene];
      if (!expr) continue;
      for (const group of groups) {
        const indices = groupIndices[group];
        if (indices.length === 0) continue;
        let sum = 0;
        let expressing = 0;
        for (const idx of indices) {
          const val = expr[idx];
          sum += val;
          if (val > 0) expressing++;
        }
        stats.push({
          gene,
          group,
          meanExpression: sum / indices.length,
          fractionExpressing: expressing / indices.length,
        });
      }
    }
    return stats;
  }, [dotplotGenes, dotplotGeneExpressions, dotplotObsData, groups]);

  const isLoading = dotplotGeneLoading || dotplotObsLoading;

  // Size chart based on data dimensions
  const chartWidth = Math.max(400, dotplotGenes.length * 60 + 200);
  const chartHeight = Math.max(300, groups.length * 30 + 120);

  return (
    <TabLayout
      sidebar={
        <>
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
          <SearchableList
            title="Obs Columns"
            items={obsColumns}
            selected={dotplotObsColumn}
            onSelect={setDotplotObsColumn}
            onClear={clearDotplotObsColumn}
            loading={dotplotObsLoading ? dotplotObsColumn : null}
            placeholder="Search obs columns..."
            height={350}
            style={{ marginTop: 16 }}
          />
        </>
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
        {isLoading ? (
          <div style={{ textAlign: "center", padding: 48 }}>
            <Spin size="large" />
          </div>
        ) : dotplotData ? (
          <div style={{ overflowX: "auto" }}>
            <Dotplot
              genes={dotplotGenes}
              groups={groups}
              data={dotplotData}
              width={chartWidth}
              height={chartHeight}
            />
          </div>
        ) : (
          <Text type="secondary">
            Select one or more genes and an obs column to view the dotplot.
          </Text>
        )}
      </Card>
    </TabLayout>
  );
}
