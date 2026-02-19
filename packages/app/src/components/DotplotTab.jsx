import { useEffect, useMemo, useRef } from "react";
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

  // Auto-select defaults on first mount
  const initialized = useRef(false);
  useEffect(() => {
    if (initialized.current) return;
    initialized.current = true;

    const defaultGenes = ["EGFR", "DAPL1"];
    for (const name of defaultGenes) {
      const match = geneNames.find((g) => g.toLowerCase() === name.toLowerCase());
      if (match && !dotplotGenes.includes(match)) toggleDotplotGene(match);
    }

    const obsMatch = obsColumns.find((c) => c.toLowerCase() === "cell_type");
    if (obsMatch && !dotplotObsColumn) setDotplotObsColumn(obsMatch);
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

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

  // Size chart based on data dimensions (groups as integers on x, genes on y)
  const chartWidth = groups.length * 20 + 136;
  const chartHeight = dotplotGenes.length * 28 + 56;

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
          <>
            <div style={{ marginBottom: 8, fontSize: 12, color: "#595959" }}>
              {dotplotObsData.length.toLocaleString()} cells across {groups.length} groups
              {dotplotObsColumn && <> grouped by <strong>{dotplotObsColumn}</strong></>}
            </div>
            <div style={{ overflowX: "auto" }}>
              <Dotplot
                genes={dotplotGenes}
                groups={groups}
                data={dotplotData}
                width={chartWidth}
                height={chartHeight}
              />
            </div>
            <div style={{ marginTop: 12, fontSize: 11, color: "#595959", lineHeight: 1.6 }}>
              <div style={{ marginBottom: 8 }}>
                <strong>X-axis legend</strong>
                <div style={{ display: "flex", flexWrap: "wrap", gap: "2px 12px", marginTop: 2 }}>
                  {groups.map((group, i) => (
                    <span key={group}><strong>{i + 1}</strong> = {group}</span>
                  ))}
                </div>
              </div>
              <div style={{ color: "#8c8c8c" }}>
                <p style={{ margin: "0 0 6px" }}>
                  A dotplot summarizes gene expression across cell groups from single-cell RNA-seq data.
                  Each dot represents one gene in one cell group, answering two questions at a glance:
                </p>
                <ul style={{ margin: "0 0 6px", paddingLeft: 18 }}>
                  <li>
                    <strong>Dot size (Fraction Expressing)</strong> — the percentage of cells in the
                    group where the gene is detected (expression &gt; 0). A larger dot means
                    the gene is active in more cells.
                  </li>
                  <li>
                    <strong>Dot color (Mean Expression)</strong> — the average expression level
                    across <em>all</em> cells in the group, including non-expressing cells.
                    Mapped to a viridis color scale (dark = low, bright = high).
                  </li>
                </ul>
                <p style={{ margin: 0 }}>
                  Hover a dot for exact values. Hover an x-axis number for the group name.
                </p>
              </div>
            </div>
          </>
        ) : (
          <Text type="secondary">
            Select one or more genes and an obs column to view the dotplot.
          </Text>
        )}
      </Card>
    </TabLayout>
  );
}
