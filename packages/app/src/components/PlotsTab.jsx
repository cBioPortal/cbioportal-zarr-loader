import { useCallback, useEffect, useMemo, useState } from "react";
import { Card, Typography, Spin, Select } from "antd";
import SearchableList from "./SearchableList";
import TabLayout from "./TabLayout";
import useAppStore from "../store/useAppStore";
import { computeBoxplotStats } from "../utils/boxplotUtils";
import { computeViolinStats } from "../utils/violinUtils";
import ViolinPlot from "./charts/ViolinPlot";

const { Text } = Typography;

export default function PlotsTab() {
  const {
    metadata,
    plotGene,
    plotGeneExpression,
    plotGeneLoading,
    plotObsColumn,
    plotObsData,
    plotObsLoading,
    setPlotGene,
    clearPlotGene,
    setPlotObsColumn,
    clearPlotObsColumn,
  } = useAppStore();

  const { geneNames, obsColumns } = metadata;
  const [filterExpression, setFilterExpression] = useState(null);
  const [containerWidth, setContainerWidth] = useState(800);
  const containerRef = useCallback((node) => {
    if (!node) return;
    const ro = new ResizeObserver(([entry]) => {
      setContainerWidth(entry.contentRect.width);
    });
    ro.observe(node);
    setContainerWidth(node.clientWidth);
  }, []);

  // Dev defaults: auto-select gene and obs column on first mount
  useEffect(() => {
    if (!plotGene && geneNames?.includes("CETN2")) setPlotGene("CETN2");
    if (!plotObsColumn && obsColumns?.includes("cell_type")) setPlotObsColumn("cell_type");
  }, []);

  // Compute top frequent rounded expression values for the exclude dropdown
  const frequentValues = useMemo(() => {
    if (!plotGeneExpression) return [];
    const counts = new Map();
    for (let i = 0; i < plotGeneExpression.length; i++) {
      const v = Math.round(plotGeneExpression[i] * 10000) / 10000;
      counts.set(v, (counts.get(v) || 0) + 1);
    }
    return [...counts.entries()]
      .sort((a, b) => b[1] - a[1])
      .slice(0, 20)
      .map(([val, count]) => ({
        value: val,
        label: `${val} (${count.toLocaleString()}x)`,
      }));
  }, [plotGeneExpression]);

  // Auto-select the most frequent expression value for filtering
  useEffect(() => {
    if (frequentValues.length > 0) {
      setFilterExpression(frequentValues[0].value);
    } else {
      setFilterExpression(null);
    }
  }, [frequentValues]);

  const data = useMemo(() => {
    if (!plotGeneExpression || !plotObsData) return null;
    const raw = Array.from(plotGeneExpression, (val, i) => ({
      [plotObsColumn]: String(plotObsData[i]),
      [plotGene]: Math.round(val * 10000) / 10000,
    }));
    if (filterExpression !== null) {
      return raw.filter((d) => d[plotGene] !== filterExpression);
    }
    return raw;
  }, [plotGeneExpression, plotObsData, plotObsColumn, plotGene, filterExpression]);

  const boxplotData = useMemo(() => {
    if (!data) return null;
    return computeBoxplotStats(data, plotObsColumn, plotGene);
  }, [data, plotObsColumn, plotGene]);

  const violinData = useMemo(() => {
    if (!data) return null;
    return computeViolinStats(data, plotObsColumn, plotGene);
  }, [data, plotObsColumn, plotGene]);

  if (data) {
    const categories = new Set(data.map((d) => d[plotObsColumn]));
    console.debug("[PlotsTab] Data built:", data.length, "points,", categories.size, "categories:", [...categories].slice(0, 10));
    console.debug("[PlotsTab] Sample data:", data.slice(0, 15));
  }

  const isLoading = plotGeneLoading || plotObsLoading;
  const hasSelections = plotGene && plotObsColumn;

  return (
    <TabLayout
      sidebar={
        <>
          <SearchableList
            title="Genes"
            items={geneNames}
            selected={plotGene}
            onSelect={setPlotGene}
            onClear={clearPlotGene}
            loading={plotGeneLoading ? plotGene : null}
            placeholder="Search genes..."
            height={350}
          />
          <SearchableList
            title="Obs Columns"
            items={obsColumns}
            selected={plotObsColumn}
            onSelect={setPlotObsColumn}
            onClear={clearPlotObsColumn}
            loading={plotObsLoading ? plotObsColumn : null}
            placeholder="Search obs columns..."
            height={350}
            style={{ marginTop: 16 }}
          />
        </>
      }
    >
      <Card
        title={hasSelections ? `${plotGene} by ${plotObsColumn}` : "Box Plot"}
        size="small"
      >
        {isLoading ? (
          <div style={{ textAlign: "center", padding: 48 }}>
            <Spin size="large" />
          </div>
        ) : data ? (
          <>
            <div style={{ marginBottom: 8 }}>
              <Text>
                {data.length.toLocaleString()} points,{" "}
                {new Set(data.map((d) => d[plotObsColumn])).size} categories
              </Text>
              <Text style={{ marginLeft: 16 }}>
                Exclude expression:{" "}
              </Text>
              <Select
                size="small"
                allowClear
                placeholder="Select value to exclude"
                value={filterExpression}
                onChange={(val) => setFilterExpression(val ?? null)}
                options={frequentValues}
                style={{ width: 200 }}
              />
            </div>
            <div ref={containerRef}>
              {violinData && (
                <ViolinPlot
                  groups={violinData.groups}
                  violins={violinData.violins}
                  boxplotStats={boxplotData?.stats}
                  showBoxplot
                  containerWidth={containerWidth}
                  height={500}
                  xLabel={plotObsColumn}
                  yLabel={plotGene}
                />
              )}
            </div>
          </>
        ) : (
          <Text type="secondary">
            Select a gene and an obs column to view the box plot.
          </Text>
        )}
      </Card>
    </TabLayout>
  );
}
