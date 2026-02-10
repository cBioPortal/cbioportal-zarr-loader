import { useEffect, useMemo, useState } from "react";
import { Card, Typography, Spin, Segmented, Select } from "antd";
import { Box, Violin } from "@ant-design/charts";
import SearchableList from "./SearchableList";
import TabLayout from "./TabLayout";
import useAppStore from "../store/useAppStore";

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
  const [maxPoints, setMaxPoints] = useState(5000);
  const [filterExpression, setFilterExpression] = useState(null);

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
      category: String(plotObsData[i]),
      expression: Math.round(val * 10000) / 10000,
    }));
    if (filterExpression !== null) {
      return raw.filter((d) => d.expression !== filterExpression);
    }
    return raw;
  }, [plotGeneExpression, plotObsData, filterExpression]);

  if (data) {
    const categories = new Set(data.map((d) => d.category));
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
                Data ready: {data.length.toLocaleString()} points,{" "}
                {new Set(data.map((d) => d.category)).size} categories.
                {" "}Showing: {Math.min(maxPoints, data.length).toLocaleString()}
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
              <div style={{ marginTop: 8 }}>
                <Text>Max points: </Text>
                <Segmented
                  size="small"
                  value={maxPoints}
                  onChange={setMaxPoints}
                  options={(() => {
                    const opts = [5000, 10000];
                    for (let v = 50000; v < data.length; v += 50000) {
                      opts.push(v);
                    }
                    opts.push(data.length);
                    return opts.map((v) => ({
                      value: v,
                      label: v === data.length ? `All (${v.toLocaleString()})` : v.toLocaleString(),
                    }));
                  })()}
                />
              </div>
            </div>
            <Box
              data={data.slice(0, maxPoints)}
              xField="category"
              yField="expression"
              boxType="boxplot"
            />
            <Violin
              data={data.slice(0, maxPoints)}
              xField="category"
              yField="expression"
            />
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
