import { useState, useEffect, useMemo } from "react";
import { Typography, Space, Select } from "antd";
import DeckGL from "@deck.gl/react";
import { ScatterplotLayer } from "@deck.gl/layers";
import { OrthographicView } from "@deck.gl/core";

const { Text } = Typography;

const COLORS = [
  [31, 119, 180],   // #1f77b4
  [255, 127, 14],   // #ff7f0e
  [44, 160, 44],    // #2ca02c
  [214, 39, 40],    // #d62728
  [148, 103, 189],  // #9467bd
  [140, 86, 75],    // #8c564b
  [227, 119, 194],  // #e377c2
  [127, 127, 127],  // #7f7f7f
  [188, 189, 34],   // #bcbd22
  [23, 190, 207],   // #17becf
  [174, 199, 232],  // #aec7e8
  [255, 187, 120],  // #ffbb78
  [152, 223, 138],  // #98df8a
  [255, 152, 150],  // #ff9896
  [197, 176, 213],  // #c5b0d5
];

export default function EmbeddingScatterplotDeckGL({
  data,
  shape,
  label,
  obsColumns,
  fetchColorData,
  maxPoints = 50000, // deck.gl can handle more points
}) {
  const [colorColumn, setColorColumn] = useState(null);
  const [colorData, setColorData] = useState(null);
  const [colorLoading, setColorLoading] = useState(false);
  const [hoverInfo, setHoverInfo] = useState(null);

  useEffect(() => {
    if (!colorColumn || !fetchColorData) {
      setColorData(null);
      return;
    }
    async function loadColorData() {
      setColorLoading(true);
      try {
        const values = await fetchColorData(colorColumn);
        setColorData(values);
      } catch (err) {
        console.error(err);
        setColorData(null);
      } finally {
        setColorLoading(false);
      }
    }
    loadColorData();
  }, [colorColumn, fetchColorData]);

  const { points, categoryColorMap, bounds } = useMemo(() => {
    if (!data || !shape) return { points: [], categoryColorMap: {}, bounds: null };

    const cols = shape[1];
    const step = Math.max(1, Math.floor(shape[0] / maxPoints));

    // Build category color map
    const categories = new Map();
    if (colorData) {
      for (let i = 0; i < shape[0]; i += step) {
        const cat = String(colorData[i]);
        if (!categories.has(cat)) {
          categories.set(cat, (categories.size % COLORS.length));
        }
      }
    }

    // Build points array and calculate bounds
    const pts = [];
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;

    for (let i = 0; i < shape[0]; i += step) {
      const x = data[i * cols];
      const y = data[i * cols + 1];

      minX = Math.min(minX, x);
      maxX = Math.max(maxX, x);
      minY = Math.min(minY, y);
      maxY = Math.max(maxY, y);

      const cat = colorData ? String(colorData[i]) : "All";
      pts.push({
        position: [x, y],
        category: cat,
        colorIndex: categories.get(cat) ?? 0,
        index: i,
      });
    }

    // Convert categories map to object for legend
    const colorMap = {};
    categories.forEach((colorIdx, cat) => {
      colorMap[cat] = COLORS[colorIdx];
    });

    return {
      points: pts,
      categoryColorMap: colorMap,
      bounds: { minX, maxX, minY, maxY },
    };
  }, [data, shape, colorData, maxPoints]);

  const initialViewState = useMemo(() => {
    if (!bounds) return { target: [0, 0], zoom: 0 };

    const centerX = (bounds.minX + bounds.maxX) / 2;
    const centerY = (bounds.minY + bounds.maxY) / 2;
    const rangeX = bounds.maxX - bounds.minX;
    const rangeY = bounds.maxY - bounds.minY;
    const maxRange = Math.max(rangeX, rangeY);

    // Calculate zoom to fit the data (rough approximation)
    const zoom = Math.log2(400 / maxRange) - 1;

    return {
      target: [centerX, centerY],
      zoom: Math.max(-5, Math.min(zoom, 10)),
    };
  }, [bounds]);

  const layers = [
    new ScatterplotLayer({
      id: "scatterplot",
      data: points,
      getPosition: (d) => d.position,
      getFillColor: (d) => colorData ? COLORS[d.colorIndex] : [24, 144, 255],
      getRadius: 3,
      radiusUnits: "pixels",
      radiusMinPixels: 1,
      radiusMaxPixels: 10,
      opacity: 0.7,
      pickable: true,
      onHover: (info) => setHoverInfo(info.object ? info : null),
    }),
  ];

  const sortedCategories = useMemo(() => {
    if (!colorData || Object.keys(categoryColorMap).length === 0) return [];
    return Object.entries(categoryColorMap)
      .sort((a, b) => {
        const countA = points.filter(p => p.category === a[0]).length;
        const countB = points.filter(p => p.category === b[0]).length;
        return countB - countA;
      });
  }, [categoryColorMap, colorData, points]);

  return (
    <>
      <Space style={{ marginBottom: 16 }} wrap>
        <Text>Color by:</Text>
        <Select
          allowClear
          placeholder="Select obs column"
          style={{ width: 200 }}
          value={colorColumn}
          onChange={setColorColumn}
          loading={colorLoading}
          options={obsColumns?.map((c) => ({ label: c, value: c })) || []}
        />
        {points.length < shape?.[0] && (
          <>
            <Text type="secondary">|</Text>
            <Text type="secondary">
              Showing {points.length.toLocaleString()} of {shape[0].toLocaleString()} points
            </Text>
          </>
        )}
      </Space>

      <div style={{ display: "flex", gap: 16 }}>
        <div
          style={{
            width: 500,
            height: 500,
            position: "relative",
            border: "1px solid #d9d9d9",
            borderRadius: 4,
          }}
        >
          <DeckGL
            views={new OrthographicView({ id: "ortho" })}
            initialViewState={initialViewState}
            controller={true}
            layers={layers}
          />
          {hoverInfo && (
            <div
              style={{
                position: "absolute",
                left: hoverInfo.x + 10,
                top: hoverInfo.y + 10,
                background: "rgba(0,0,0,0.8)",
                color: "white",
                padding: "4px 8px",
                borderRadius: 4,
                fontSize: 12,
                pointerEvents: "none",
              }}
            >
              <div>x: {hoverInfo.object.position[0].toFixed(4)}</div>
              <div>y: {hoverInfo.object.position[1].toFixed(4)}</div>
              {colorData && <div>{colorColumn}: {hoverInfo.object.category}</div>}
            </div>
          )}

          {/* Axis labels */}
          <div
            style={{
              position: "absolute",
              bottom: 4,
              left: "50%",
              transform: "translateX(-50%)",
              fontSize: 12,
              color: "#666",
            }}
          >
            {label}_1
          </div>
          <div
            style={{
              position: "absolute",
              left: 4,
              top: "50%",
              transform: "translateY(-50%) rotate(-90deg)",
              fontSize: 12,
              color: "#666",
            }}
          >
            {label}_2
          </div>
        </div>

        {/* Legend */}
        {colorColumn && sortedCategories.length > 1 && (
          <div style={{ maxHeight: 500, overflow: "auto", fontSize: 12 }}>
            {sortedCategories.map(([cat, color]) => (
              <div key={cat} style={{ display: "flex", alignItems: "center", gap: 4, marginBottom: 2 }}>
                <div
                  style={{
                    width: 12,
                    height: 12,
                    borderRadius: 2,
                    backgroundColor: `rgb(${color.join(",")})`,
                  }}
                />
                <span>{cat}</span>
              </div>
            ))}
          </div>
        )}
      </div>
    </>
  );
}
