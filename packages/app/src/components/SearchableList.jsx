import { useState, useMemo } from "react";
import { Card, Input } from "antd";

export default function SearchableList({
  title,
  items,
  selected,
  onSelect,
  loading = null,
  multiSelect = false,
  placeholder = "Search...",
  height = 300,
  width = 220,
  style = {},
}) {
  const [searchText, setSearchText] = useState("");

  const filteredItems = useMemo(() => {
    if (!items) return [];
    if (!searchText) return items;
    const search = searchText.toLowerCase();
    return items.filter((name) => name.toLowerCase().includes(search));
  }, [items, searchText]);

  if (!items || items.length === 0) return null;

  const isSelected = multiSelect
    ? (item) => selected.includes(item)
    : (item) => selected === item;

  const selectedCount = multiSelect ? selected.length : (selected ? 1 : 0);

  return (
    <Card
      size="small"
      title={`${title} (${multiSelect ? selectedCount + "/" : ""}${items.length.toLocaleString()})`}
      style={{ width, height, ...style }}
      styles={{
        body: {
          padding: 0,
          height: "calc(100% - 38px)",
          display: "flex",
          flexDirection: "column",
          overflow: "hidden",
        },
      }}
    >
      <div style={{ padding: 8, borderBottom: "1px solid #f0f0f0" }}>
        <Input.Search
          placeholder={placeholder}
          size="small"
          value={searchText}
          onChange={(e) => setSearchText(e.target.value)}
          allowClear
        />
      </div>
      <div style={{ flex: 1, overflow: "auto" }}>
        {filteredItems.map((item) => (
          <div
            key={item}
            style={{
              padding: "4px 12px",
              cursor: loading === item ? "wait" : "pointer",
              backgroundColor: isSelected(item) ? "#e6f4ff" : undefined,
              fontSize: 12,
              whiteSpace: "nowrap",
              overflow: "hidden",
              textOverflow: "ellipsis",
              opacity: loading === item ? 0.5 : 1,
            }}
            onClick={() => onSelect(item)}
          >
            {item}
          </div>
        ))}
        {filteredItems.length === 0 && searchText && (
          <div style={{ padding: "8px 12px", color: "#999", fontSize: 12 }}>
            No matches found
          </div>
        )}
      </div>
      {searchText && (
        <div
          style={{
            padding: "4px 12px",
            fontSize: 11,
            color: "#999",
            borderTop: "1px solid #f0f0f0",
          }}
        >
          {filteredItems.length.toLocaleString()} matches
        </div>
      )}
    </Card>
  );
}
