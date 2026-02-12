import { useState, useMemo } from "react";
import { Card, Input, Button } from "antd";

export default function SearchableList({
  title,
  items,
  selected,
  onSelect,
  onClear,
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

  const ITEM_HEIGHT = 25;
  const HEADER_HEIGHT = 38;
  const SEARCH_HEIGHT = 41;
  const SEARCH_THRESHOLD = 10;
  const showSearch = items.length > SEARCH_THRESHOLD;
  const contentHeight = HEADER_HEIGHT + (showSearch ? SEARCH_HEIGHT : 0) + items.length * ITEM_HEIGHT;
  const autoHeight = Math.min(contentHeight, height);

  const isSelected = multiSelect
    ? (item) => selected.includes(item)
    : (item) => selected === item;

  const selectedCount = multiSelect ? selected.length : (selected ? 1 : 0);

  return (
    <Card
      size="small"
      title={`${title} (${multiSelect ? selectedCount + "/" : ""}${items.length.toLocaleString()})`}
      extra={onClear && selectedCount > 0 ? (
        <Button type="link" size="small" onClick={onClear} style={{ padding: 0 }}>
          Clear
        </Button>
      ) : null}
      style={{ width, height: autoHeight, ...style }}
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
      {showSearch && (
        <div style={{ padding: 8, borderBottom: "1px solid #f0f0f0" }}>
          <Input.Search
            placeholder={placeholder}
            size="small"
            value={searchText}
            onChange={(e) => setSearchText(e.target.value)}
            allowClear
          />
        </div>
      )}
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
