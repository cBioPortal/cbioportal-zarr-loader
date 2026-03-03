import { List } from "react-window";

const ROW_HEIGHT = 28;

function Row({ index, style, items, selected, loading, onSelect }) {
  const item = items[index];
  const isSelected = item === selected;
  const isLoading = item === loading;

  return (
    <div
      style={{
        ...style,
        padding: "0 16px",
        fontSize: 13,
        lineHeight: `${ROW_HEIGHT}px`,
        cursor: isLoading ? "wait" : "pointer",
        backgroundColor: isSelected ? "#e6f4ff" : undefined,
        opacity: isLoading ? 0.5 : 1,
        whiteSpace: "nowrap",
        overflow: "hidden",
        textOverflow: "ellipsis",
        boxSizing: "border-box",
      }}
      onClick={() => onSelect(item)}
      title={item}
    >
      {item}
    </div>
  );
}

export default function VirtualizedList({
  items,
  selected,
  onSelect,
  loading = null,
  height = 300,
}) {
  if (!items || items.length === 0) return null;

  return (
    <List
      rowComponent={Row}
      rowCount={items.length}
      rowHeight={ROW_HEIGHT}
      rowProps={{ items, selected, loading, onSelect }}
      style={{ height }}
    />
  );
}
