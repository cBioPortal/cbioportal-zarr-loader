import useAppStore from "../../store/useAppStore";
import SearchableList from "../ui/SearchableList";

export default function TooltipColumnList({ height = 300, width = 220, style = {} }) {
  const { metadata, tooltipColumns, tooltipColumnLoading, toggleTooltipColumn, clearTooltipColumns } = useAppStore();
  const { obsColumns } = metadata || {};

  return (
    <SearchableList
      title="Tooltip"
      items={obsColumns}
      selected={tooltipColumns}
      onSelect={toggleTooltipColumn}
      onClear={clearTooltipColumns}
      loading={tooltipColumnLoading}
      multiSelect
      placeholder="Search columns..."
      height={height}
      width={width}
      style={style}
    />
  );
}
