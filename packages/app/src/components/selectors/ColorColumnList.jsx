import useAppStore from "../../store/useAppStore";
import SearchableList from "../ui/SearchableList";

export default function ColorColumnList({ height = 300, width = 220, style = {} }) {
  const { metadata, colorColumn, colorLoading, setColorColumn } = useAppStore();
  const { obsColumns } = metadata || {};

  const handleClear = () => setColorColumn(null);

  return (
    <SearchableList
      title="Color By"
      items={obsColumns}
      selected={colorColumn}
      onSelect={setColorColumn}
      onClear={handleClear}
      loading={colorLoading ? colorColumn : null}
      placeholder="Search columns..."
      height={height}
      width={width}
      style={style}
    />
  );
}
