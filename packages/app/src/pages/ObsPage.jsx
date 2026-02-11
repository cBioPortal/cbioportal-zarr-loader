import useAppStore from "../store/useAppStore";
import ColumnExplorerView from "../views/ColumnExplorerView";

export default function ObsPage() {
  const {
    metadata,
    obsColumnsSelected,
    obsColumnsData,
    obsColumnLoading,
    obsColumnTime,
    obsIndex,
    toggleObsColumn,
    clearObsColumns,
  } = useAppStore();

  const { obsColumns } = metadata;

  return (
    <ColumnExplorerView
      columns={obsColumns}
      selectedColumns={obsColumnsSelected}
      columnsData={obsColumnsData}
      index={obsIndex}
      loading={obsColumnLoading}
      time={obsColumnTime}
      onToggleColumn={toggleObsColumn}
      onClearAll={clearObsColumns}
    />
  );
}
