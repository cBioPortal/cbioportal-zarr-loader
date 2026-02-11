import useAppStore from "../store/useAppStore";
import ColumnExplorerView from "../views/ColumnExplorerView";

export default function VarPage() {
  const {
    metadata,
    varColumnsSelected,
    varColumnsData,
    varColumnLoading,
    varColumnTime,
    varIndex,
    toggleVarColumn,
    clearVarColumns,
  } = useAppStore();

  const { varColumns } = metadata;

  return (
    <ColumnExplorerView
      columns={varColumns}
      selectedColumns={varColumnsSelected}
      columnsData={varColumnsData}
      index={varIndex}
      loading={varColumnLoading}
      time={varColumnTime}
      onToggleColumn={toggleVarColumn}
      onClearAll={clearVarColumns}
    />
  );
}
