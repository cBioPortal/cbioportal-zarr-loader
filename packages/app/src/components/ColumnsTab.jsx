import { useState } from "react";
import { Segmented } from "antd";
import ColumnExplorer from "./ColumnExplorer";

export default function ColumnsTab({
  obsColumns,
  obsColumnsSelected,
  obsColumnsData,
  obsIndex,
  obsColumnLoading,
  obsColumnTime,
  toggleObsColumn,
  clearObsColumns,
  varColumns,
  varColumnsSelected,
  varColumnsData,
  varIndex,
  varColumnLoading,
  varColumnTime,
  toggleVarColumn,
  clearVarColumns,
}) {
  const [activeGroup, setActiveGroup] = useState("obs");

  return (
    <div>
      <Segmented
        value={activeGroup}
        onChange={setActiveGroup}
        options={[
          { label: `obs (${obsColumns.length})`, value: "obs" },
          { label: `var (${varColumns.length})`, value: "var" },
        ]}
        style={{ marginBottom: 12 }}
      />
      {activeGroup === "obs" ? (
        <ColumnExplorer
          columns={obsColumns}
          selectedColumns={obsColumnsSelected}
          columnsData={obsColumnsData}
          index={obsIndex}
          loading={obsColumnLoading}
          time={obsColumnTime}
          onToggleColumn={toggleObsColumn}
          onClearAll={clearObsColumns}
        />
      ) : (
        <ColumnExplorer
          columns={varColumns}
          selectedColumns={varColumnsSelected}
          columnsData={varColumnsData}
          index={varIndex}
          loading={varColumnLoading}
          time={varColumnTime}
          onToggleColumn={toggleVarColumn}
          onClearAll={clearVarColumns}
        />
      )}
    </div>
  );
}
