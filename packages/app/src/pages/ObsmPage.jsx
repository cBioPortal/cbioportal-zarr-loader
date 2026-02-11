import { useEffect } from "react";
import useAppStore from "../store/useAppStore";
import EmbeddingExplorerView from "../views/EmbeddingExplorerView";
import SearchableList from "../components/ui/SearchableList";

export default function ObsmPage() {
  const {
    metadata,
    selectedObsm,
    obsmData,
    obsmLoading,
    obsmTime,
    fetchObsm,
  } = useAppStore();

  const { obsmKeys } = metadata;

  // Auto-fetch UMAP embedding on mount
  useEffect(() => {
    if (!selectedObsm && obsmKeys.length > 0) {
      const umapKey = obsmKeys.find((k) => /umap/i.test(k));
      if (umapKey) {
        fetchObsm(umapKey);
      }
    }
  }, [obsmKeys, selectedObsm, fetchObsm]);

  return (
    <EmbeddingExplorerView
      selectedObsm={selectedObsm}
      obsmData={obsmData}
      obsmLoading={obsmLoading}
      obsmTime={obsmTime}
      onReload={() => fetchObsm(selectedObsm)}
      sidebar={
        <SearchableList
          title="Keys"
          items={obsmKeys}
          selected={selectedObsm}
          onSelect={fetchObsm}
          loading={obsmLoading ? selectedObsm : null}
          placeholder="Search keys..."
          height={200}
        />
      }
    />
  );
}
