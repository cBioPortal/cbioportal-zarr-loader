import useAppStore from "../../store/useAppStore";
import SearchableList from "../ui/SearchableList";

export default function GeneList({ height = 300, width = 220, style = {} }) {
  const { metadata, selectedGene, setSelectedGene, clearGeneSelection } = useAppStore();
  const { geneNames } = metadata || {};

  return (
    <SearchableList
      title="Genes"
      items={geneNames}
      selected={selectedGene}
      onSelect={setSelectedGene}
      onClear={clearGeneSelection}
      placeholder="Search genes..."
      height={height}
      width={width}
      style={style}
    />
  );
}
