import useAppStore from "../store/useAppStore";
import SearchableList from "./SearchableList";

export default function GeneList({ height = 300, width = 220, style = {} }) {
  const { metadata, selectedGene, setSelectedGene } = useAppStore();
  const { geneNames } = metadata || {};

  return (
    <SearchableList
      title="Genes"
      items={geneNames}
      selected={selectedGene}
      onSelect={setSelectedGene}
      placeholder="Search genes..."
      height={height}
      width={width}
      style={style}
    />
  );
}
