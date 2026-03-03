import { useState, useEffect, useRef } from "react";

export default function useMediaQuery(query) {
  const mqlRef = useRef(null);
  if (mqlRef.current === null || mqlRef.current.media !== query) {
    mqlRef.current = window.matchMedia(query);
  }

  const [matches, setMatches] = useState(() => mqlRef.current.matches);

  useEffect(() => {
    const mql = mqlRef.current;
    setMatches(mql.matches);

    const handler = (e) => setMatches(e.matches);
    mql.addEventListener("change", handler);
    return () => mql.removeEventListener("change", handler);
  }, [query]);

  return matches;
}
