import localFlags from "../../feature-flags.json";

export async function fetchFeatureFlags() {
  const url = import.meta.env.VITE_FEATURE_FLAGS_URL;
  if (!url) return localFlags;

  try {
    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    return await res.json();
  } catch {
    return localFlags;
  }
}
