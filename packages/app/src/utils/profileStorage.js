const STORAGE_KEY = "czl:profileHistory";
const MAX_SESSIONS = 50;

export function getProfileHistory() {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return [];
    const parsed = JSON.parse(raw);
    if (!Array.isArray(parsed)) return [];
    return parsed;
  } catch {
    return [];
  }
}

export function saveProfileSession(url, nObs, nVar, entries) {
  const history = getProfileHistory();
  history.unshift({
    url,
    timestamp: Date.now(),
    nObs,
    nVar,
    entries,
  });
  if (history.length > MAX_SESSIONS) history.length = MAX_SESSIONS;
  localStorage.setItem(STORAGE_KEY, JSON.stringify(history));
}

export function removeProfileSession(index) {
  const history = getProfileHistory();
  history.splice(index, 1);
  localStorage.setItem(STORAGE_KEY, JSON.stringify(history));
}

export function clearProfileHistory() {
  localStorage.removeItem(STORAGE_KEY);
}
