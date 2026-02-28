export const MEASURE_PREFIX = "czl:";

export interface ProfileEntry {
  id: number;
  method: string;
  key: string;
  cacheHit: boolean;
  startTime: number;
  duration: number;
}

export interface MeasureDetail {
  key: string;
  cacheHit: boolean;
}

export class ProfileCollector {
  entries: ProfileEntry[] = [];
  enabled: boolean = true;
  #nextId = 1;
  #listeners: Set<() => void> = new Set();
  #observer: PerformanceObserver;

  constructor() {
    this.#observer = new PerformanceObserver((list) => {
      if (!this.enabled) return;
      for (const entry of list.getEntries()) {
        if (!entry.name.startsWith(MEASURE_PREFIX)) continue;
        const detail = (entry as PerformanceMeasure).detail as MeasureDetail;
        if (!detail) continue;
        this.entries.push({
          id: this.#nextId++,
          method: detail.key.split(":")[0],
          key: detail.key,
          cacheHit: detail.cacheHit,
          startTime: entry.startTime,
          duration: entry.duration,
        });
      }
      this.#notify();
    });
    this.#observer.observe({ type: "measure", buffered: false });
  }

  clear(): void {
    this.entries = [];
    this.#nextId = 1;
    this.#notify();
  }

  subscribe(listener: () => void): () => void {
    this.#listeners.add(listener);
    return () => {
      this.#listeners.delete(listener);
    };
  }

  toJSON(): ProfileEntry[] {
    return this.entries;
  }

  dispose(): void {
    this.#observer.disconnect();
    this.#listeners.clear();
  }

  #notify(): void {
    for (const listener of this.#listeners) {
      listener();
    }
  }
}

let measureSeq = 0;

/** Place a start mark and return a finish callback that records the measure. */
export function startMeasure(key: string, cacheHit: boolean): () => void {
  const seq = ++measureSeq;
  const startMark = `${MEASURE_PREFIX}${key}:${seq}:start`;
  const measureName = `${MEASURE_PREFIX}${key}:${seq}`;
  performance.mark(startMark);
  return () => {
    performance.measure(measureName, {
      start: startMark,
      detail: { key, cacheHit } satisfies MeasureDetail,
    });
  };
}
