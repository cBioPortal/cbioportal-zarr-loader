/// <reference types="vite/client" />

declare module '*?worker' {
  const WorkerConstructor: new () => Worker
  export default WorkerConstructor
}

declare module '@cbioportal-zarr-loader/profiler' {
  import type { FC } from 'react'

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  export const ProfileBar: FC<{
    profiler?: any
    onSave?: (entries: unknown[]) => void
  }>
  export const ProfilePage: FC
  export const PROFILE_BAR_HEIGHT: number
  export function saveProfileSession(
    datasetUrl: string | null,
    nObs: number,
    nVar: number,
    entries: unknown[],
  ): void
}
