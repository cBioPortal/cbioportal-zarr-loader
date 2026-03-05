import { useEffect, useMemo, useRef } from 'react'
import { useSearchParams, useNavigate } from 'react-router-dom'
import { InputNumber, Layout, Switch, Typography, Select, Spin } from 'antd'
import { DeckGL } from '@deck.gl/react'
import { OrthographicView } from '@deck.gl/core'
import { ScatterplotLayer } from '@deck.gl/layers'
import { CollisionFilterExtension } from '@deck.gl/extensions'
import { _StatsWidget as StatsWidget } from '@deck.gl/widgets'
import { ProfileBar, PROFILE_BAR_HEIGHT, saveProfileSession } from '@cbioportal-zarr-loader/profiler'
import useAppStore from '../store/useAppStore'
import ColorBySection from '../components/ColorBySection'
import { loadDatasets, saveDatasets } from '../utils/datasets'

const { Sider, Content } = Layout

const WIDGETS = [new StatsWidget({ type: 'deck', framesPerUpdate: 5, placement: 'top-left' })]

// Fallback color when no color buffer is ready yet
const FALLBACK_COLOR: [number, number, number, number] = [100, 150, 255, 77]

function urlLabel(url: string): string {
  return url.replace(/\/+$/, '').split('/').pop() ?? url
}

function Sidebar() {
  const navigate = useNavigate()
  const datasetUrl = useAppStore((s) => s.datasetUrl)
  const loading = useAppStore((s) => s.loading)
  const nObs = useAppStore((s) => s.nObs)
  const nVar = useAppStore((s) => s.nVar)
  const obsmKeys = useAppStore((s) => s.obsmKeys)
  const selectedEmbedding = useAppStore((s) => s.selectedEmbedding)
  const setSelectedEmbedding = useAppStore((s) => s.setSelectedEmbedding)

  // Build dataset options from localStorage, ensuring current URL is included
  const datasetOptions = useMemo(() => {
    const saved = loadDatasets()
    if (datasetUrl && !saved.includes(datasetUrl)) {
      const updated = [datasetUrl, ...saved]
      saveDatasets(updated)
      return updated
    }
    return saved
  }, [datasetUrl])

  const sectionStyle = { padding: '12px 16px', borderBottom: '1px solid #f0f0f0' } as const
  const labelStyle = { fontSize: 12, fontWeight: 600, color: '#666', textTransform: 'uppercase', marginBottom: 8 } as const

  return (
    <Sider width={280} theme="light" style={{ borderRight: '1px solid #f0f0f0', overflow: 'auto' }}>
      <div style={sectionStyle}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline', ...labelStyle }}>
          <span>Dataset</span>
          {nObs != null && (
            <Typography.Text type="secondary" style={{ fontSize: 11, fontWeight: 400, textTransform: 'none' }}>
              {nObs.toLocaleString()} cells &middot; {(nVar ?? 0).toLocaleString()} genes
            </Typography.Text>
          )}
        </div>
        <Select
          style={{ width: '100%' }}
          size="small"
          placeholder="Select dataset"
          value={datasetUrl}
          onChange={(url: string) => navigate(`/view?url=${encodeURIComponent(url)}`)}
          options={datasetOptions.map((url) => ({ label: urlLabel(url), value: url }))}
        />
      </div>

      <div style={sectionStyle}>
        <div style={labelStyle}>Embedding</div>
        {loading ? (
          <Spin size="small" />
        ) : (
          <Select
            style={{ width: '100%' }}
            size="small"
            placeholder="Select embedding"
            value={selectedEmbedding}
            onChange={setSelectedEmbedding}
            options={obsmKeys.map((key) => ({ label: key, value: key }))}
            disabled={obsmKeys.length === 0}
          />
        )}
      </div>

      <ColorBySection />
      <RenderingControls />
    </Sider>
  )
}

function RenderingControls() {
  const pointRadius = useAppStore((s) => s.pointRadius)
  const setPointRadius = useAppStore((s) => s.setPointRadius)
  const opacity = useAppStore((s) => s.opacity)
  const setOpacity = useAppStore((s) => s.setOpacity)
  const antialiasing = useAppStore((s) => s.antialiasing)
  const setAntialiasing = useAppStore((s) => s.setAntialiasing)
  const collisionEnabled = useAppStore((s) => s.collisionEnabled)
  const setCollisionEnabled = useAppStore((s) => s.setCollisionEnabled)
  const collisionRadiusScale = useAppStore((s) => s.collisionRadiusScale)
  const setCollisionRadiusScale = useAppStore((s) => s.setCollisionRadiusScale)

  return (
    <div style={{ padding: '12px 16px' }}>
      <div style={{ fontSize: 12, fontWeight: 600, color: '#666', textTransform: 'uppercase', marginBottom: 8 }}>
        Rendering
      </div>

      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 8 }}>
        <Typography.Text type="secondary" style={{ fontSize: 12 }}>Point radius (px)</Typography.Text>
        <InputNumber min={0.5} max={20} step={0.5} size="small" value={pointRadius} onChange={(v) => v != null && setPointRadius(v)} style={{ width: 70 }} />
      </div>

      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 8 }}>
        <Typography.Text type="secondary" style={{ fontSize: 12 }}>Opacity</Typography.Text>
        <InputNumber min={0.01} max={1} step={0.05} size="small" value={opacity} onChange={(v) => v != null && setOpacity(v)} style={{ width: 70 }} />
      </div>

      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 8 }}>
        <Typography.Text type="secondary" style={{ fontSize: 12 }}>Antialiasing</Typography.Text>
        <Switch size="small" checked={antialiasing} onChange={setAntialiasing} />
      </div>

      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 8 }}>
        <Typography.Text type="secondary" style={{ fontSize: 12 }}>Collision detection</Typography.Text>
        <Switch size="small" checked={collisionEnabled} onChange={setCollisionEnabled} />
      </div>

      {collisionEnabled && (
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
          <Typography.Text type="secondary" style={{ fontSize: 12 }}>Collision scale</Typography.Text>
          <InputNumber min={0.5} max={10} step={0.5} size="small" value={collisionRadiusScale} onChange={(v) => v != null && setCollisionRadiusScale(v)} style={{ width: 70 }} />
        </div>
      )}
    </div>
  )
}

function ColorBufferSpinner() {
  const colorBufferLoading = useAppStore((s) => s.colorBufferLoading)
  if (!colorBufferLoading) return null
  return (
    <div style={{ position: 'absolute', bottom: 16, right: 16, zIndex: 1 }}>
      <Spin size="small" />
    </div>
  )
}

function Visualization() {
  const embeddingData = useAppStore((s) => s.embeddingData)
  const embeddingLoading = useAppStore((s) => s.embeddingLoading)
  const colorBuffer = useAppStore((s) => s.colorBuffer)
  const pointRadius = useAppStore((s) => s.pointRadius)
  const antialiasing = useAppStore((s) => s.antialiasing)
  const collisionEnabled = useAppStore((s) => s.collisionEnabled)
  const collisionRadiusScale = useAppStore((s) => s.collisionRadiusScale)
  const containerRef = useRef<HTMLDivElement>(null)

  // Derive initial view state from data bounds + container size
  const initialViewState = useMemo(() => {
    if (!embeddingData?.bounds) return { target: [0, 0, 0] as [number, number, number], zoom: 1 }

    const { minX, maxX, minY, maxY } = embeddingData.bounds
    const centerX = (minX + maxX) / 2
    const centerY = (minY + maxY) / 2
    const dataWidth = maxX - minX || 1
    const dataHeight = maxY - minY || 1

    // Use container dimensions if available, otherwise a reasonable default
    const el = containerRef.current
    const viewWidth = el?.clientWidth || 800
    const viewHeight = el?.clientHeight || 600

    const padding = 1.1 // 10% padding
    const zoom = Math.log2(Math.min(
      viewWidth / (dataWidth * padding),
      viewHeight / (dataHeight * padding),
    ))

    const minZoom = zoom - 0.25
    const maxZoom = zoom + 4

    const target: [number, number, number] = [centerX, centerY, 0]
    return { target, zoom, minZoom, maxZoom }
  }, [embeddingData])

  // Memoize layer data object — only recreate when position or color buffer changes
  const layerData = useMemo(() => {
    if (!embeddingData) return null
    const attributes: Record<string, { value: Float32Array | Uint8Array; size: number }> = {
      getPosition: { value: embeddingData.positions, size: 2 },
    }
    if (colorBuffer) {
      attributes.getFillColor = { value: colorBuffer, size: 4 }
    }
    return { length: embeddingData.numPoints, attributes }
  }, [embeddingData, colorBuffer])

  const layers = layerData
    ? [
      new ScatterplotLayer({
        id: 'scatterplot',
        data: layerData,
        dataComparator: (a, b) => a === b,
        // Constant fallback when color buffer hasn't arrived yet
        ...(!colorBuffer && { getFillColor: FALLBACK_COLOR }),
        updateTriggers: {
          getFillColor: [colorBuffer],
        },
        getRadius: pointRadius,
        radiusUnits: 'pixels' as const,
        antialiasing,
        ...(collisionEnabled && {
          extensions: [new CollisionFilterExtension()],
          collisionTestProps: { radiusScale: collisionRadiusScale },
        }),
      }),
    ]
    : []

  // Key forces deck.gl to re-initialize when view state changes (new embedding)
  const deckKey = useMemo(
    () => embeddingData ? `${embeddingData.bounds.minX}-${embeddingData.bounds.maxX}` : 'empty',
    [embeddingData],
  )

  return (
    <div ref={containerRef} style={{ position: 'relative', width: '100%', height: '100%' }}>
      {embeddingLoading && (
        <div style={{ position: 'absolute', top: 16, left: '50%', transform: 'translateX(-50%)', zIndex: 1 }}>
          <Spin />
        </div>
      )}
      <DeckGL
        key={deckKey}
        views={new OrthographicView()}
        initialViewState={initialViewState}
        controller
        layers={layers}
        widgets={WIDGETS}
      />
    </div>
  )
}

function ProfileBarWrapper() {
  const adata = useAppStore((s) => s.adata)
  const datasetUrl = useAppStore((s) => s.datasetUrl)

  return (
    <ProfileBar
      profiler={adata?.profiler}
      onSave={(entries: unknown[]) => {
        if (adata) saveProfileSession(datasetUrl, adata.nObs, adata.nVar, entries)
      }}
    />
  )
}

function View() {
  const [searchParams] = useSearchParams()
  const datasetUrl = searchParams.get('url')
  const openDataset = useAppStore((s) => s.openDataset)

  useEffect(() => {
    if (datasetUrl) openDataset(datasetUrl)
  }, [datasetUrl, openDataset])

  return (
    <>
      <Layout style={{ height: '100%', flex: 1, paddingBottom: PROFILE_BAR_HEIGHT }}>
        <Sidebar />
        <Content style={{ position: 'relative' }}>
          <Visualization />
          <ColorBufferSpinner />
        </Content>
      </Layout>
      <ProfileBarWrapper />
    </>
  )
}

export default View
