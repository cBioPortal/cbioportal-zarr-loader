import { memo, useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { useSearchParams, useNavigate, Link } from 'react-router-dom'
import { Collapse, InputNumber, Layout, Switch, Typography, Select, Spin } from 'antd'
import { BgColorsOutlined, DatabaseOutlined, DotChartOutlined, HolderOutlined, LeftOutlined, RightOutlined, SettingOutlined } from '@ant-design/icons'
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

const LEFT_SIDEBAR_WIDTH = 300
const RIGHT_SIDEBAR_WIDTH = 300
const SIDEBAR_COLLAPSED_WIDTH = 60
// Snap breakpoints for the right sidebar — drag releases snap to nearest
const RIGHT_SNAP_POINTS = [SIDEBAR_COLLAPSED_WIDTH, 200, RIGHT_SIDEBAR_WIDTH, 400, 550]

function snapToNearest(value: number): number {
  let closest = RIGHT_SNAP_POINTS[0]
  let minDist = Math.abs(value - closest)
  for (let i = 1; i < RIGHT_SNAP_POINTS.length; i++) {
    const dist = Math.abs(value - RIGHT_SNAP_POINTS[i])
    if (dist < minDist) { closest = RIGHT_SNAP_POINTS[i]; minDist = dist }
  }
  return closest
}

function useRightSidebarDrag(onDragEnd: (snapped: number) => void) {
  const currentSnapped = useRef(RIGHT_SIDEBAR_WIDTH)
  const ghostRef = useRef<HTMLDivElement | null>(null)
  const setSnappedRef = useCallback((w: number) => { currentSnapped.current = w }, [])

  const onMouseDown = useCallback((e: React.MouseEvent) => {
    e.preventDefault()
    const startX = e.clientX
    const startWidth = currentSnapped.current
    let latestWidth = startWidth
    let lastSnap = snapToNearest(startWidth)

    // Ghost line follows cursor
    const ghost = document.createElement('div')
    Object.assign(ghost.style, {
      position: 'fixed',
      top: '0',
      width: '2px',
      height: '100vh',
      background: '#1677ff',
      opacity: '0.4',
      zIndex: '9999',
      pointerEvents: 'none',
      left: `${e.clientX}px`,
    })

    // Shaded snap region — anchored to right edge, width = snap target
    const snapRegion = document.createElement('div')
    Object.assign(snapRegion.style, {
      position: 'fixed',
      top: '0',
      right: '0',
      width: `${lastSnap}px`,
      height: '100vh',
      background: '#1677ff',
      opacity: '0.06',
      zIndex: '9998',
      pointerEvents: 'none',
      transition: 'width 150ms ease',
      borderLeft: '2px solid rgba(22, 119, 255, 0.3)',
    })

    document.body.appendChild(snapRegion)
    document.body.appendChild(ghost)
    ghostRef.current = ghost

    const onMouseMove = (ev: MouseEvent) => {
      const delta = startX - ev.clientX
      latestWidth = Math.max(0, startWidth + delta)
      ghost.style.left = `${ev.clientX}px`

      const snap = snapToNearest(latestWidth)
      if (snap !== lastSnap) {
        lastSnap = snap
        snapRegion.style.width = `${snap}px`
      }
    }

    const onMouseUp = () => {
      document.removeEventListener('mousemove', onMouseMove)
      document.removeEventListener('mouseup', onMouseUp)
      document.body.style.cursor = ''
      document.body.style.userSelect = ''
      ghost.remove()
      snapRegion.remove()
      ghostRef.current = null
      onDragEnd(snapToNearest(latestWidth))
    }

    document.body.style.cursor = 'col-resize'
    document.body.style.userSelect = 'none'
    document.addEventListener('mousemove', onMouseMove)
    document.addEventListener('mouseup', onMouseUp)
  }, [onDragEnd])

  return { onMouseDown, setSnappedRef }
}

const WIDGETS = [new StatsWidget({ type: 'deck', framesPerUpdate: 5, placement: 'top-left' })]

// Fallback color when no color buffer is ready yet
const FALLBACK_COLOR: [number, number, number, number] = [100, 150, 255, 77]

function urlLabel(url: string): string {
  return url.replace(/\/+$/, '').split('/').pop() ?? url
}

const sectionStyle = { padding: '12px 16px', borderBottom: '1px solid #f0f0f0' } as const
const labelStyle = { fontSize: 12, fontWeight: 600, color: '#666', textTransform: 'uppercase', marginBottom: 8 } as const

const collapsedIconStyle: React.CSSProperties = {
  fontSize: 16,
  color: '#999',
  padding: '12px 0',
  display: 'flex',
  justifyContent: 'center',
}

function CollapsedSidebar({ onExpand }: { onExpand: () => void }) {
  return (
    <div onClick={onExpand} style={{ cursor: 'pointer', display: 'flex', flexDirection: 'column', alignItems: 'center', height: '100%', paddingTop: 12 }}>
      <Link to="/" style={{ textDecoration: 'none', marginBottom: 8 }} onClick={(e) => e.stopPropagation()}>
        <Typography.Text strong style={{ fontSize: 14 }}>hp</Typography.Text>
      </Link>
      <div style={collapsedIconStyle}><DatabaseOutlined /></div>
      <div style={collapsedIconStyle}><DotChartOutlined /></div>
      <div style={collapsedIconStyle}><BgColorsOutlined /></div>
      <div style={{ flex: 1 }} />
      <div style={{ ...collapsedIconStyle, borderTop: '1px solid #f0f0f0', width: '100%', paddingBottom: 12 }}><SettingOutlined /></div>
    </div>
  )
}

function BrandingHeader() {
  return (
    <div style={{ padding: '12px 16px', borderBottom: '1px solid #f0f0f0', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
      <Link to="/" style={{ textDecoration: 'none' }}>
        <Typography.Title level={5} style={{ margin: 0 }}>highperformer</Typography.Title>
      </Link>
    </div>
  )
}

function LeftSidebarContent() {
  const navigate = useNavigate()
  const datasetUrl = useAppStore((s) => s.datasetUrl)
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

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
      <div style={{ flex: 1, overflow: 'auto' }}>
        <BrandingHeader />

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
          <Select
            style={{ width: '100%' }}
            size="small"
            placeholder="Select embedding"
            value={selectedEmbedding}
            onChange={setSelectedEmbedding}
            options={obsmKeys.map((key) => ({ label: key, value: key }))}
            disabled={obsmKeys.length === 0}
          />
        </div>

        <ColorBySection />
      </div>

      <div style={{ borderTop: '1px solid #f0f0f0' }}>
        <RenderingControls />
      </div>
    </div>
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
    <Collapse
      ghost
      size="small"
      items={[{
        key: 'rendering',
        label: <span style={{ fontSize: 12, fontWeight: 600, color: '#666', textTransform: 'uppercase' }}>Rendering</span>,
        children: (
          <>
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
          </>
        ),
      }]}
    />
  )
}

function CanvasLoadingOverlay() {
  const embeddingLoading = useAppStore((s) => s.embeddingLoading)
  const colorBufferLoading = useAppStore((s) => s.colorBufferLoading)

  const message = embeddingLoading
    ? 'Loading embedding…'
    : colorBufferLoading
      ? 'Updating colors…'
      : null

  if (!message) return null

  return (
    <div style={{ position: 'absolute', top: '50%', left: '50%', transform: 'translate(-50%, -50%)', zIndex: 1, textAlign: 'center' }}>
      <Spin />
      <div style={{ marginTop: 8, fontSize: 12, color: '#999' }}>{message}</div>
    </div>
  )
}

function Visualization() {
  const embeddingData = useAppStore((s) => s.embeddingData)
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

const MemoizedVisualization = memo(Visualization)

const edgeTabStyle: React.CSSProperties = {
  position: 'absolute',
  top: '50%',
  transform: 'translateY(-50%)',
  zIndex: 10,
  width: 16,
  height: 48,
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  background: '#fff',
  border: '1px solid #e8e8e8',
  borderRadius: 4,
  cursor: 'pointer',
  fontSize: 10,
  color: '#999',
  boxShadow: '0 1px 3px rgba(0,0,0,0.08)',
}

function EdgeTab({ onClick, onMouseDown, side, icon, cursor }: {
  onClick?: () => void
  onMouseDown?: (e: React.MouseEvent) => void
  side: 'left' | 'right'
  icon: React.ReactNode
  cursor?: string
}) {
  return (
    <div
      role="button"
      tabIndex={0}
      onClick={onClick}
      onMouseDown={onMouseDown}
      onKeyDown={onClick ? (e) => { if (e.key === 'Enter' || e.key === ' ') onClick() } : undefined}
      style={{
        ...edgeTabStyle,
        ...(side === 'left' ? { left: -9 } : { right: -9 }),
        ...(cursor && { cursor }),
      }}
    >
      {icon}
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
      renderLink={(children) => <Link to="/profile">{children}</Link>}
    />
  )
}

function View() {
  const [searchParams] = useSearchParams()
  const datasetUrl = searchParams.get('url')
  const openDataset = useAppStore((s) => s.openDataset)
  const [leftCollapsed, setLeftCollapsed] = useState(false)
  const [rightWidth, setRightWidth] = useState(RIGHT_SIDEBAR_WIDTH)

  const { onMouseDown: onDragStart, setSnappedRef } = useRightSidebarDrag(setRightWidth)

  // Keep the drag hook's ref in sync with the current snapped width
  useEffect(() => { setSnappedRef(rightWidth) }, [rightWidth, setSnappedRef])

  useEffect(() => {
    if (datasetUrl) openDataset(datasetUrl)
  }, [datasetUrl, openDataset])

  return (
    <Layout style={{ height: '100vh', background: '#fff' }}>
      <Layout style={{ flex: 1, overflow: 'hidden', paddingBottom: PROFILE_BAR_HEIGHT }}>
        <Sider
          width={LEFT_SIDEBAR_WIDTH}
          collapsible
          collapsed={leftCollapsed}
          collapsedWidth={SIDEBAR_COLLAPSED_WIDTH}
          trigger={null}
          theme="light"
          style={{ borderRight: leftCollapsed ? undefined : '1px solid #f0f0f0', overflow: 'auto', background: '#fff', transition: 'width 200ms ease' }}
        >
          {leftCollapsed ? <CollapsedSidebar onExpand={() => setLeftCollapsed(false)} /> : <LeftSidebarContent />}
        </Sider>
        <Content style={{ position: 'relative' }}>
          <EdgeTab
            side="left"
            onClick={() => setLeftCollapsed((c) => !c)}
            icon={leftCollapsed ? <RightOutlined /> : <LeftOutlined />}
          />
          <MemoizedVisualization />
          <CanvasLoadingOverlay />
          <EdgeTab
            side="right"
            onMouseDown={onDragStart}
            icon={<HolderOutlined />}
            cursor="col-resize"
          />
        </Content>
        <Sider
          width={rightWidth}
          theme="light"
          style={{ borderLeft: rightWidth <= SIDEBAR_COLLAPSED_WIDTH ? undefined : '1px solid #f0f0f0', background: '#fff', transition: 'width 200ms ease', position: 'relative' }}
        >
          {/* Drag handle — thin strip on the left edge of the right sidebar */}
          <div
            onMouseDown={onDragStart}
            style={{
              position: 'absolute',
              top: 0,
              left: 0,
              width: 4,
              height: '100%',
              cursor: 'col-resize',
              zIndex: 10,
            }}
          />
          {/* Right sidebar — placeholder for future content */}
        </Sider>
      </Layout>
      <ProfileBarWrapper />
    </Layout>
  )
}

export default View
