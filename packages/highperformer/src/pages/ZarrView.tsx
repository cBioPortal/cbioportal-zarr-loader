import { useEffect, useMemo, useState } from 'react'
import { useSearchParams } from 'react-router-dom'
import { Layout, Tree, Typography, Spin, Alert, Tag, Card } from 'antd'
import { FolderOutlined, FileOutlined, InfoCircleOutlined, CopyOutlined } from '@ant-design/icons'
import { ZarrStore } from '@cbioportal-cell-explorer/zarrstore'
import ChunkShapeViz from '../components/ChunkShapeViz'
import { fetchShardIndex, type ShardIndex } from '../utils/shardIndex'

const { Content } = Layout

interface TreeNode {
  title: React.ReactNode
  key: string
  icon?: React.ReactNode
  children?: TreeNode[]
  isLeaf?: boolean
  /** Raw metadata for this node (e.g. .zarray or .zattrs contents) */
  metadata?: Record<string, unknown>
}

interface NodeMeta {
  type: 'group' | 'array'
  metadata: Record<string, unknown>
  // Array-specific fields extracted for display
  dtype?: string
  shape?: number[]
  chunks?: number[]
  innerChunks?: number[]
}

/**
 * Parse consolidated metadata into a flat map of path → NodeMeta.
 * Handles both v2 and v3 formats.
 *
 * v2 keys: "obs/.zgroup", "obs/.zattrs", "obs/leiden/.zarray"
 * v3 keys: "obs", "obs/leiden" with value { node_type, shape, data_type, ... }
 */
function parseConsolidated(
  consolidatedMetadata: Record<string, unknown>,
  zarrVersion: 2 | 3,
): Map<string, NodeMeta> {
  const nodeMap = new Map<string, NodeMeta>()

  if (zarrVersion === 3) {
    for (const [key, value] of Object.entries(consolidatedMetadata)) {
      const meta = value as Record<string, unknown>
      const nodeType = meta.node_type as string | undefined
      const isArray = nodeType === 'array'
      const chunkGrid = meta.chunk_grid as { configuration?: { chunk_shape?: number[] } } | undefined
      const codecs = meta.codecs as { name: string; configuration?: { chunk_shape?: number[] } }[] | undefined
      const shardingCodec = codecs?.find(c => c.name === 'sharding_indexed')
      const innerChunks = shardingCodec?.configuration?.chunk_shape

      nodeMap.set(key, {
        type: isArray ? 'array' : 'group',
        metadata: meta,
        ...(isArray && {
          dtype: meta.data_type != null ? String(meta.data_type) : undefined,
          shape: meta.shape as number[] | undefined,
          chunks: chunkGrid?.configuration?.chunk_shape,
          ...(innerChunks && { innerChunks }),
        }),
      })
    }
  } else {
    // v2: collect .zarray, .zattrs, .zgroup per path
    const v2Map = new Map<string, { zarray?: Record<string, unknown>; zattrs?: Record<string, unknown>; zgroup?: Record<string, unknown> }>()

    for (const [key, value] of Object.entries(consolidatedMetadata)) {
      let nodePath: string
      let metaType: 'zarray' | 'zattrs' | 'zgroup'

      if (key.endsWith('/.zarray')) {
        nodePath = key.slice(0, -'/.zarray'.length)
        metaType = 'zarray'
      } else if (key.endsWith('/.zattrs')) {
        nodePath = key.slice(0, -'/.zattrs'.length)
        metaType = 'zattrs'
      } else if (key.endsWith('/.zgroup')) {
        nodePath = key.slice(0, -'/.zgroup'.length)
        metaType = 'zgroup'
      } else {
        continue
      }

      if (!v2Map.has(nodePath)) v2Map.set(nodePath, {})
      v2Map.get(nodePath)![metaType] = value as Record<string, unknown>
    }

    for (const [path, meta] of v2Map) {
      const isArray = !!meta.zarray
      const combined: Record<string, unknown> = {}
      if (meta.zgroup) combined['.zgroup'] = meta.zgroup
      if (meta.zattrs) combined['.zattrs'] = meta.zattrs
      if (meta.zarray) combined['.zarray'] = meta.zarray

      const zarray = meta.zarray as { dtype?: string; shape?: number[]; chunks?: number[] } | undefined

      nodeMap.set(path, {
        type: isArray ? 'array' : 'group',
        metadata: combined,
        ...(isArray && zarray && {
          dtype: zarray.dtype != null ? String(zarray.dtype) : undefined,
          shape: zarray.shape,
          chunks: zarray.chunks,
        }),
      })
    }
  }

  return nodeMap
}

/** Build an antd Tree hierarchy from a parsed node map. */
function buildTree(nodeMap: Map<string, NodeMeta>): TreeNode[] {
  const rootChildren: TreeNode[] = []
  const treeNodeLookup = new Map<string, TreeNode>()

  const sortedPaths = [...nodeMap.keys()].sort()

  for (const path of sortedPaths) {
    const meta = nodeMap.get(path)!
    const isArray = meta.type === 'array'
    const parts = path.split('/')
    const name = parts[parts.length - 1]

    const treeNode: TreeNode = {
      title: name,
      key: path,
      icon: isArray ? <FileOutlined /> : <FolderOutlined />,
      children: [],
      isLeaf: isArray,
      metadata: meta.metadata,
    }

    treeNodeLookup.set(path, treeNode)

    // Find parent
    const parentPath = parts.slice(0, -1).join('/')
    if (parentPath && treeNodeLookup.has(parentPath)) {
      treeNodeLookup.get(parentPath)!.children!.push(treeNode)
    } else {
      rootChildren.push(treeNode)
    }
  }

  // Clean up empty children arrays (for leaf display)
  for (const node of treeNodeLookup.values()) {
    if (node.children && node.children.length === 0) {
      delete node.children
    }
  }

  return rootChildren
}

/** Flatten tree to build a key→metadata lookup for the detail panel. */
function buildMetadataMap(nodes: TreeNode[]): Map<string, Record<string, unknown>> {
  const map = new Map<string, Record<string, unknown>>()
  function walk(list: TreeNode[]) {
    for (const node of list) {
      if (node.metadata) {
        map.set(node.key, node.metadata)
      }
      if (node.children) {
        walk(node.children)
      }
    }
  }
  walk(nodes)
  return map
}

function ZarrView() {
  const [searchParams] = useSearchParams()
  const url = searchParams.get('url')

  const [store, setStore] = useState<ZarrStore | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [selectedKey, setSelectedKey] = useState<string | null>(null)
  const [metadataKey, setMetadataKey] = useState<string | null>(null)

  useEffect(() => {
    if (!url) return

    let cancelled = false
    setLoading(true)
    setError(null)
    setStore(null)
    setSelectedKey(null)

    ZarrStore.open(url)
      .then((s) => {
        if (!cancelled) {
          setStore(s)
          setLoading(false)
        }
      })
      .catch((err) => {
        if (!cancelled) {
          setError(err instanceof Error ? err.message : String(err))
          setLoading(false)
        }
      })

    return () => { cancelled = true }
  }, [url])

  const { treeData, nodeMetaMap } = useMemo(() => {
    if (!store?.consolidatedMetadata) return { treeData: [] as TreeNode[], nodeMetaMap: new Map<string, NodeMeta>() }
    const nodeMap = parseConsolidated(store.consolidatedMetadata, store.zarrVersion)
    return { treeData: buildTree(nodeMap), nodeMetaMap: nodeMap }
  }, [store])

  const metadataMap = useMemo(() => buildMetadataMap(treeData), [treeData])

  const selectedMetadata = selectedKey ? metadataMap.get(selectedKey) ?? null : null
  const selectedNodeMeta = selectedKey ? nodeMetaMap.get(selectedKey) ?? null : null

  // Shard selection and index fetching
  const [selectedShard, setSelectedShard] = useState<number | null>(null)
  const [shardIndex, setShardIndex] = useState<ShardIndex | null>(null)
  const [shardIndexLoading, setShardIndexLoading] = useState(false)

  // Shard scan results keyed by array path — persists across node switches
  interface ShardScanEntry {
    indexes: Map<number, ShardIndex>
    failed: number
    loading: boolean
    abort: AbortController | null
  }
  const [shardScans, setShardScans] = useState<Map<string, ShardScanEntry>>(new Map())

  // Derive current node's scan state
  const currentScan = selectedKey ? shardScans.get(selectedKey) : undefined
  const allShardIndexes = currentScan?.indexes ?? new Map<number, ShardIndex>()
  const heatmapFailed = currentScan?.failed ?? 0
  const heatmapLoading = currentScan?.loading ?? false

  // Helper to update a single scan entry
  const updateScan = (key: string, update: Partial<ShardScanEntry>) => {
    setShardScans((prev) => {
      const existing = prev.get(key) ?? { indexes: new Map(), failed: 0, loading: false, abort: null }
      const next = new Map(prev)
      next.set(key, { ...existing, ...update })
      return next
    })
  }

  // Reset shard selection when array changes (scan data is preserved in map)
  useEffect(() => {
    if (selectedNodeMeta?.innerChunks) {
      setSelectedShard(0)
    } else {
      setSelectedShard(null)
    }
    setShardIndex(null)
  }, [selectedKey, selectedNodeMeta])

  // Fetch shard index for the selected shard only
  useEffect(() => {
    setShardIndex(null)
    if (!url || !selectedNodeMeta?.innerChunks || !selectedNodeMeta.chunks || !selectedKey || selectedShard == null) return

    // If we already have it from a previous scan, reuse it
    const cached = allShardIndexes.get(selectedShard)
    if (cached) {
      setShardIndex(cached)
      return
    }

    const chunks = selectedNodeMeta.chunks
    const innerChunks = selectedNodeMeta.innerChunks
    const ndim = chunks.length
    const shardsAlongCols = ndim > 1 ? Math.ceil(selectedNodeMeta.shape![1] / chunks[1]) : 1
    const row = Math.floor(selectedShard / shardsAlongCols)
    const col = selectedShard % shardsAlongCols
    const coords = ndim > 1 ? [row, col] : [row]

    let cancelled = false
    setShardIndexLoading(true)

    fetchShardIndex(url, selectedKey, chunks, innerChunks, coords)
      .then((idx) => {
        if (!cancelled) {
          setShardIndex(idx)
          setShardIndexLoading(false)
          // Add to scan map so it colors progressively
          updateScan(selectedKey!, {
            indexes: new Map(allShardIndexes).set(selectedShard!, idx),
          })
        }
      })
      .catch(() => {
        if (!cancelled) setShardIndexLoading(false)
      })

    return () => { cancelled = true }
  }, [url, selectedKey, selectedNodeMeta, selectedShard, allShardIndexes])

  // Scan all shard indexes (triggered by button)
  const scanAllShards = () => {
    if (!url || !selectedNodeMeta?.innerChunks || !selectedNodeMeta.chunks || !selectedNodeMeta.shape || !selectedKey) return

    const scanKey = selectedKey
    const chunks = selectedNodeMeta.chunks
    const innerChunks = selectedNodeMeta.innerChunks
    const shape = selectedNodeMeta.shape
    const ndim = shape.length > 1 ? 2 : 1
    const shardsAlongRows = Math.ceil(shape[0] / chunks[0])
    const shardsAlongCols = ndim > 1 ? Math.ceil(shape[1] / chunks[1]) : 1
    const totalShards = shardsAlongRows * shardsAlongCols

    const controller = new AbortController()
    updateScan(scanKey, { loading: true, abort: controller, failed: 0 })

    const CONCURRENCY = 4
    const BATCH_SIZE = totalShards <= 50 ? totalShards : Math.min(200, Math.max(50, Math.ceil(totalShards / 20)))
    const results = new Map<number, ShardIndex>()
    let failedCount = 0

    async function run() {
      for (let batchStart = 0; batchStart < totalShards; batchStart += BATCH_SIZE) {
        if (controller.signal.aborted) return

        const batchEnd = Math.min(batchStart + BATCH_SIZE, totalShards)
        let nextIdx = batchStart

        async function worker() {
          while (nextIdx < batchEnd) {
            if (controller.signal.aborted) return
            const idx = nextIdx++
            const row = Math.floor(idx / shardsAlongCols)
            const col = idx % shardsAlongCols
            const coords = ndim > 1 ? [row, col] : [row]
            try {
              const si = await fetchShardIndex(url!, scanKey, chunks, innerChunks, coords)
              if (!controller.signal.aborted) {
                results.set(idx, si)
              }
            } catch {
              failedCount++
            }
          }
        }

        const batchWorkers = Math.min(CONCURRENCY, batchEnd - batchStart)
        await Promise.all(Array.from({ length: batchWorkers }, () => worker()))

        if (!controller.signal.aborted) {
          updateScan(scanKey, { indexes: new Map(results), failed: failedCount })
          await new Promise((r) => requestAnimationFrame(r))
        }
      }

      if (!controller.signal.aborted) {
        updateScan(scanKey, { indexes: new Map(results), failed: failedCount, loading: false })
      }
    }

    run()
  }

  // Compute per-shard total compressed bytes for the tiling heatmap
  const shardBytesMap = useMemo(() => {
    if (allShardIndexes.size === 0) return undefined
    const map = new Map<number, number>()
    for (const [idx, si] of allShardIndexes) {
      const total = si.entries.reduce((sum, e) => sum + (e.nbytes > 0n ? Number(e.nbytes) : 0), 0)
      map.set(idx, total)
    }
    return map
  }, [allShardIndexes])

  if (!url) {
    return (
      <Layout style={{ minHeight: '100vh', background: '#fff' }}>
        <Content style={{ maxWidth: 960, margin: '0 auto', padding: '32px 24px' }}>
          <Alert type="warning" title="No URL provided. Add a ?url= query parameter." />
        </Content>
      </Layout>
    )
  }

  return (
    <Layout style={{ minHeight: '100vh', background: '#fff' }}>
      <Content style={{ padding: '32px 24px', overflow: 'auto' }}>
        <Card size="small" style={{ marginBottom: 24 }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: 8, flexWrap: 'wrap' }}>
            <Typography.Text strong style={{ fontSize: 14 }}>Zarr Inspector</Typography.Text>
            {store && <Tag style={{ margin: 0 }}>Zarr v{store.zarrVersion}</Tag>}
            {store && (
              <Tag color={store.consolidatedMetadata ? 'green' : 'orange'} style={{ margin: 0 }}>
                {store.consolidatedMetadata ? 'consolidated' : 'no consolidated metadata'}
              </Tag>
            )}
            {typeof store?.attrs?.['encoding-type'] === 'string' && (
              <Tag color="blue" style={{ margin: 0 }}>
                {String(store.attrs['encoding-type'])}
                {store.attrs['encoding-version'] ? ` v${String(store.attrs['encoding-version'])}` : ''}
              </Tag>
            )}
          </div>
          <Typography.Text
            type="secondary"
            copyable
            style={{ fontSize: 13, display: 'block', marginTop: 4, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap', maxWidth: '100%' }}
          >
            {url}
          </Typography.Text>
        </Card>

        {loading && (
          <div style={{ textAlign: 'center', padding: '48px 0' }}>
            <Spin size="large" />
            <div style={{ marginTop: 12, color: '#999' }}>Loading Zarr store...</div>
          </div>
        )}

        {error && (
          <Alert type="error" title="Failed to open Zarr store" description={error} showIcon />
        )}

        {store && !store.consolidatedMetadata && (
          <>
            <Alert type="info" title="No consolidated metadata found. Tree view requires consolidated metadata (.zmetadata or zarr.json)." showIcon style={{ marginBottom: 16 }} />
            {Object.keys(store.attrs).length > 0 && (
              <div>
                <Typography.Title level={5}>Root Attributes</Typography.Title>
                <pre style={{
                  background: '#fafafa',
                  border: '1px solid #f0f0f0',
                  borderRadius: 6,
                  padding: 16,
                  fontSize: 13,
                  overflow: 'auto',
                  maxHeight: 400,
                }}>
                  {JSON.stringify(store.attrs, null, 2)}
                </pre>
              </div>
            )}
          </>
        )}

        {treeData.length > 0 && (
          <Layout style={{ background: '#fff' }}>
            <Layout.Sider
              width={360}
              style={{
                background: '#fff',
                borderRight: '1px solid #f0f0f0',
                overflow: 'auto',
                height: 'calc(100vh - 200px)',
                position: 'sticky',
                top: 32,
              }}
            >
              <div style={{ padding: '0 12px 12px 0' }}>
                <Typography.Title level={5} style={{ marginBottom: 8 }}>
                  Structure
                </Typography.Title>
                <Tree
                  showIcon
                  defaultExpandAll={false}
                  treeData={treeData}
                  selectedKeys={selectedKey ? [selectedKey] : []}
                  style={{ whiteSpace: 'nowrap' }}
                  onSelect={(keys) => {
                    setSelectedKey(keys.length > 0 ? String(keys[0]) : null)
                  }}
                  titleRender={(node) => {
                    const key = node.key as string
                    const meta = nodeMetaMap.get(key)
                    const isArray = meta?.type === 'array'
                    return (
                      <span style={{ display: 'inline-flex', alignItems: 'center', gap: 4 }}>
                        <span>{String(node.title)}</span>
                        {isArray && meta?.dtype && <Tag color="blue" style={{ fontSize: 11, margin: 0 }}>{meta.dtype}</Tag>}
                        {isArray && meta?.shape && <Tag color="green" style={{ fontSize: 11, margin: 0 }}>{JSON.stringify(meta.shape)}</Tag>}
                        {isArray && meta?.chunks && (
                          <Tag color="orange" style={{ fontSize: 11, margin: 0 }}>
                            {meta.innerChunks ? 'shard' : 'chunks'}: {JSON.stringify(meta.chunks)}
                          </Tag>
                        )}
                        {isArray && meta?.innerChunks && (
                          <Tag color="volcano" style={{ fontSize: 11, margin: 0 }}>
                            inner chunks: {JSON.stringify(meta.innerChunks)}
                          </Tag>
                        )}
                        <InfoCircleOutlined
                          style={{
                            fontSize: 12,
                            color: metadataKey === key ? '#1890ff' : '#bbb',
                            marginLeft: 4,
                            cursor: 'pointer',
                          }}
                          onClick={(e) => {
                            e.stopPropagation()
                            setMetadataKey(metadataKey === key ? null : key)
                          }}
                        />
                      </span>
                    )
                  }}
                />
              </div>
            </Layout.Sider>

            <Content style={{ padding: '0 0 0 24px', minWidth: 0 }}>
              {selectedMetadata ? (
                <>
                  <Typography.Title level={5} style={{ marginBottom: 8 }}>
                    {selectedKey}
                  </Typography.Title>
                  {selectedNodeMeta && (!selectedNodeMeta.shape?.length || !selectedNodeMeta.chunks?.length) && (
                    <Typography.Text type="secondary">
                      {selectedNodeMeta.type === 'group'
                        ? 'This is a group node. Select a child array to view chunk visualization.'
                        : 'This array has a scalar shape — no chunk layout to visualize.'}
                    </Typography.Text>
                  )}
                  {selectedNodeMeta && selectedNodeMeta.shape && selectedNodeMeta.shape.length > 0 && selectedNodeMeta.chunks && selectedNodeMeta.chunks.length > 0 && (
                    <>
                      <ChunkShapeViz
                        shape={selectedNodeMeta.shape}
                        chunks={selectedNodeMeta.chunks}
                        innerChunks={selectedNodeMeta.innerChunks}
                        dtype={selectedNodeMeta.dtype}
                        shardIndex={shardIndex}
                        selectedShard={selectedShard}
                        onShardSelect={(idx) => {
                          setSelectedShard(idx)
                        }}
                        shardBytesMap={shardBytesMap}
                        onScanAll={scanAllShards}
                        heatmapLoading={heatmapLoading}
                        heatmapFetched={allShardIndexes.size}
                        heatmapFailed={heatmapFailed}
                        onCancelScan={() => {
                          if (selectedKey) {
                            currentScan?.abort?.abort()
                            updateScan(selectedKey, { loading: false, abort: null })
                          }
                        }}
                      />
                      {shardIndexLoading && !heatmapLoading && (
                        <div style={{ color: '#999', fontSize: 12, marginTop: 4 }}>
                          <Spin size="small" /> Loading shard index...
                        </div>
                      )}
                    </>
                  )}
                </>
              ) : (
                <Typography.Text type="secondary">Select a node to view chunk visualization</Typography.Text>
              )}

              {/* Metadata JSON — toggled separately via info icon */}
              {metadataKey && metadataMap.get(metadataKey) && (
                <>
                  <div style={{ display: 'flex', alignItems: 'center', gap: 8, marginTop: 16, marginBottom: 8 }}>
                    <Typography.Title level={5} style={{ margin: 0 }}>
                      Metadata — {metadataKey}
                    </Typography.Title>
                    <CopyOutlined
                      style={{ fontSize: 14, color: '#888', cursor: 'pointer' }}
                      onClick={() => navigator.clipboard.writeText(JSON.stringify(metadataMap.get(metadataKey!), null, 2))}
                    />
                  </div>
                  <pre style={{
                    background: '#fafafa',
                    border: '1px solid #f0f0f0',
                    borderRadius: 6,
                    padding: 16,
                    fontSize: 13,
                    overflow: 'auto',
                    maxHeight: 400,
                  }}>
                    {JSON.stringify(metadataMap.get(metadataKey), null, 2)}
                  </pre>
                </>
              )}
            </Content>
          </Layout>
        )}
      </Content>
    </Layout>
  )
}

export default ZarrView
