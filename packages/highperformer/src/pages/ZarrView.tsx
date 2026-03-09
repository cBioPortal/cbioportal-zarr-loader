import { useEffect, useMemo, useState } from 'react'
import { useSearchParams, Link } from 'react-router-dom'
import { Layout, Tree, Typography, Spin, Alert, Tag } from 'antd'
import { FolderOutlined, FileOutlined } from '@ant-design/icons'
import { ZarrStore } from '@cbioportal-cell-explorer/zarrstore'

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

      nodeMap.set(key, {
        type: isArray ? 'array' : 'group',
        metadata: meta,
        ...(isArray && {
          dtype: meta.data_type != null ? String(meta.data_type) : undefined,
          shape: meta.shape as number[] | undefined,
          chunks: chunkGrid?.configuration?.chunk_shape,
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

    let title: React.ReactNode = name
    if (isArray) {
      title = (
        <span>
          {name}
          {meta.dtype && <Tag color="blue" style={{ marginLeft: 6, fontSize: 11 }}>{meta.dtype}</Tag>}
          {meta.shape && <Tag color="green" style={{ marginLeft: 2, fontSize: 11 }}>{JSON.stringify(meta.shape)}</Tag>}
          {meta.chunks && <Tag color="orange" style={{ marginLeft: 2, fontSize: 11 }}>chunks: {JSON.stringify(meta.chunks)}</Tag>}
        </span>
      )
    }

    const treeNode: TreeNode = {
      title,
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

  const treeData = useMemo(() => {
    if (!store?.consolidatedMetadata) return []
    const nodeMap = parseConsolidated(store.consolidatedMetadata, store.zarrVersion)
    return buildTree(nodeMap)
  }, [store])

  const metadataMap = useMemo(() => buildMetadataMap(treeData), [treeData])

  const selectedMetadata = selectedKey ? metadataMap.get(selectedKey) ?? null : null

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
      <Content style={{ maxWidth: 1400, margin: '0 auto', padding: '32px 24px', overflow: 'auto' }}>
        <div style={{ marginBottom: 24 }}>
          <Link to="/" style={{ fontSize: 14 }}>&larr; Home</Link>
          <Typography.Title level={4} style={{ margin: '12px 0 4px' }}>
            Zarr Inspector
          </Typography.Title>
          <Typography.Text type="secondary" copyable style={{ fontSize: 13, wordBreak: 'break-all' }}>
            {url}
          </Typography.Text>
          {store && (
            <div style={{ marginTop: 4 }}>
              <Tag>Zarr v{store.zarrVersion}</Tag>
            </div>
          )}
        </div>

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
          <div style={{ display: 'flex', gap: 24, alignItems: 'flex-start' }}>
            <div style={{ flex: '0 0 50%', overflow: 'auto', maxHeight: 'calc(100vh - 200px)' }}>
              <Typography.Title level={5} style={{ marginBottom: 8 }}>
                Structure
              </Typography.Title>
              <Tree
                showIcon
                defaultExpandAll={false}
                treeData={treeData}
                selectedKeys={selectedKey ? [selectedKey] : []}
                onSelect={(keys) => {
                  setSelectedKey(keys.length > 0 ? String(keys[0]) : null)
                }}
              />
            </div>

            <div style={{ flex: 1, minWidth: 0, position: 'sticky', top: 32 }}>
              {selectedMetadata ? (
                <>
                  <Typography.Title level={5} style={{ marginBottom: 8 }}>
                    {selectedKey}
                  </Typography.Title>
                  <pre style={{
                    background: '#fafafa',
                    border: '1px solid #f0f0f0',
                    borderRadius: 6,
                    padding: 16,
                    fontSize: 13,
                    overflow: 'auto',
                    maxHeight: 'calc(100vh - 250px)',
                  }}>
                    {JSON.stringify(selectedMetadata, null, 2)}
                  </pre>
                </>
              ) : (
                <Typography.Text type="secondary">Select a node to view its metadata</Typography.Text>
              )}
            </div>
          </div>
        )}
      </Content>
    </Layout>
  )
}

export default ZarrView
