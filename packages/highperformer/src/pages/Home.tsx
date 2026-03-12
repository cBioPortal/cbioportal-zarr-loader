import { useMemo, useState, useEffect } from 'react'
import { Link } from 'react-router-dom'
import { Input, Button, List, Tag, Tooltip, Typography } from 'antd'
import { ApartmentOutlined, CheckCircleOutlined, CopyOutlined, DeleteOutlined, ExclamationCircleOutlined, LinkOutlined, LoadingOutlined } from '@ant-design/icons'
import { loadDatasets, saveDatasets } from '../utils/datasets'
import { probeStore, isLocalUrl } from '../utils/datasetProbe'

const ENABLE_ZARR_VIEW = import.meta.env.VITE_ENABLE_ZARR_VIEW === 'true'

interface ProbeResult {
  status: 'pending' | 'ok' | 'error'
  version?: number
}

const StatusIcon = ({ status }: { status: ProbeResult['status'] }) => {
  if (status === 'pending') return <LoadingOutlined style={{ fontSize: 14, color: '#d9d9d9' }} />
  if (status === 'ok') return <CheckCircleOutlined style={{ fontSize: 14, color: '#52c41a' }} />
  return <ExclamationCircleOutlined style={{ fontSize: 14, color: '#ff4d4f' }} />
}

function probeTooltip(result: ProbeResult): string {
  if (result.status === 'pending') return 'Checking...'
  if (result.status === 'error') return 'Unreachable — check CORS, URL, or permissions'
  return `Accessible (Zarr v${result.version})`
}

interface DatasetListProps {
  title: string
  datasets: string[]
  probeResults: Map<string, ProbeResult>
  onRemove: (url: string) => void
}

function DatasetList({ title, datasets, probeResults, onRemove }: DatasetListProps) {
  if (datasets.length === 0) return null

  return (
    <div style={{ marginBottom: 24 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 8 }}>
        <Typography.Text strong>{title}</Typography.Text>
        <Tooltip title="Copy all URLs">
          <Button
            type="text"
            size="small"
            icon={<CopyOutlined />}
            onClick={() => navigator.clipboard.writeText(datasets.join('\n'))}
          />
        </Tooltip>
      </div>
      <List
        bordered
        dataSource={datasets}
        renderItem={(item) => {
          const result = probeResults.get(item) ?? { status: 'pending' as const }
          return (
            <List.Item
              actions={[
                ...(ENABLE_ZARR_VIEW ? [
                  <Tooltip key="inspect" title="Inspect Zarr structure">
                    <Link to={`/zarr_view?url=${encodeURIComponent(item)}`}>
                      <Button type="text" icon={<ApartmentOutlined />} />
                    </Link>
                  </Tooltip>,
                ] : []),
                <Tooltip key="copy" title="Copy zarr URL">
                  <Button
                    type="text"
                    icon={<LinkOutlined />}
                    onClick={() => navigator.clipboard.writeText(item)}
                  />
                </Tooltip>,
                <Button
                  key="delete"
                  type="text"
                  danger
                  icon={<DeleteOutlined />}
                  onClick={() => onRemove(item)}
                />,
              ]}
            >
              <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
                <Tooltip title={probeTooltip(result)}>
                  <StatusIcon status={result.status} />
                </Tooltip>
                <Link to={`/view?url=${encodeURIComponent(item)}`}>
                  <Typography.Text>{item}</Typography.Text>
                </Link>
                {result.status === 'ok' && result.version && (
                  <Tag color="default" style={{ fontSize: 11, lineHeight: '18px', margin: 0 }}>v{result.version}</Tag>
                )}
                {result.status === 'error' && (
                  <Tag color="error" style={{ fontSize: 11, lineHeight: '18px', margin: 0 }}>unreachable</Tag>
                )}
              </div>
            </List.Item>
          )
        }}
      />
    </div>
  )
}

function Home() {
  const [url, setUrl] = useState('')
  const [datasets, setDatasets] = useState<string[]>(loadDatasets)
  const [probeResults, setProbeResults] = useState<Map<string, ProbeResult>>(new Map())

  useEffect(() => {
    saveDatasets(datasets)
  }, [datasets])

  // Probe all datasets on mount and when list changes
  useEffect(() => {
    const controller = new AbortController()

    for (const ds of datasets) {
      setProbeResults((prev) => {
        if (prev.has(ds)) return prev
        const next = new Map(prev)
        next.set(ds, { status: 'pending' })
        return next
      })

      probeStore(ds, controller.signal)
        .then((result) => {
          if (controller.signal.aborted) return
          setProbeResults((prev) => new Map(prev).set(ds, result.ok
            ? { status: 'ok', version: result.version }
            : { status: 'error' },
          ))
        })
        .catch(() => {
          if (controller.signal.aborted) return
          setProbeResults((prev) => new Map(prev).set(ds, { status: 'error' }))
        })
    }

    // Clean up probes for removed datasets
    setProbeResults((prev) => {
      const dsSet = new Set(datasets)
      let changed = false
      const next = new Map(prev)
      for (const key of next.keys()) {
        if (!dsSet.has(key)) { next.delete(key); changed = true }
      }
      return changed ? next : prev
    })

    return () => controller.abort()
  }, [datasets])

  const localDatasets = useMemo(() => datasets.filter(isLocalUrl), [datasets])
  const remoteDatasets = useMemo(() => datasets.filter((d) => !isLocalUrl(d)), [datasets])

  const handleAdd = () => {
    const urls = url.split(/[\n,]+/).map((u) => u.trim()).filter(Boolean)
    const unique = urls.filter((u) => !datasets.includes(u))
    if (unique.length > 0) {
      setDatasets((prev) => [...unique, ...prev])
    }
    setUrl('')
  }

  const handleRemove = (target: string) => {
    setDatasets((prev) => prev.filter((d) => d !== target))
  }

  return (
    <div style={{ maxWidth: 960 }}>
      <h2>Datasets</h2>
      <div style={{ display: 'flex', gap: 8, marginBottom: 24 }}>
        <Input.TextArea
          placeholder="Paste one or more .zarr URLs (one per line or comma-separated)"
          value={url}
          onChange={(e) => setUrl(e.target.value)}
          onPressEnter={(e) => { if (!e.shiftKey) { e.preventDefault(); handleAdd() } }}
          autoSize={{ minRows: 1, maxRows: 4 }}
          style={{ flex: 1 }}
        />
        <Button type="primary" onClick={handleAdd} style={{ alignSelf: 'flex-end' }}>
          Add
        </Button>
      </div>

      <DatasetList
        title="Remote"
        datasets={remoteDatasets}
        probeResults={probeResults}
        onRemove={handleRemove}
      />
      <DatasetList
        title="Local"
        datasets={localDatasets}
        probeResults={probeResults}
        onRemove={handleRemove}
      />

      {datasets.length === 0 && (
        <Typography.Text type="secondary">No datasets added yet</Typography.Text>
      )}
    </div>
  )
}

export default Home
