import { useEffect, useState } from 'react'
import { Modal, Tabs } from 'antd'
import type { ReactNode } from 'react'

interface ChartModalProps {
  title: string
  open: boolean
  onClose: () => void
  chart: ReactNode
  table: ReactNode
  defaultTab?: 'chart' | 'table'
}

export default function ChartModal({ title, open, onClose, chart, table, defaultTab = 'chart' }: ChartModalProps) {
  const [activeTab, setActiveTab] = useState(defaultTab)

  useEffect(() => {
    if (open) setActiveTab(defaultTab)
  }, [open, defaultTab])

  return (
    <Modal
      title={title}
      open={open}
      onCancel={onClose}
      footer={null}
      width={960}
      styles={{ body: { maxHeight: '75vh', overflow: 'auto' } }}
    >
      <Tabs
        size="small"
        activeKey={activeTab}
        onChange={(key) => setActiveTab(key as 'chart' | 'table')}
        items={[
          { key: 'chart', label: 'Chart', children: chart },
          { key: 'table', label: 'Table', children: table },
        ]}
      />
    </Modal>
  )
}
