import { Modal, Tabs } from 'antd'
import type { ReactNode } from 'react'

interface ChartModalProps {
  title: string
  open: boolean
  onClose: () => void
  chart: ReactNode
  table: ReactNode
}

export default function ChartModal({ title, open, onClose, chart, table }: ChartModalProps) {
  return (
    <Modal
      title={title}
      open={open}
      onCancel={onClose}
      footer={null}
      width={720}
      styles={{ body: { maxHeight: '70vh', overflow: 'auto' } }}
    >
      <Tabs
        size="small"
        items={[
          { key: 'chart', label: 'Chart', children: chart },
          { key: 'table', label: 'Table', children: table },
        ]}
      />
    </Modal>
  )
}
