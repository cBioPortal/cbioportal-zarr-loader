import { Select, Tag, Typography } from 'antd'

const MAX_SELECTIONS = 5
const TAG_COLORS = ['blue', 'green', 'purple', 'orange', 'cyan'] as const

interface VariablePickerProps {
  label: string
  options: string[]
  selected: string[]
  onAdd: (name: string) => void
  onRemove: (name: string) => void
  labelMap?: Map<string, string> | null
  loading?: Set<string>
}

export default function VariablePicker({ label, options, selected, onAdd, onRemove, labelMap, loading }: VariablePickerProps) {
  const atLimit = selected.length >= MAX_SELECTIONS

  return (
    <div style={{ marginBottom: 12 }}>
      <Typography.Text type="secondary" style={{ fontSize: 11, fontWeight: 600, textTransform: 'uppercase', display: 'block', marginBottom: 4 }}>
        {label} {atLimit && <span style={{ fontWeight: 400, textTransform: 'none' }}>(max {MAX_SELECTIONS})</span>}
      </Typography.Text>

      <Select
        mode="multiple"
        style={{ width: '100%' }}
        size="small"
        placeholder={atLimit ? `Max ${MAX_SELECTIONS} reached` : `Search ${label.toLowerCase()}...`}
        showSearch
        allowClear
        maxTagCount="responsive"
        value={selected}
        onChange={(values: string[]) => {
          const added = values.filter((v) => !selected.includes(v))
          const removed = selected.filter((v) => !values.includes(v))
          for (const name of added) onAdd(name)
          for (const name of removed) onRemove(name)
        }}
        onSelect={(value: string) => {
          if (atLimit) return
          if (!selected.includes(value)) onAdd(value)
        }}
        onDeselect={(value: string) => {
          onRemove(value)
        }}
        tagRender={({ label: tagLabel, value, closable, onClose }) => {
          const idx = selected.indexOf(value as string)
          const color = TAG_COLORS[idx % TAG_COLORS.length]
          const isLoading = loading?.has(value as string)
          return (
            <Tag
              color={color}
              closable={closable}
              onClose={onClose}
              style={{ marginRight: 4, fontSize: 11 }}
            >
              {isLoading ? `${tagLabel}...` : tagLabel}
            </Tag>
          )
        }}
        filterOption={(input, option) => {
          const lower = input.toLowerCase()
          return (option?.label as string ?? '').toLowerCase().includes(lower)
            || (option?.value as string ?? '').toLowerCase().includes(lower)
        }}
        options={options.map((name) => {
          const display = labelMap?.get(name)
          return {
            label: display ?? name,
            value: name,
            disabled: atLimit && !selected.includes(name),
          }
        })}
      />
    </div>
  )
}
