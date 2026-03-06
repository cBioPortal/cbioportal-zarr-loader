import { Select, Tag, Typography } from 'antd'

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
  const available = options.filter((o) => !selected.includes(o))

  return (
    <div style={{ marginBottom: 12 }}>
      <Typography.Text type="secondary" style={{ fontSize: 11, fontWeight: 600, textTransform: 'uppercase', display: 'block', marginBottom: 4 }}>
        {label}
      </Typography.Text>

      <Select
        style={{ width: '100%' }}
        size="small"
        placeholder={`Search ${label.toLowerCase()}...`}
        showSearch
        value={undefined}
        onChange={(value: string) => onAdd(value)}
        filterOption={(input, option) =>
          (option?.label as string ?? '').toLowerCase().includes(input.toLowerCase())
        }
        options={available.map((name) => {
          const display = labelMap?.get(name)
          return {
            label: display ? `${display} (${name})` : name,
            value: name,
          }
        })}
      />

      {selected.length > 0 && (
        <div style={{ marginTop: 6, display: 'flex', flexWrap: 'wrap', gap: 4 }}>
          {selected.map((name) => {
            const display = labelMap?.get(name) ?? name
            const isLoading = loading?.has(name)
            return (
              <Tag
                key={name}
                closable
                onClose={() => onRemove(name)}
                style={{ fontSize: 11 }}
              >
                {isLoading ? `${display}...` : display}
              </Tag>
            )
          })}
        </div>
      )}
    </div>
  )
}
