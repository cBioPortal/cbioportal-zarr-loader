import useAppStore from '../store/useAppStore'

const MAX_DISPLAY = 15

export default function CategoricalLegend() {
  const categoryMap = useAppStore((s) => s.categoryMap)

  if (categoryMap.length === 0) return null

  const displayed = categoryMap.slice(0, MAX_DISPLAY)
  const remaining = categoryMap.length - MAX_DISPLAY

  return (
    <div style={{ marginTop: 8 }}>
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px 12px' }}>
        {displayed.map(({ label, color }) => (
          <div key={label} style={{ display: 'flex', alignItems: 'center', gap: 4, fontSize: 11 }}>
            <span
              style={{
                display: 'inline-block',
                width: 10,
                height: 10,
                backgroundColor: `rgb(${color[0]}, ${color[1]}, ${color[2]})`,
                borderRadius: 2,
                flexShrink: 0,
              }}
            />
            <span style={{ whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis', maxWidth: 100 }}>
              {label}
            </span>
          </div>
        ))}
      </div>
      {remaining > 0 && (
        <div style={{ fontSize: 11, color: '#888', marginTop: 4 }}>
          ...and {remaining} more
        </div>
      )}
    </div>
  )
}
