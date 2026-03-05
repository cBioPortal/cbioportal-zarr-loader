import { useState } from 'react'
import { Segmented, Select, AutoComplete, Alert, Button, Popover, Space } from 'antd'
import { SettingOutlined } from '@ant-design/icons'
import useAppStore from '../store/useAppStore'
import type { ColorMode } from '../store/useAppStore'
import { COLOR_SCALES } from '../utils/colors'
import CategoricalLegend from './CategoricalLegend'
import ContinuousLegend from './ContinuousLegend'

const scaleOptions = Object.keys(COLOR_SCALES).map((name) => ({
  value: name,
  label: name.charAt(0).toUpperCase() + name.slice(1),
}))

function ScaleSettingsButton() {
  const colorScaleName = useAppStore((s) => s.colorScaleName)
  const setColorScaleName = useAppStore((s) => s.setColorScaleName)
  const varColumns = useAppStore((s) => s.varColumns)
  const geneLabelColumn = useAppStore((s) => s.geneLabelColumn)
  const setGeneLabelColumn = useAppStore((s) => s.setGeneLabelColumn)
  const [open, setOpen] = useState(false)

  const labelColumnOptions = [
    { value: '__none__', label: '(None — use index)' },
    ...varColumns.map((col) => ({ value: col, label: col })),
  ]

  return (
    <Popover
      content={
        <div style={{ width: 200 }}>
          <div style={{ fontSize: 12, marginBottom: 4, marginTop: 12 }}>Gene Label Column</div>
          <Select
            value={geneLabelColumn ?? '__none__'}
            onChange={(v) => setGeneLabelColumn(v === '__none__' ? null : v)}
            options={labelColumnOptions}
            style={{ width: '100%' }}
            size="small"
          />
          <div style={{ fontSize: 12, marginBottom: 4 }}>Color Scale</div>
          <Select
            value={colorScaleName}
            onChange={(v) => { setColorScaleName(v); setOpen(false) }}
            options={scaleOptions}
            style={{ width: '100%' }}
            size="small"
          />
        </div>
      }
      trigger="click"
      placement="bottomRight"
      open={open}
      onOpenChange={setOpen}
    >
      <Button type="default" size="small" icon={<SettingOutlined />} />
    </Popover>
  )
}

export default function ColorBySection() {
  const colorMode = useAppStore((s) => s.colorMode)
  const setColorMode = useAppStore((s) => s.setColorMode)
  const obsColumnNames = useAppStore((s) => s.obsColumnNames)
  const varNames = useAppStore((s) => s.varNames)
  const selectedObsColumn = useAppStore((s) => s.selectedObsColumn)
  const selectedGene = useAppStore((s) => s.selectedGene)
  const selectObsColumn = useAppStore((s) => s.selectObsColumn)
  const clearObsColumn = useAppStore((s) => s.clearObsColumn)
  const selectGene = useAppStore((s) => s.selectGene)
  const clearGene = useAppStore((s) => s.clearGene)
  const categoryWarning = useAppStore((s) => s.categoryWarning)
  const geneLabelMap = useAppStore((s) => s.geneLabelMap)

  const [searchText, setSearchText] = useState('')

  const columnOptions = obsColumnNames.map((name) => ({ value: name, label: name }))

  const geneOptions = varNames
    .map((varIndex) => {
      const symbol = geneLabelMap?.get(varIndex)
      return {
        value: varIndex,
        label: symbol ?? varIndex,
      }
    })
    .filter((opt) =>
      !searchText || opt.label.toLowerCase().includes(searchText.toLowerCase()),
    )

  // Show the display label for the selected gene
  const selectedGeneDisplay = selectedGene
    ? (geneLabelMap?.get(selectedGene) ?? selectedGene)
    : undefined

  return (
    <div style={{ padding: '12px 16px', borderBottom: '1px solid #f0f0f0' }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 8 }}>
        <div style={{ fontSize: 12, fontWeight: 600, color: '#666', textTransform: 'uppercase' }}>Color By</div>
        <Segmented
          value={colorMode}
          onChange={(value) => setColorMode(value as ColorMode)}
          size="small"
          options={[{ value: 'category', label: 'Category' }, { value: 'gene', label: 'Gene' }]}
        />
      </div>

      {colorMode === 'category' && (
        <div style={{ marginBottom: 8 }}>
          <Select
            showSearch={{ optionFilterProp: 'label' }}
            allowClear
            placeholder="Select column..."
            value={selectedObsColumn}
            onChange={(v) => v ? selectObsColumn(v) : clearObsColumn()}
            options={columnOptions}
            style={{ width: '100%' }}
            size="small"
          />
          {categoryWarning && (
            <Alert
              title={categoryWarning}
              type="warning"
              showIcon
              style={{ marginTop: 8, fontSize: 12 }}
            />
          )}
        </div>
      )}

      {colorMode === 'gene' && (
        <div style={{ marginBottom: 8 }}>
          <Space.Compact style={{ width: '100%' }}>
            <AutoComplete
              options={geneOptions}
              allowClear
              placeholder="Search gene..."
              value={searchText || selectedGeneDisplay}
              showSearch={{ onSearch: setSearchText }}
              onSelect={(value: string) => { selectGene(value); setSearchText('') }}
              onClear={() => { clearGene(); setSearchText('') }}
              style={{ flex: 1 }}
              size="small"
            />
            <ScaleSettingsButton />
          </Space.Compact>
        </div>
      )}

      {colorMode === 'category' && <CategoricalLegend />}
      {colorMode === 'gene' && <ContinuousLegend />}
    </div>
  )
}
