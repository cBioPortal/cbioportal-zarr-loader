import { useState, useEffect } from 'react'
import { Link } from 'react-router-dom'
import { Input, Button, List, Typography, Space } from 'antd'
import { DeleteOutlined } from '@ant-design/icons'
import { loadDatasets, saveDatasets } from '../utils/datasets'

function Home() {
  const [url, setUrl] = useState('')
  const [datasets, setDatasets] = useState<string[]>(loadDatasets)

  useEffect(() => {
    saveDatasets(datasets)
  }, [datasets])

  const handleAdd = () => {
    const trimmed = url.trim()
    if (!trimmed || datasets.includes(trimmed)) return
    setDatasets((prev) => [trimmed, ...prev])
    setUrl('')
  }

  const handleRemove = (target: string) => {
    setDatasets((prev) => prev.filter((d) => d !== target))
  }

  return (
    <div style={{ maxWidth: 700 }}>
      <h2>Datasets</h2>
      <Space.Compact style={{ width: '100%', marginBottom: 24 }}>
        <Input
          placeholder="Paste a .zarr dataset URL"
          value={url}
          onChange={(e) => setUrl(e.target.value)}
          onPressEnter={handleAdd}
        />
        <Button type="primary" onClick={handleAdd}>
          Add
        </Button>
      </Space.Compact>

      <List
        bordered
        dataSource={datasets}
        locale={{ emptyText: 'No datasets added yet' }}
        renderItem={(item) => (
          <List.Item
            actions={[
              <Button
                key="delete"
                type="text"
                danger
                icon={<DeleteOutlined />}
                onClick={() => handleRemove(item)}
              />,
            ]}
          >
            <Link to={`/view?url=${item}`}>
              <Typography.Text>{item}</Typography.Text>
            </Link>
          </List.Item>
        )}
      />
    </div>
  )
}

export default Home
