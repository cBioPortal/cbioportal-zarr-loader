import { Layout, Breadcrumb, Button } from 'antd'
import { BrowserRouter, Routes, Route, Link, useLocation } from 'react-router-dom'
import { ProfilePage } from '@cbioportal-zarr-loader/profiler'
import Home from './pages/Home'
import View from './pages/View'

const { Header, Content } = Layout

const breadcrumbNameMap: Record<string, string> = {
  '/view': 'View',
  '/profile': 'Profile',
}

function AppBreadcrumb() {
  const location = useLocation()
  const pathSnippets = location.pathname.split('/').filter(Boolean)

  const items = [
    { title: <Link to="/">Home</Link> },
    ...pathSnippets.map((_, index) => {
      const url = `/${pathSnippets.slice(0, index + 1).join('/')}`
      const isLast = index === pathSnippets.length - 1
      return {
        title: isLast ? breadcrumbNameMap[url] || url : <Link to={url}>{breadcrumbNameMap[url] || url}</Link>,
      }
    }),
  ]

  return <Breadcrumb style={{ marginBottom: 16 }} items={items} />
}

function App() {
  return (
    <BrowserRouter>
      <Layout style={{ height: '100vh' }}>
        <Header style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', background: '#fff', borderBottom: '1px solid #f0f0f0' }}>
          <Link to="/" style={{ fontSize: 18, fontWeight: 600, color: 'inherit', textDecoration: 'none' }}>
            highperformer
          </Link>
          <Link to="/profile">
            <Button type="text">Profile History</Button>
          </Link>
        </Header>
        <Content style={{ display: 'flex', flexDirection: 'column', flex: 1, overflow: 'hidden' }}>
          <div style={{ padding: '16px 24px 0' }}>
            <AppBreadcrumb />
          </div>
          <div style={{ flex: 1, display: 'flex', flexDirection: 'column', overflow: 'hidden' }}>
            <Routes>
              <Route path="/" element={<div style={{ padding: '0 24px 24px' }}><Home /></div>} />
              <Route path="/view" element={<View />} />
              <Route path="/profile" element={<div style={{ padding: '0 24px 24px', overflow: 'auto', flex: 1 }}><ProfilePage /></div>} />
            </Routes>
          </div>
        </Content>
      </Layout>
    </BrowserRouter>
  )
}

export default App
