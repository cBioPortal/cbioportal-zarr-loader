import { Layout } from 'antd'
import { BrowserRouter, Routes, Route } from 'react-router-dom'
import { ProfilePage } from '@cbioportal-zarr-loader/profiler'
import Home from './pages/Home'
import View from './pages/View'

const { Content } = Layout

function App() {
  return (
    <BrowserRouter>
      <Routes>
        <Route path="/" element={
          <Layout style={{ minHeight: '100vh', background: '#fff' }}>
            <Content style={{ maxWidth: 960, margin: '0 auto', padding: '32px 24px' }}>
              <Home />
            </Content>
          </Layout>
        } />
        <Route path="/view" element={<View />} />
        <Route path="/profile" element={
          <Layout style={{ minHeight: '100vh', background: '#fff' }}>
            <Content style={{ maxWidth: 960, margin: '0 auto', padding: '32px 24px', overflow: 'auto' }}>
              <ProfilePage />
            </Content>
          </Layout>
        } />
      </Routes>
    </BrowserRouter>
  )
}

export default App
