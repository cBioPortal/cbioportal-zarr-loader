import { Layout } from 'antd'
import { BrowserRouter, Routes, Route } from 'react-router-dom'
import { ProfilePage } from '@cbioportal-cell-explorer/profiler'
import Home from './pages/Home'
import View from './pages/View'

const { Content } = Layout

const ENABLE_PROFILER = import.meta.env.VITE_ENABLE_PROFILER === 'true'

function App() {
  return (
    <BrowserRouter basename={import.meta.env.BASE_URL}>
      <Routes>
        <Route path="/" element={
          <Layout style={{ minHeight: '100vh', background: '#fff' }}>
            <Content style={{ maxWidth: 960, margin: '0 auto', padding: '32px 24px' }}>
              <Home />
            </Content>
          </Layout>
        } />
        <Route path="/view" element={<View />} />
        {ENABLE_PROFILER && (
          <Route path="/profile" element={
            <Layout style={{ minHeight: '100vh', background: '#fff' }}>
              <Content style={{ maxWidth: 960, margin: '0 auto', padding: '32px 24px', overflow: 'auto' }}>
                <ProfilePage />
              </Content>
            </Layout>
          } />
        )}
      </Routes>
    </BrowserRouter>
  )
}

export default App
