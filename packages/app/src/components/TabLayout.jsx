import { Row, Col } from "antd";

export default function TabLayout({ sidebar, children }) {
  return (
    <Row gutter={[16, 16]}>
      <Col xs={24} md={6}>
        {sidebar}
      </Col>
      <Col xs={24} md={18}>
        {children}
      </Col>
    </Row>
  );
}
