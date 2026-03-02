import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import { analyzer } from "vite-bundle-analyzer";

export default defineConfig({
  base: process.env.VITE_BASE_URL || "/",
  resolve: {
    dedupe: ["react", "react-dom", "react-router"],
  },
  build: {
    rollupOptions: {
      output: {
        manualChunks: {
          "vendor-react": ["react", "react-dom", "react-router", "zustand"],
          "vendor-antd": ["antd", "@ant-design/icons", "@ant-design/cssinjs"],
          "vendor-deckgl": [
            "deck.gl",
            "@deck.gl/core",
            "@deck.gl/layers",
            "@deck.gl/react",
            "@luma.gl/core",
            "@luma.gl/engine",
            "@luma.gl/webgl",
            "@luma.gl/shadertools",
            "@luma.gl/constants",
          ],
          "vendor-visx": ["@visx/visx"],
          "vendor-zod": ["zod"],
        },
      },
    },
  },
  plugins: [
    react(),
    process.env.ANALYZE && analyzer({ openAnalyzer: true }),
    {
      name: "docs-bypass",
      configureServer(server) {
        // Serve /docs/* as static files instead of rewriting to index.html
        server.middlewares.use((req, _res, next) => {
          if (req.url?.startsWith("/docs")) {
            req.headers.accept = "";
          }
          next();
        });
      },
    },
  ],
  test: {
    environment: "jsdom",
    setupFiles: "./src/test/setup.js",
  },
});
