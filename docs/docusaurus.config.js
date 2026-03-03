import { themes as prismThemes } from 'prism-react-renderer';
import path from 'path';

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'cbioportal-zarr-loader',
  tagline: 'AnnData ↔ Zarr bridge',
  favicon: 'img/favicon.ico',

  url: 'https://anndata-zarr.example.com',
  baseUrl: process.env.DOCS_BASE_URL || '/docs/',
  trailingSlash: false,

  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          path: path.resolve(__dirname, './content'),
          routeBasePath: '/',
          sidebarPath: path.resolve(__dirname, './sidebars.js'),
        },
        blog: false,
        theme: {
          customCss: path.resolve(__dirname, './src/css/custom.css'),
        },
      }),
    ],
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      navbar: {
        title: 'cbioportal-zarr-loader',
        items: [
          {
            type: 'docSidebar',
            sidebarId: 'api',
            position: 'left',
            label: 'API',
          },
          {
            type: 'docSidebar',
            sidebarId: 'guides',
            position: 'left',
            label: 'Guides',
          },
          {
            type: 'docSidebar',
            sidebarId: 'progress',
            position: 'left',
            label: 'Progress',
          },
          {
            href: 'https://github.com/cBioPortal/cbioportal-zarr-loader',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      prism: {
        theme: prismThemes.github,
        darkTheme: prismThemes.dracula,
      },
    }),
};

export default config;
