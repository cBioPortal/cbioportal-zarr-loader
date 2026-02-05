import {themes as prismThemes} from 'prism-react-renderer';

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'anndata-zarr',
  tagline: 'AnnData ↔ Zarr bridge',
  favicon: 'img/favicon.ico',

  url: 'https://anndata-zarr.example.com',
  baseUrl: '/',

  onBrokenLinks: 'throw',

  markdown: {
    hooks: {
      onBrokenMarkdownLinks: 'warn',
    },
  },

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
          path: './content',
          routeBasePath: '/',
          sidebarPath: './sidebars.js',
        },
        blog: false,
        theme: {
          customCss: './src/css/custom.css',
        },
      }),
    ],
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      navbar: {
        title: 'anndata-zarr',
        items: [
          {
            type: 'docSidebar',
            sidebarId: 'progress',
            position: 'left',
            label: 'Progress',
          },
          {
            type: 'docSidebar',
            sidebarId: 'notes',
            position: 'left',
            label: 'Notes',
          },
          {
            href: 'https://github.com/anndata-zarr/anndata-zarr',
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
