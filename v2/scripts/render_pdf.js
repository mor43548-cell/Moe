// Render v2/report/Lotto_Randomness_Report_v2.html to PDF with headless Chromium (Playwright).
const { chromium } = require('playwright');
const path = require('path');
(async () => {
  const root = path.resolve(__dirname, '../..');
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto('file://' + path.join(root, 'v2/report/Lotto_Randomness_Report_v2.html'));
  await page.pdf({ path: path.join(root, 'v2/report/Lotto_Randomness_Report_v2.pdf'), format: 'A4', landscape: true,
    margin: { top: '12mm', bottom: '12mm', left: '10mm', right: '10mm' }, printBackground: true,
    displayHeaderFooter: true, headerTemplate: '<span></span>',
    footerTemplate: '<div style="font-size:8px;width:100%;text-align:center">Lotto randomness assessment v2 · page <span class="pageNumber"></span>/<span class="totalPages"></span></div>' });
  await browser.close();
  console.log('PDF written');
})();
