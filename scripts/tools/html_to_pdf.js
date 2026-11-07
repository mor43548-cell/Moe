// Render report/randomness_report.html to PDF (A4) with the pre-installed Playwright Chromium.
const path = require('path');
const { chromium } = require('playwright');
(async () => {
  const root = path.resolve(__dirname, '..', '..');
  const b = await chromium.launch(); const p = await b.newPage();
  await p.goto('file://' + path.join(root, 'report', 'randomness_report.html'));
  await p.pdf({ path: path.join(root, 'report', 'randomness_report.pdf'), format: 'A4', printBackground: true,
                margin: { top: '12mm', bottom: '12mm', left: '10mm', right: '10mm' } });
  await b.close();
})();
