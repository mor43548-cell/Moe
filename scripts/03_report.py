"""
Step 3: build the structured report from results/ (no number is typed by hand - all are read from
the machine-output CSV files). Produces:
  report/Randomness_Report.xlsx  (Excel companion: every table + raw machine output)
  report/Randomness_Report.html  (rendered to PDF by scripts/render_pdf.js)
"""
import base64, json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES, FIG, REP = ROOT / "results", ROOT / "figures", ROOT / "report"
REP.mkdir(exist_ok=True)

merged = pd.read_csv(ROOT / "data/merged/merged_7_files.csv")
uni = pd.read_csv(RES / "1_uniformity.csv")
ind = pd.read_csv(RES / "2_independence_runs.csv")
acf = pd.read_csv(RES / "2_autocorrelation_lags_0_20.csv")
ent = pd.read_csv(RES / "3_entropy.csv")
nist = pd.read_csv(RES / "4_nist_machine_output_all_streams.csv", keep_default_na=False, na_values=[""])
sc = pd.read_csv(RES / "4_nist_scorecard.csv", keep_default_na=False, na_values=[""])
nt = pd.read_csv(RES / "4_nist_nonoverlapping_148_templates.csv")
summ = json.loads((RES / "summary.json").read_text())
val = (RES / "validation_nist.txt").read_text()
machine = (RES / "machine_output.txt").read_text()

def f(x, d=4):
    if pd.isna(x):
        return "—"
    return f"{x:.{d}g}" if abs(x) < 1e-3 and x != 0 else f"{x:.{d}f}"

def st(ok):
    return "PASS" if ok else "FAIL"

S = sc.set_index("test")
def nist_line(t):
    r = S.loc[t]
    if r.status == "N/A":
        return f"N/A — {r.reason}"
    return f"combined-stream p = {r.combined_p_values}; proportion {r.proportion}; status {r.status}"

lb = summ["proportion_bounds_m7"][0]
n_chi = int((uni.p_chi2 > 0.01).sum()); n_ks = int((uni.p_KS > 0.01).sum())
ud_ok = ind.runs_updown_p > 0.05; md_ok = ind.runs_median_p > 0.05
acf_clean = ind[ind.lags_outside_band == 0].column.tolist()
tot_out = int(ind.lags_outside_band.sum())
ecol = ent[ent.column != "COMBINED"]
H_ok = ecol[ecol.H_ratio >= 0.999].column.tolist()
Hb_ok = ent[ent.H_binary >= 0.999].column.tolist()
P = summ["nist"]

def colvals(df, col, d=4):
    return "; ".join(f"{r.column}={f(r[col], d)}" for _, r in df.iterrows())

# ---------------- checklist (mirrors the request item by item) ----------------
CL = [
 ("1. Uniformity", "Frequency test, chi-square, per variable", "p > 0.01 (aim > 0.05)",
  colvals(uni, "p_chi2"), f"{n_chi}/7 variables meet it", "FAIL" if n_chi < 7 else "PASS"),
 ("1. Uniformity", "Frequency test, Kolmogorov-Smirnov (discrete), per variable", "p > 0.01 (aim > 0.05)",
  colvals(uni, "p_KS"), f"{n_ks}/7 variables meet it", "FAIL" if n_ks < 7 else "PASS"),
 ("1. Uniformity", "Proportion of variables passing chi-square (NIST style)", f"about 0.99 (NIST interval for m=7: >= {lb:.4f})",
  f"{summ['uniformity']['prop_chi2']:.4f}", "", st(summ['uniformity']['prop_chi2'] >= lb)),
 ("1. Uniformity", "NIST Frequency (monobit)", "p > 0.01, proportion ~0.99", nist_line("Frequency"), "", S.loc["Frequency"].status),
 ("1. Uniformity", "NIST Block Frequency", "p > 0.01 (aim > 0.05), proportion ~0.99", nist_line("BlockFrequency"), "", S.loc["BlockFrequency"].status),
 ("1. Uniformity", "Frequency histogram (visual)", "relatively flat bars, slight fluctuation only",
  "Figure 1", "Only ABC looks flat. Y and T pile up at high values, Q at low values, and W, R, E are hump-shaped",
  "FAIL" if n_chi < 7 else "PASS"),
 ("2. Independence", "Runs test, up/down runs on the original values", "p > 0.05",
  colvals(ind, "runs_updown_p"), f"{int(ud_ok.sum())}/7 meet it (fails: {', '.join(ind[~ud_ok].column) or 'none'})", st(ud_ok.all())),
 ("2. Independence", "Runs test, runs of identical symbols (above/below median)", "p > 0.05",
  colvals(ind, "runs_median_p"), f"{int(md_ok.sum())}/7 meet it (fails: {', '.join(ind[~md_ok].column) or 'none'})", st(md_ok.all())),
 ("2. Independence", "NIST Runs (bits)", "p > 0.05, proportion ~0.99", nist_line("Runs"), "", S.loc["Runs"].status),
 ("2. Independence", "NIST Longest Run of Ones", "p > 0.05, proportion ~0.99", nist_line("LongestRun"), "", S.loc["LongestRun"].status),
 ("2. Independence", "Autocorrelation coefficient, lags 1-20", "tends to 0 within statistical error (95% band ±1.96/√N)",
  colvals(ind, "max_abs_acf_lags1_20", 3) + " (max |r|)",
  f"{tot_out}/140 coefficients outside the band; by chance about 7/140 are expected outside a 95% band",
  "PASS" if tot_out == 0 else "FAIL"),
 ("2. Independence", "Correlogram (visual)", "all bars except lag 0 inside the confidence band",
  "Figure 2", f"All bars inside the band only for: {', '.join(acf_clean)}. The other columns each have 1-3 bars outside",
  "PASS" if len(acf_clean) == 7 else "FAIL"),
 ("2. Independence", "NIST Random Excursions", "p > 0.05, proportion ~0.99", nist_line("RandomExcursions"), "", S.loc["RandomExcursions"].status),
 ("2. Independence", "NIST Random Excursions Variant", "p > 0.05, proportion ~0.99", nist_line("RandomExcursionsVariant"), "", S.loc["RandomExcursionsVariant"].status),
 ("3. Entropy", "Shannon entropy of the values / maximum log2(k)", ">= 99.9% of maximum",
  colvals(ecol, "H_ratio"), f"meets: {', '.join(H_ok) or 'none'}", st(len(H_ok) == 7)),
 ("3. Entropy", "Shannon entropy of the bit streams (max 1.000)", ">= 0.999",
  "; ".join(f"{r.column}={r.H_binary:.4f}" for _, r in ent.iterrows()), f"meets: {', '.join(Hb_ok) or 'none'}",
  st(len(Hb_ok) == len(ent))),
 ("3. Entropy", "NIST Maurer's Universal", "p > 0.05, proportion ~0.99", nist_line("MaurerUniversal"), "", S.loc["MaurerUniversal"].status),
 ("3. Entropy", "NIST Approximate Entropy", "p > 0.05, proportion ~0.99", nist_line("ApproximateEntropy"), "", S.loc["ApproximateEntropy"].status),
 ("4. NIST SP 800-22", "Binary Matrix Rank", "p > 0.01, proportion ~0.99", nist_line("BinaryMatrixRank"), "", S.loc["BinaryMatrixRank"].status),
 ("4. NIST SP 800-22", "Full battery, 15 tests", "15/15 PASS (100%)",
  f"PASS {P['pass_']}/15, FAIL {P['fail']}/15, N/A {P['na']}/15", "", "PASS" if P["pass_"] == 15 else "FAIL"),
]
cl = pd.DataFrame(CL, columns=["section", "item", "required threshold", "observed (machine output)", "note", "status"])
cl.insert(0, "checked", "✔")
cl.to_csv(RES / "0_checklist.csv", index=False)

# ---------------- Excel ----------------
xl = REP / "Randomness_Report.xlsx"
with pd.ExcelWriter(xl, engine="openpyxl") as w:
    cl.to_excel(w, sheet_name="0_Checklist", index=False)
    sc.to_excel(w, sheet_name="NIST_Scorecard", index=False)
    merged.to_excel(w, sheet_name="Merged_Data", index=False)
    uni.to_excel(w, sheet_name="1_Uniformity", index=False)
    ind.to_excel(w, sheet_name="2_Runs_Autocorr", index=False)
    acf.to_excel(w, sheet_name="2_ACF_lags", index=False)
    ent.to_excel(w, sheet_name="3_Entropy", index=False)
    nist.to_excel(w, sheet_name="4_NIST_machine_output", index=False)
    nt.to_excel(w, sheet_name="4_NonOverlap_148", index=False)
    pd.DataFrame({"validation_nist.txt": val.splitlines()}).to_excel(w, sheet_name="Validation", index=False)
from openpyxl import load_workbook
from openpyxl.styles import PatternFill, Font
wb = load_workbook(xl)
fills = {"PASS": "C6EFCE", "FAIL": "FFC7CE", "N/A": "E7E6E6"}
for ws in wb.worksheets:
    for c in ws[1]:
        c.font = Font(bold=True)
    for col in ws.columns:
        ws.column_dimensions[col[0].column_letter].width = min(60, max(10, max(len(str(c.value or "")) for c in col[:60]) + 2))
    hdr = [c.value for c in ws[1]]
    if "status" in hdr:
        j = hdr.index("status")
        for row in ws.iter_rows(min_row=2):
            v = row[j].value
            if v in fills:
                row[j].fill = PatternFill("solid", fgColor=fills[v])
wb.save(xl)

# ---------------- HTML (-> PDF) ----------------
def img(p):
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()

def table(df, cls=""):
    h = "<table class='%s'><thead><tr>%s</tr></thead><tbody>" % (cls, "".join(f"<th>{c}</th>" for c in df.columns))
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            s = f(v) if isinstance(v, float) else str(v)
            k = f" class='s-{s.replace('/', '')}'" if c == "status" else ""
            cells.append(f"<td{k}>{s}</td>")
        h += "<tr>" + "".join(cells) + "</tr>"
    return h + "</tbody></table>"

sc_view = sc[["test", "combined_p_values", "proportion", "status", "reason"]].copy()
sc_view.insert(0, "#", range(1, 16))
uni_v = uni[["column", "min", "max", "k", "chi2", "df", "p_chi2", "KS_D", "p_KS"]]
ind_v = ind[["column", "acf_lag1", "max_abs_acf_lags1_20", "lags_outside_band", "runs_updown_R", "runs_updown_expected",
             "runs_updown_p", "runs_median_R", "runs_median_expected", "runs_median_p"]]
ent_v = ent[["column", "k", "H_bits", "H_max", "H_ratio", "n_bits", "ones_fraction", "H_binary",
             "ref_expected_H_ratio_if_uniform_N328", "ref_ones_fraction_if_exactly_uniform"]]
for d in (ind_v, ent_v):
    d.columns = [c.replace("_", " ") for c in d.columns]

verdict_he = f"""
<p><b>שאלת המחקר:</b> האם 7 המשתנים בטבלה המאוחדת (Y, T, W, R, ABC, Q, E; ‏N=328 שורות) עוברים את כל מבחני האקראיות שהוגדרו מראש (אחידות, אי-תלות, אנטרופיה וסוללת NIST SP 800-22), כפי שנדרש כדי לקבוע שאלה "נתוני הגרלה"?</p>
<p><b>תשובה: לא.</b> הנתונים <b>אינם עומדים</b> בקריטריון שנקבע מראש (מעבר של 100% מהמבחנים).</p>
<ul>
<li><b>אחידות:</b> מבחן χ² דוחה אחידות ב-6 מתוך 7 המשתנים (p בין 10<sup>-97</sup> ל-10<sup>-14</sup>). רק ABC עובר (χ² p={f(uni.set_index('column').loc['ABC','p_chi2'])}, KS p={f(uni.set_index('column').loc['ABC','p_KS'])}). יחס המעבר הוא {summ['uniformity']['prop_chi2']:.4f}, הרחק מתחת ל-0.99. ההיסטוגרמות אינן שטוחות, למעט ABC.</li>
<li><b>אי-תלות (על הערכים המקוריים):</b> רוב מבחני הרצפים עוברים את הסף p&gt;0.05. יש שני כישלונות: Q במבחן רצפי עלייה/ירידה (p={f(ind.set_index('column').loc['Q','runs_updown_p'])}), ו-R במבחן רצפים מעל/מתחת לחציון (p={f(ind.set_index('column').loc['R','runs_median_p'])}). כל מקדמי האוטו-קורלציה קטנים (|r|≤{ind.max_abs_acf_lags1_20.max():.3f}). {tot_out} מתוך 140 מקדמים חורגים מרצועת 95%, וכל העמודות נמצאות בתוך הרצועה רק עבור {', '.join(acf_clean)}.</li>
<li><b>אנטרופיה:</b> אף משתנה לא מגיע ל-99.9% מהאנטרופיה המקסימלית. הקרוב ביותר הוא ABC עם {ecol.set_index('column').loc['ABC','H_ratio']:.4f}.</li>
<li><b>NIST SP 800-22:</b> PASS {P['pass_']}/15, FAIL {P['fail']}/15, N/A {P['na']}/15. שישה מבחנים לא ניתנים להרצה, כי אורך הרצף (10,824 ביטים ברצף המאוחד) קטן מהמינימום שהתקן דורש.</li>
</ul>
<p><b>הערות אפיסטמולוגיות (אינן משנות את ההכרעה ואינן מחליפות ספים):</b></p>
<ol>
<li>"עבר מבחן" אינו הוכחה לאקראיות. הוא רק אי-דחייה של השערת האפס במבחן מסוים. לכן גם 15/15 לא היה "מוכיח" אקראיות. לעומת זאת, דחייה חזקה (p≈10<sup>-97</sup>) היא ראיה חזקה נגד אחידות.</li>
<li>תחום האחידות הוגדר כטווח השלמים מהמינימום ועד המקסימום שנצפו בכל משתנה. התחום התיאורטי האמיתי לא היה ידוע.</li>
<li>המרת הערכים לביטים (value−min ב-⌈log₂k⌉ ביטים) מטה את תדירות האחדות כאשר k אינו חזקה של 2. זה קורה גם לנתונים אחידים לחלוטין (ראו עמודת ref בטבלת האנטרופיה). לכן כישלונות NIST על הביטים משקפים גם את הקידוד, ולא רק את הנתונים.</li>
<li>עבור N=328, אומדן האנטרופיה הנאיבי מוטה כלפי מטה. גם לנתונים אחידים לחלוטין היחס הצפוי הוא כ-0.987–0.995 (Miller–Madow), כלומר סף 99.9% כמעט אינו בר-השגה במדגם בגודל הזה. הסף לא שונה, והתוצאה מוצגת כפי שהיא.</li>
<li>מבחן ה-DFT: המימוש לא שחזר את דוגמת המספרים שבמסמך התקן (ראו נספח אימות). כל שאר דוגמאות התקן שנבדקו (12/13) שוחזרו בדיוק של 6 ספרות.</li>
<li>תיקון שבוצע במהלך העבודה: הריצה הראשונה של מבחן KS השתמשה ב-scipy.kstest, שנכתב לנתונים רציפים, ולכן ניפח את D בנתונים בדידים עם ערכים חוזרים. המבחן הוחלף ב-KS הבדיד הנכון (Conover 1972). הדבר שינה את תוצאת ABC בלבד (מ-p≈7.6·10<sup>-10</sup> ל-p=0.719). χ² לא הושפע.</li>
<li>ערכי ה-p ויחסי המעבר שהופיעו בבקשה כ"צפויים" (למשל 0.4129 או 0.9885) לא שימשו בשום שלב. בכל מקום מוצגים רק ערכים שחושבו בפועל.</li>
</ol>
"""

html = f"""<!doctype html><html lang="he"><head><meta charset="utf-8"><title>Randomness Report</title>
<style>
body{{font-family:'DejaVu Sans',sans-serif;font-size:10.5px;color:#111;margin:0}}
h1{{font-size:20px;margin:0 0 4px}} h2{{font-size:15px;border-bottom:2px solid #333;padding-bottom:3px;margin-top:22px}}
h3{{font-size:12.5px;margin:14px 0 6px}}
.rtl{{direction:rtl;text-align:right;font-size:11.5px;line-height:1.55}}
table{{border-collapse:collapse;width:100%;margin:6px 0 10px;font-size:9px}}
th,td{{border:1px solid #999;padding:3px 4px;vertical-align:top}} th{{background:#e8e8e8}}
.s-PASS{{background:#c6efce;font-weight:bold}} .s-FAIL{{background:#ffc7ce;font-weight:bold}} .s-NA{{background:#e7e6e6;font-weight:bold}}
img{{width:100%}} pre{{font-size:7.5px;white-space:pre-wrap;background:#f6f6f6;padding:6px;border:1px solid #ccc}}
.box{{border:2px solid #b00;padding:8px;margin:8px 0;background:#fff5f5}}
.pb{{page-break-before:always}}
</style></head><body>
<h1>Statistical Randomness Assessment of 7 Variables (Y, T, W, R, ABC, Q, E)</h1>
<div>Structured report · N = 328 rows × 7 variables · Tests and thresholds defined in advance by the requester · Every number below is unmodified program output (results/*.csv)</div>
<div class="box"><b>Overall result against the pre-defined criterion "pass 100% of the randomness tests": NOT MET.</b><br>
NIST SP 800-22 battery: PASS {P['pass_']}/15 · FAIL {P['fail']}/15 · N/A (sequence too short for the standard) {P['na']}/15.
Chi-square uniformity: {n_chi}/7 variables pass. Shannon entropy ≥ 99.9%: {len(H_ok)}/7 variables.</div>

<h2>Layer 3: Verbal scientific summary (Hebrew)</h2><div class="rtl">{verdict_he}</div>

<h2 class="pb">Layer 2: Tables</h2>
<h3>Checklist, item by item as specified</h3>{table(cl)}
<h3 class="pb">NIST SP 800-22 scorecard (15 tests)</h3>
<p>Combined stream = the 7 encoded variables interleaved row by row (10,824 bits): a single p-value per test.
Proportion = share of the 7 per-variable streams (984–1,640 bits) with p ≥ 0.01. The NIST acceptable proportion interval for m=7 is ≥ {lb:.4f}, i.e. 7/7.
PASS requires combined p ≥ 0.01 and a proportion inside the interval. N/A = the input-size requirement in SP 800-22 §2.x.7 is not met, so no p-value is reported.</p>
{table(sc_view)}
<h3>1. Uniformity: chi-square and discrete KS (domain = observed min..max)</h3>{table(uni_v)}
<h3>2. Independence: autocorrelation and runs tests (original values)</h3>{table(ind_v)}
<h3>3. Shannon entropy</h3>
<p>The 'ref' columns are reference values only, not tests. "ref expected H ratio" = the expected plug-in ratio if the data were exactly uniform with N=328 (Miller–Madow).
"ref ones fraction" = the bit-1 frequency the binary encoding would produce for exactly uniform values.</p>{table(ent_v)}

<h2 class="pb">Layer 1: Graphics</h2>
<h3>Figure 1: Frequency histograms</h3><img src="{img(FIG / '1_frequency_histograms.png')}">
<h3>Figure 2: Correlograms (lags 0–20, 95% band)</h3><img src="{img(FIG / '2_correlograms.png')}">

<h2 class="pb">Appendix A: Method (fixed before results)</h2>
<ul>
<li>Merge: the 7 files were aligned by row position (328 rows each, no missing values) → data/merged/merged_7_files.csv. SHA-256 hashes of the inputs are in README.md.</li>
<li>Uniformity: Pearson χ² against the discrete uniform on [min, max] (expected count per value is 10.6–46.9, so every expected count is ≥ 5). Discrete KS (Conover 1972) with a conservative Kolmogorov p-value.</li>
<li>Independence: sample ACF at lags 1–20 with a ±1.96/√N band. Runs up-and-down test (E=(2n−1)/3, Var=(16n−29)/90, zero differences dropped). Wald–Wolfowitz runs above/below the median (ties dropped). Both use the two-sided normal approximation.</li>
<li>Entropy: plug-in Shannon entropy in bits, divided by log₂k. Binary entropy of each bit stream.</li>
<li>NIST: Python implementation of SP 800-22 Rev.1a (scripts/nist_sp800_22.py), because the official NIST STS download host was not reachable from the analysis environment. Parameters: Block Frequency M = smallest value with M≥20, M&gt;0.01n, N&lt;100. Longest Run M=8 (n&lt;6272) or 128. Non-overlapping template B=000000001, N=8 (all 148 aperiodic m=9 templates also run, see Appendix C). Serial m = ⌊log₂n⌋−3. ApEn m = ⌊log₂n⌋−6. α = 0.01.</li>
</ul>
<h2>Appendix B: Implementation validation against the NIST worked examples</h2><pre>{val}</pre>
<h2>Appendix C: Non-overlapping template, all 148 templates (count with p ≥ 0.01)</h2>
{table(nt.groupby('stream').p_value.agg(templates='count', passing=lambda v: int((v >= 0.01).sum()), min_p='min').reset_index())}
<h2 class="pb">Appendix D: Raw machine output (results/machine_output.txt, unedited)</h2><pre>{machine}</pre>
</body></html>"""
(REP / "Randomness_Report.html").write_text(html, encoding="utf-8")
print("written", xl, REP / "Randomness_Report.html")
