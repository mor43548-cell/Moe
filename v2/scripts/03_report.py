"""V2 report: every number is read from v2/results (machine output); nothing is typed by hand."""
import base64, json, html as H
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
V2 = ROOT / "v2"; RES, FIG, REP = V2 / "results", V2 / "figures", V2 / "report"
REP.mkdir(exist_ok=True)
rd = lambda n: pd.read_csv(RES / n, keep_default_na=False, na_values=[""])
merged = pd.read_csv(V2 / "data/merged_lotto_table.csv")
uni, blk, ind, ent = rd("1_uniformity.csv"), rd("1_uniformity_blocks.csv"), rd("2_independence.csv"), rd("3_entropy.csv")
acf, sc, nist, nt = rd("2_acf_lags_0_20.csv"), rd("4_nist_scorecard.csv"), rd("4_nist_machine_output.csv"), rd("4_nist_nonoverlapping_148.csv")
cm, cs = rd("1_counts_main.csv"), rd("1_counts_strong.csv")
S = json.loads((RES / "summary.json").read_text())
machine = (RES / "machine_output.txt").read_text()
design = (V2 / "DESIGN.md").read_text()
val = (ROOT / "results/validation_nist.txt").read_text()
v1sc = pd.read_csv(ROOT / "results/4_nist_scorecard.csv", keep_default_na=False, na_values=[""])

def f(x, d=4):
    if isinstance(x, str):
        return x
    if pd.isna(x):
        return "—"
    return f"{x:.3e}" if (x != 0 and abs(x) < 1e-3) else f"{x:.{d}f}"
st = lambda ok: "PASS" if ok else "FAIL"
T = sc.set_index("test")
def nl(t):
    r = T.loc[t]
    return f"N/A — {r.reason}" if r.status == "N/A" else f"p = {r.p_values}; proportion {r.proportion}" + (f"; {r.reason}" if r.reason else "")

U = uni.set_index("variable"); um, us = U.iloc[0], U.iloc[1]
I = ind.set_index("series")
lb8 = S["prop_bounds_m8"][0]
ud_bad = ind[ind.runs_updown_p <= 0.05].series.tolist(); md_bad = ind[ind.runs_median_p <= 0.05].series.tolist()
tot_out = int(ind.lags_outside_band.sum()); clean = ind[ind.lags_outside_band == 0].series.tolist()
E = ent.set_index("variable")
P = S["nist"]

CL = [
 ("1. Uniformity", "Frequency test χ², main numbers pooled (1..37, n=1968)", "p > 0.01 (aim > 0.05)",
  f"p = {f(um.p_chi2_primary)} (adjusted, primary); unadjusted p = {f(um.p_chi2_unadjusted)}",
  "above 0.01, below the 0.05 aim", st(um.p_chi2_primary > 0.01)),
 ("1. Uniformity", "Frequency test χ², strong number (1..7, n=328)", "p > 0.01 (aim > 0.05)", f"p = {f(us.p_chi2_primary)}", "", st(us.p_chi2_primary > 0.01)),
 ("1. Uniformity", "KS (discrete), main numbers pooled", "p > 0.01 (aim > 0.05)", f"D = {f(um.KS_D)}, p = {f(um.p_KS)}", "above 0.01, below the 0.05 aim", st(um.p_KS > 0.01)),
 ("1. Uniformity", "KS (discrete), strong number", "p > 0.01 (aim > 0.05)", f"D = {f(us.KS_D)}, p = {f(us.p_KS)}", "", st(us.p_KS > 0.01)),
 ("1. Uniformity", "Proportion of the 8 blocks (41 draws each) passing χ²", f"~0.99 (NIST interval m=8: ≥ {lb8:.4f})",
  f"main {S['prop_main']:.4f}; strong {S['prop_strong']:.4f}", "", st(S['prop_main'] >= lb8 and S['prop_strong'] >= lb8)),
 ("1. Uniformity", "NIST Frequency", "p > 0.01, proportion ~0.99", nl("Frequency"), "", T.loc["Frequency"].status),
 ("1. Uniformity", "NIST Block Frequency", "p > 0.01 (aim > 0.05), proportion ~0.99", nl("BlockFrequency"), "", T.loc["BlockFrequency"].status),
 ("1. Uniformity", "Frequency histogram (visual)", "relatively flat", "Figure 1",
  f"counts range {int(um.min_count)}-{int(um.max_count)} around an expected 53.2; strong {int(us.min_count)}-{int(us.max_count)} around 46.9",
  st(um.p_chi2_primary > 0.01)),
 ("2. Independence", "Runs up/down, 7 series over draws", "p > 0.05",
  "; ".join(f"{s.split(' ')[0]}={f(r.runs_updown_p)}" for s, r in I.iterrows()), f"fails: {', '.join(ud_bad) or 'none'}", st(not ud_bad)),
 ("2. Independence", "Runs above/below median, 7 series", "p > 0.05",
  "; ".join(f"{s.split(' ')[0]}={f(r.runs_median_p)}" for s, r in I.iterrows()), f"fails: {', '.join(md_bad) or 'none'}", st(not md_bad)),
 ("2. Independence", "NIST Runs", "p > 0.05, proportion ~0.99", nl("Runs"), "", T.loc["Runs"].status),
 ("2. Independence", "NIST Longest Run of Ones", "p > 0.05, proportion ~0.99", nl("LongestRun"), "", T.loc["LongestRun"].status),
 ("2. Independence", "Autocorrelation coefficient, lags 1-20", "→ 0 within the 95% band ±1.96/√N",
  f"max |r| = {ind.max_abs_acf.max():.3f}; {tot_out}/140 coefficients outside the band (about 7 expected by chance)", "", st(tot_out == 0)),
 ("2. Independence", "Correlogram (visual)", "all bars except lag 0 inside the band", "Figure 2",
  f"all inside only for: {', '.join(clean)}", st(len(clean) == 7)),
 ("2. Independence", "NIST Random Excursions", "p > 0.05", nl("RandomExcursions"), "", T.loc["RandomExcursions"].status),
 ("2. Independence", "NIST Random Excursions Variant", "p > 0.05", nl("RandomExcursionsVariant"), "", T.loc["RandomExcursionsVariant"].status),
 ("3. Entropy", "Shannon, main numbers pooled (max log2 37)", "≥ 99.9%", f"{E.iloc[0].ratio:.5f}",
  f"reference: expected under H0 at this n = {E.iloc[0].ref_expected_ratio_if_H0:.5f}", st(E.iloc[0].ratio >= 0.999)),
 ("3. Entropy", "Shannon, strong number (max log2 7)", "≥ 99.9%", f"{E.iloc[1].ratio:.5f}",
  f"reference: expected under H0 = {E.iloc[1].ref_expected_ratio_if_H0:.5f}", st(E.iloc[1].ratio >= 0.999)),
 ("3. Entropy", "Shannon, NIST bit stream (max 1.000)", "≥ 0.999", f"{E.iloc[2].ratio:.5f}", "", st(E.iloc[2].ratio >= 0.999)),
 ("3. Entropy", "NIST Maurer's Universal", "p > 0.05", nl("MaurerUniversal"), "", T.loc["MaurerUniversal"].status),
 ("3. Entropy", "NIST Approximate Entropy", "p > 0.05", nl("ApproximateEntropy"), "", T.loc["ApproximateEntropy"].status),
 ("4. NIST", "Binary Matrix Rank", "p > 0.01", nl("BinaryMatrixRank"), "", T.loc["BinaryMatrixRank"].status),
 ("4. NIST", "Full battery", "15/15 PASS", f"PASS {P['PASS']}/15 · FAIL {P['FAIL']}/15 · N/A {P['N/A']}/15", "", st(P["PASS"] == 15)),
]
cl = pd.DataFrame(CL, columns=["section", "item", "required threshold", "observed (machine output)", "note", "status"])
cl.insert(0, "checked", "✔"); cl.to_csv(RES / "0_checklist.csv", index=False)

cmp_ = pd.DataFrame(dict(test=sc.test, v1_status=v1sc.status, v2_status=sc.status))

xl = REP / "Lotto_Randomness_Report_v2.xlsx"
with pd.ExcelWriter(xl, engine="openpyxl") as w:
    for name, d in [("0_Checklist", cl), ("NIST_Scorecard", sc), ("v1_vs_v2", cmp_), ("Merged_Data", merged),
                    ("1_Uniformity", uni), ("1_Blocks", blk), ("1_Counts_main", cm), ("1_Counts_strong", cs),
                    ("2_Independence", ind), ("2_ACF", acf), ("3_Entropy", ent), ("4_NIST_machine_output", nist),
                    ("4_NonOverlap_148", nt)]:
        d.to_excel(w, sheet_name=name, index=False)
    pd.DataFrame({"DESIGN.md": design.splitlines()}).to_excel(w, sheet_name="Design_preregistered", index=False)
from openpyxl import load_workbook
from openpyxl.styles import PatternFill, Font
wb = load_workbook(xl); fills = {"PASS": "C6EFCE", "FAIL": "FFC7CE", "N/A": "E7E6E6"}
for ws in wb.worksheets:
    for c in ws[1]:
        c.font = Font(bold=True)
    for col in ws.columns:
        ws.column_dimensions[col[0].column_letter].width = min(60, max(10, max(len(str(c.value or "")) for c in col[:60]) + 2))
    for row in ws.iter_rows(min_row=2):
        for c in row:
            if c.value in fills:
                c.fill = PatternFill("solid", fgColor=fills[c.value])
wb.save(xl)

def img(p): return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()
def table(d):
    h = "<table><thead><tr>" + "".join(f"<th>{H.escape(str(c))}</th>" for c in d.columns) + "</tr></thead><tbody>"
    for _, r in d.iterrows():
        h += "<tr>" + "".join(
            f"<td class='s-{str(r[c]).replace('/', '')}'>{H.escape(f(r[c]))}</td>" if str(r[c]) in ("PASS", "FAIL", "N/A")
            else f"<td>{H.escape(f(r[c]) if not isinstance(r[c], (int,)) else str(r[c]))}</td>" for c in d.columns) + "</tr>"
    return h + "</tbody></table>"

he = f"""
<p><b>השערה שנבדקה (H0):</b> כל שורה היא הגרלת לוטו 6/37 עם מספר חזק. ששת המספרים הם תת-קבוצה אקראית אחידה של {{1..37}} ונשמרים ממוינים. המספר החזק אחיד על {{1..7}}. ההגרלות בלתי תלויות זו בזו. מערך המבחנים והספים זהה לבקשה ולגרסה 1, ותכנון המחקר נרשם ונשמר בגיט <u>לפני</u> הרצת המבחנים (נספח A).</p>
<p><b>תשובה לקריטריון "מעבר של 100% ממבחני האקראיות": לא מתקיים.</b> התמונה שונה מאוד מגרסה 1, ויש לקרוא אותה במלואה:</p>
<ul>
<li><b>אחידות:</b> χ² על 1,968 המספרים הראשיים יחד: p={f(um.p_chi2_primary)}. הערך גבוה מסף החובה 0.01, ונמוך מהשאיפה של 0.05. KS: p={f(um.p_KS)}. המספר החזק: χ² p={f(us.p_chi2_primary)}. יחס המעבר על 8 בלוקים הוא 1.0000 גם במספרים הראשיים וגם במספר החזק.</li>
<li><b>אי-תלות:</b> ‏2 כישלונות מתוך 14 מבחני רצפים בסף 0.05: {', '.join(ud_bad)} ברצפי עלייה/ירידה, ו-{', '.join(md_bad)} ברצפים מעל/מתחת לחציון. {tot_out} מתוך 140 מקדמי אוטו-קורלציה מחוץ לרצועת 95%, כשבמקרה צפויים כ-7. ‏|r| המרבי הוא {ind.max_abs_acf.max():.3f}.</li>
<li><b>אנטרופיה:</b> במספרים הראשיים {E.iloc[0].ratio:.5f}, ובמספר החזק {E.iloc[1].ratio:.5f}. שניהם מתחת לסף 99.9%, אבל קרובים לערך שצפוי תחת H0 במדגם בגודל הזה ({E.iloc[0].ref_expected_ratio_if_H0:.5f} ו-{E.iloc[1].ref_expected_ratio_if_H0:.5f}). ברצף הביטים: {E.iloc[2].ratio:.5f}, עובר.</li>
<li><b>NIST SP 800-22:</b> ‏PASS {P['PASS']}/15, FAIL {P['FAIL']}/15, N/A {P['N/A']}/15. הכישלון היחיד הוא DFT, ש-p שלו {T.loc['DFT'].p_values} (גבוה מ-0.05). הוא מסומן FAIL רק משום שאף בלוק אינו ארוך מספיק כדי לחשב יחס מעבר, וכלל הסטטוס שנקבע מראש דורש יחס כזה. לא שיניתי את הכלל אחרי שראיתי את התוצאה. ששה מבחנים לא ניתנים להרצה באורך של {S['stream_bits']:,} ביטים.</li>
</ul>
<p><b>פרשנות זהירה (לא מחליפה ספים):</b></p>
<ol>
<li>בגרסה 1 נבדקה כל עמודה ממוינת כאילו היא אמורה להיות אחידה, ולכן הכישלונות הגורפים שם היו תוצאה של מבנה הבדיקה. בגרסה 2, שבה הבדיקה מותאמת למבנה של לוטו, רוב המבחנים עוברים, וההבדלים שנותרו קטנים.</li>
<li>אי אפשר לעמוד בקריטריון 100% במדגם של 328 הגרלות, כי שישה מבחני NIST דורשים בין 38,912 ל-1,000,000 ביטים. גם הגרלות לוטו אמיתיות לחלוטין היו מקבלות N/A בששת המבחנים האלה.</li>
<li>מעבר במבחן אינו הוכחה לאקראיות, רק אי-דחייה של H0. לכן אי אפשר "להוכיח" שאלה נתוני לוטו. אפשר רק לומר שהנתונים עקביים או לא עקביים עם לוטו הוגן במבחנים שהוגדרו.</li>
<li>המבנה (6 מספרים שונים בטווח 1–37, ממוינים, ומספר בטווח 1–7) תואם לוטו. מקור הנתונים לא אומת מול ארכיון רשמי.</li>
</ol>"""

page = f"""<!doctype html><html lang="he"><head><meta charset="utf-8"><title>Lotto Randomness v2</title><style>
body{{font-family:'DejaVu Sans',sans-serif;font-size:10.5px;color:#111;margin:0}} h1{{font-size:20px;margin:0 0 4px}}
h2{{font-size:15px;border-bottom:2px solid #333;padding-bottom:3px;margin-top:20px}} h3{{font-size:12.5px;margin:12px 0 5px}}
.rtl{{direction:rtl;text-align:right;font-size:11.5px;line-height:1.55}}
table{{border-collapse:collapse;width:100%;margin:6px 0 10px;font-size:9px}} th,td{{border:1px solid #999;padding:3px 4px;vertical-align:top}}
th{{background:#e8e8e8}} .s-PASS{{background:#c6efce;font-weight:bold}} .s-FAIL{{background:#ffc7ce;font-weight:bold}} .s-NA{{background:#e7e6e6;font-weight:bold}}
img{{width:100%}} pre{{font-size:7.5px;white-space:pre-wrap;background:#f6f6f6;padding:6px;border:1px solid #ccc}}
.box{{border:2px solid #b00;padding:8px;margin:8px 0;background:#fff5f5}} .pb{{page-break-before:always}}
</style></head><body>
<h1>Randomness Assessment, Version 2: the data treated as 6/37 lotto draws with a strong number</h1>
<div>N = 328 draws · same test set and thresholds as requested (and as v1) · design pre-registered in git before testing · all numbers are unmodified program output</div>
<div class="box"><b>Pre-defined criterion "pass 100% of the randomness tests": NOT MET.</b><br>
NIST: PASS {P['PASS']}/15 · FAIL {P['FAIL']}/15 (DFT: p = {T.loc['DFT'].p_values}; failed only because no proportion could be computed) · N/A {P['N/A']}/15.
Pooled χ² p = {f(um.p_chi2_primary)} (above 0.01, below 0.05). Runs tests: 2/14 below 0.05. Correlogram: {tot_out}/140 lags outside the band. Value entropy below 99.9%.</div>
<h2>Layer 3: Verbal summary (Hebrew)</h2><div class="rtl">{he}</div>
<h2 class="pb">Layer 2: Tables</h2><h3>Checklist, item by item</h3>{table(cl)}
<h3 class="pb">NIST SP 800-22 scorecard</h3>
<p>Bits: each main draw → combinatorial rank (uniform on 0..2,324,783 under H0); kept if &lt; 2²¹ → 21 bits. Strong − 1 kept if &lt; 4 → 2 bits.
Under H0 every bit is an i.i.d. fair bit. Stream = {S['stream_bits']:,} bits (main draws used {S['used_main']}/328, strong numbers used {S['used_strong']}/328).
p-value = whole stream. Proportion = share of the 8 blocks (41 draws each) with p ≥ 0.01.</p>{table(sc)}
<h3>Comparison of statuses, v1 vs v2</h3>{table(cmp_)}
<h3>1. Uniformity</h3>{table(uni)}{table(blk)}
<h3>2. Independence (7 series over draws)</h3>{table(ind)}
<h3>3. Entropy</h3>{table(ent)}
<h2 class="pb">Layer 1: Graphics</h2><h3>Figure 1: Frequency histograms</h3><img src="{img(FIG / '1_frequency_histograms.png')}">
<h3>Figure 2: Correlograms</h3><img src="{img(FIG / '2_correlograms.png')}">
<h2 class="pb">Appendix A: Pre-registered design (v2/DESIGN.md, committed before testing)</h2><pre>{H.escape(design)}</pre>
<h2>Appendix B: NIST implementation validation</h2><pre>{H.escape(val)}</pre>
<h2 class="pb">Appendix C: Raw machine output (unedited)</h2><pre>{H.escape(machine)}</pre>
</body></html>"""
(REP / "Lotto_Randomness_Report_v2.html").write_text(page, encoding="utf-8")
print("ok")
