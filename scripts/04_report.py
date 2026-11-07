"""Step 4: build the structured report (HTML + PDF) and the companion Excel workbook
from output/machine_output/03_results.json. No number is recomputed or altered here:
every figure shown is read from the machine output of step 3.
"""
import base64, json, math, pathlib, datetime
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
MO = ROOT / "output" / "machine_output"; FIG = ROOT / "output" / "figures"; REP = ROOT / "report"
REP.mkdir(exist_ok=True)
R = json.loads((MO / "03_results.json").read_text())
V = json.loads((MO / "02_validation.json").read_text())
CH = json.loads((ROOT / "data" / "merged" / "merge_checks.json").read_text())
tbl = pd.read_csv(ROOT / "data" / "merged" / "merged_table.csv", index_col="row")
COLS = R["columns"]; U = R["uniformity"]; I = R["independence"]; E = R["entropy"]
NI = {r["name"]: r for r in R["nist"]}
BE = R["bit_encoding"]

f4 = lambda p: "—" if p is None else (f"{p:.4f}" if p >= 1e-4 else f"{p:.2e}")
def verdict(ok):
    return "עומד בסף" if ok else "לא עומד בסף"

# -------- derived descriptive references (arithmetic only, no new test) --------
K, N = 37, len(tbl) * 6
mm_ratio = (math.log2(K) - (K - 1) / (2 * N * math.log(2))) / math.log2(K)       # Miller-Madow
p1 = BE["ones_fraction_if_numbers_perfectly_uniform_1_37"]
hb_enc = -(p1 * math.log2(p1) + (1 - p1) * math.log2(1 - p1))
obs = U["pooled"]["chi2"]["observed"]; e_ = U["pooled"]["chi2"]["expected_per_value"]
band = 1.96 * math.sqrt(e_ * (1 - 1 / K))
out_band = [i + 1 for i, o in enumerate(obs) if abs(o - e_) > band]
col_out = {c: I[c]["lags_outside_band"] for c in COLS}
n_out = sum(len(v) for v in col_out.values())

# -------- checklist (section-by-section, exactly as specified) -------------------
CL = []
def add(sec, item, thr, val, ok, note=""):
    CL.append(dict(section=sec, item=item, threshold=thr, measured=val,
                   status=("N/A" if ok is None else "PASS" if ok else "FAIL"), note=note))

S1 = "1. אחידות (Uniformity)"
c, k = U["pooled"]["chi2"], U["pooled"]["ks"]
add(S1, "מבחן תדירות χ² — כל 1,968 המספרים (1..37)", "p > 0.05", f"χ²={c['chi2']:.3f}, df=36, p={f4(c['p_value'])}", c["p_value"] > 0.05)
add(S1, "מבחן KS (בדיד) — כל 1,968 המספרים", "p > 0.05", f"D={k['D']:.4f}, p={f4(k['p_value'])}", k["p_value"] > 0.05)
for col in COLS:
    add(S1, f"מבחן תדירות χ² — עמודה {col}", "p > 0.05", f"χ²={U[col]['chi2']['chi2']:.1f}, p={f4(U[col]['chi2']['p_value'])}", U[col]["chi2"]["p_value"] > 0.05,
        "כל עמודה היא מקום קבוע ברשימה ממוינת (ראו §מבנה הנתונים)")
    add(S1, f"מבחן KS — עמודה {col}", "p > 0.05", f"D={U[col]['ks']['D']:.4f}, p={f4(U[col]['ks']['p_value'])}", U[col]["ks"]["p_value"] > 0.05)
for nm in ["Frequency (Monobit)", "Block Frequency"]:
    r = NI[nm]
    add(S1, f"NIST {nm}", "p > 0.01 (שאיפה > 0.05); Proportion ≈ 0.99; PASS", f"p={f4(r['min_p'])}, Proportion={r['proportion']:.0f}/1", r["status"] == "PASS")
add(S1, "היסטוגרמת תדרים (ויזואלי)", "עמודות שטוחות יחסית סביב הצפי",
    f"{len(out_band)}/37 עמודות מחוץ לרצועת ±1.96σ (מספרים {', '.join(map(str, out_band))}); טווח {min(obs)}–{max(obs)} מול צפי {e_:.2f}", None,
    "סף ויזואלי — מוצג בגרף 1 ללא הכרעה מספרית")

S2 = "2. אי-תלות (Independence)"
for col in COLS:
    add(S2, f"אוטו-קורלציה — עמודה {col} (lags 1–20)", "מקדמים ≈ 0, בתוך ±1.96/√n",
        f"max|r|={I[col]['max_abs_acf']:.3f}; מחוץ לרצועה: {col_out[col] or 'אין'}", len(col_out[col]) == 0)
add(S2, "אוטו-קורלציה — רצף שורה-אחר-שורה (lags 1–30)", "מקדמים ≈ 0, בתוך ±1.96/√n",
    f"max|r|={I['stream']['max_abs_acf']:.3f}; {len(I['stream']['lags_outside_band'])}/30 מחוץ לרצועה", len(I['stream']['lags_outside_band']) == 0,
    "נובע מהמיון בתוך כל שורה")
for col in COLS + ["stream"]:
    lbl = "רצף שורה-אחר-שורה" if col == "stream" else f"עמודה {col}"
    ru, rm = I[col]["runs_updown"], I[col]["runs_median"]
    add(S2, f"מבחן רצפים — עליות/ירידות — {lbl}", "p > 0.05", f"runs={ru['runs']} (צפי {ru['expected']:.1f}), p={f4(ru['p_value'])}", ru["p_value"] > 0.05)
    add(S2, f"מבחן רצפים — מעל/מתחת לחציון — {lbl}", "p > 0.05", f"runs={rm['runs']} (צפי {rm['expected']:.1f}), p={f4(rm['p_value'])}", rm["p_value"] > 0.05)
    if "runs_identical" in I[col]:
        ri = I[col]["runs_identical"]
        add(S2, f"רצפים של ערכים זהים עוקבים — {lbl}", "p > 0.05", f"זוגות={ri['consecutive_identical_pairs']} (צפי {ri['expected']:.1f}), p={f4(ri['p_value'])}", ri["p_value"] > 0.05)
for nm in ["Runs", "Longest Run of Ones", "Random Excursions", "Random Excursions Variant"]:
    r = NI[nm]
    add(S2, f"NIST {nm}", "p > 0.05 / > 0.01; Proportion ≈ 0.99; PASS",
        r["status"] if not r["applicable"] else f"p={f4(r['min_p'])}, Proportion={r['proportion']:.0f}/1",
        None if not r["applicable"] else r["status"] == "PASS", r["note"])

S3 = "3. אנטרופיה (Entropy)"
add(S3, "אנטרופיית שאנון — כל 1,968 המספרים", "≥ 99.9% מ-log₂37", f"H={E['pooled']['H_bits']:.5f} / {E['pooled']['H_max']:.5f} = {E['pooled']['ratio']*100:.3f}%",
    E["pooled"]["ratio"] >= 0.999, f"לשם השוואה: התוחלת המקורבת (Miller–Madow) למקור אחיד מושלם ב-n=1,968 היא {mm_ratio*100:.3f}%")
for col in COLS:
    add(S3, f"אנטרופיית שאנון — עמודה {col}", "≥ 99.9% מ-log₂37", f"{E[col]['ratio']*100:.2f}%", E[col]["ratio"] >= 0.999)
add(S3, "אנטרופיית שאנון — רצף הביטים (11,808)", "≥ 0.999 (ביטים)", f"{E['bits']['H_bits']:.5f}", E["bits"]["ratio"] >= 0.999,
    f"תחת קידוד 6 ביט של מקור 1..37 אחיד מושלם התוחלת היא {hb_enc:.5f}")
for nm in ["Maurer's Universal Statistical", "Approximate Entropy"]:
    r = NI[nm]
    add(S3, f"NIST {nm}", "p > 0.01; Proportion ≈ 0.99; PASS",
        r["status"] if not r["applicable"] else f"p={f4(r['min_p'])}, Proportion={r['proportion']:.0f}/1",
        None if not r["applicable"] else r["status"] == "PASS", r["note"])

S4 = "4. סוללת NIST SP 800-22"
for r in R["nist"]:
    add(S4, f"NIST {r['name']}", "כל p ≥ 0.01; PASS",
        r["status"] if not r["applicable"] else f"min p={f4(r['min_p'])}; {r['n_pvalues_ge_0.01']}/{r['n_pvalues']} ערכי p ≥ 0.01",
        None if not r["applicable"] else r["status"] == "PASS", r["note"])
n_pass = sum(r["status"] == "PASS" for r in R["nist"])
n_fail = sum(r["status"] == "FAIL" for r in R["nist"])
n_na = sum(r["status"] == "NOT APPLICABLE" for r in R["nist"])
add(S4, "עמידה ב-100% (15/15) מהסוללה", "15/15 PASS", f"{n_pass} PASS · {n_fail} FAIL · {n_na} NOT APPLICABLE", n_pass == 15)

cl = pd.DataFrame(CL)

# -------- scorecard --------------------------------------------------------------
sc = pd.DataFrame([dict(
    test=r["name"], params=json.dumps(r["params"], ensure_ascii=False),
    p_value=("—" if not r["p_values"] else f4(r["p_values"][0]) if len(r["p_values"]) == 1 else f"min {f4(r['min_p'])}"),
    n_subtests=r["n_pvalues"], n_ge_001=r["n_pvalues_ge_0.01"], n_gt_005=r["n_pvalues_gt_0.05"],
    proportion=("—" if r["proportion"] is None else f"{int(r['proportion'])}/1"),
    status=r["status"], note=r["note"]) for r in R["nist"]])

# -------- Excel --------------------------------------------------------------------
allp = []
for r in R["nist"]:
    for i, p in enumerate(r["p_values"]):
        allp.append(dict(test=r["name"], index=i + 1, p_value=p, ge_0_01=p >= 0.01))
xl = REP / "randomness_report.xlsx"
with pd.ExcelWriter(xl, engine="openpyxl") as w:
    pd.DataFrame({"field": list(CH.keys()), "value": [json.dumps(v) if isinstance(v, dict) else v for v in CH.values()]}).to_excel(w, sheet_name="0_data_checks", index=False)
    tbl.to_excel(w, sheet_name="1_merged_table")
    cl.to_excel(w, sheet_name="2_checklist", index=False)
    sc.to_excel(w, sheet_name="3_NIST_scorecard", index=False)
    pd.DataFrame(allp).to_excel(w, sheet_name="4_NIST_all_pvalues", index=False)
    pd.DataFrame([dict(object=k, chi2=U[k]["chi2"]["chi2"], df=36, chi2_p=U[k]["chi2"]["p_value"], KS_D=U[k]["ks"]["D"], KS_p=U[k]["ks"]["p_value"],
                       min_count=U[k]["chi2"]["min_obs"], max_count=U[k]["chi2"]["max_obs"]) for k in ["pooled"] + COLS]).to_excel(w, sheet_name="5_uniformity", index=False)
    pd.DataFrame({"number": range(1, 38), **{k: U[k]["chi2"]["observed"] for k in ["pooled"] + COLS}}).to_excel(w, sheet_name="6_frequencies", index=False)
    acf = pd.DataFrame({k: I[k]["acf"][:21] for k in COLS}); acf.index.name = "lag"; acf.to_excel(w, sheet_name="7_acf_columns")
    pd.DataFrame({"lag": range(31), "acf_stream": I["stream"]["acf"]}).to_excel(w, sheet_name="8_acf_stream", index=False)
    pd.DataFrame([dict(object=k, **{f"updown_{a}": b for a, b in I[k]["runs_updown"].items()}, **{f"median_{a}": b for a, b in I[k]["runs_median"].items()},
                       **{f"identical_{a}": b for a, b in I[k].get("runs_identical", {}).items()}) for k in COLS + ["stream"]]).to_excel(w, sheet_name="9_runs", index=False)
    pd.DataFrame([dict(object=k, **v) for k, v in E.items()]).to_excel(w, sheet_name="10_entropy", index=False)
    pd.DataFrame(V).to_excel(w, sheet_name="11_impl_validation", index=False)

# -------- HTML ---------------------------------------------------------------------
img = lambda f: "data:image/png;base64," + base64.b64encode((FIG / f).read_bytes()).decode()
badge = lambda s: {"PASS": '<span class="b pass">PASS</span>', "FAIL": '<span class="b fail">FAIL</span>'}.get(s, f'<span class="b na">{s}</span>')
def rows_for(sec):
    return "".join(f"<tr><td>{x['item']}</td><td class='ltr'>{x['threshold']}</td><td class='ltr'>{x['measured']}</td><td>{badge(x['status'])}</td><td class='note'>{x['note']}</td></tr>"
                   for x in CL if x["section"] == sec)
sec_tbl = lambda sec: f"<table class='cl'><thead><tr><th>פריט</th><th>סף שנקבע</th><th>ערך שנמדד (פלט מכונה)</th><th>סטטוס</th><th>הערה</th></tr></thead><tbody>{rows_for(sec)}</tbody></table>"
sc_rows = "".join(f"<tr><td class='ltr'>{i+1}</td><td class='ltr'>{r['test']}</td><td class='ltr'>{r['params']}</td><td class='ltr'>{r['p_value']}</td>"
                  f"<td class='ltr'>{r['n_ge_001']}/{r['n_subtests']}</td><td class='ltr'>{r['proportion']}</td><td>{badge(r['status'])}</td><td class='note ltr'>{r['note']}</td></tr>"
                  for i, r in sc.iterrows())
val_ok = sum(1 for v in V if v["match"] is True); val_n = sum(1 for v in V if v["match"] is not None)
val_bad = [v["case"] for v in V if v["match"] is False]
n_cl = sum(1 for x in CL if x["status"] != "N/A"); n_cl_pass = sum(1 for x in CL if x["status"] == "PASS")

html = f"""<!doctype html><html lang="he" dir="rtl"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>דו״ח בדיקת אקראיות</title>
<style>
:root{{--fg:#1a1a1a;--muted:#555;--line:#d0d0d0;--bg:#fff;--head:#f2f4f7;--pass:#1e7b34;--fail:#b3261e;--na:#6b6b6b}}
body{{font-family:"Noto Sans Hebrew","Arial",sans-serif;color:var(--fg);background:var(--bg);margin:0 auto;max-width:1000px;padding:24px 16px;line-height:1.5;font-size:13px}}
h1{{font-size:22px;margin:0 0 4px}} h2{{font-size:17px;border-bottom:2px solid var(--fg);padding-bottom:3px;margin-top:28px}} h3{{font-size:14px;margin:18px 0 6px}}
.meta{{color:var(--muted);font-size:12px}} .ltr{{direction:ltr;unicode-bidi:embed;text-align:left}}
table{{border-collapse:collapse;width:100%;margin:8px 0;font-size:11.5px}} th,td{{border:1px solid var(--line);padding:4px 6px;vertical-align:top}} th{{background:var(--head)}}
.note{{color:var(--muted);font-size:10.5px}} .b{{font-weight:700;padding:1px 6px;border-radius:3px;color:#fff;font-size:10.5px;white-space:nowrap}}
.pass{{background:var(--pass)}} .fail{{background:var(--fail)}} .na{{background:var(--na)}}
.box{{border:1px solid var(--fg);padding:10px 14px;margin:10px 0;background:#fafafa}} img{{max-width:100%;border:1px solid var(--line)}}
pre{{direction:ltr;text-align:left;background:#f6f6f6;padding:8px;font-size:10px;white-space:pre-wrap;word-break:break-all}}
@media print{{h2{{break-before:auto}} table,img{{break-inside:avoid}}}}
</style></head><body>
<h1>דו״ח בדיקת אקראיות — טבלת 6 העמודות (Y, T, W, E, R, Q)</h1>
<div class="meta">הופק: {datetime.date.today().isoformat()} · כל המספרים בדו״ח נקראים ישירות מפלט המכונה (<span class="ltr">output/machine_output/03_results.json</span>) · קוד ונתונים: תיקיית <span class="ltr">scripts/</span></div>

<div class="box"><b>שורה תחתונה:</b> הקריטריון שהוגדר — מעבר של 100% ממבחני האקראיות (15/15 בסוללת NIST) — <b>לא מתקיים</b>.
בסוללת NIST: {n_pass} PASS, {n_fail} FAIL, {n_na} NOT APPLICABLE. ברשימת הבדיקה המלאה: {n_cl_pass} מתוך {n_cl} פריטים הניתנים להכרעה עומדים בסף.
מבחן התדירות על כל 1,968 המספרים (χ², p={f4(U['pooled']['chi2']['p_value'])}) ומבחן KS (p={f4(U['pooled']['ks']['p_value'])}) עומדים בסף p&gt;0.05;
אנטרופיית שאנון ({E['pooled']['ratio']*100:.3f}%) נמוכה מסף 99.9%.</div>

<h2>0. הנתונים</h2>
<ul>
<li><b>התקבלו {CH['files_received']} קבצים</b> (Y, T, W, E, R, Q), ולא 7 כפי שצוין בבקשה. הניתוח בוצע על 6 הקבצים שהתקבלו בלבד.</li>
<li>טבלה מאוחדת: {CH['rows']} שורות × 6 עמודות, ללא ערכים חסרים, כולם שלמים. ערך מינימלי {CH['global_min']}, מקסימלי {CH['global_max']}. הטבלה יוצאה לפני כל ניתוח: <span class="ltr">data/merged/merged_table.csv / .xlsx</span>.</li>
<li><b>מבנה הנתונים (תיאור, לא מבחן):</b> בכל {CH['rows_with_6_distinct_values']} השורות יש 6 ערכים שונים, ובכל {CH['rows_strictly_ascending_Q<W<E<R<T<Y']} השורות מתקיים Q&lt;W&lt;E&lt;R&lt;T&lt;Y. כלומר כל שורה היא קבוצה של 6 מספרים שונים מתוך 1..37 בסדר ממוין, וכל עמודה היא מקום קבוע ברשימה הממוינת (Q הקטן ביותר, Y הגדול ביותר). לכן התפלגות כל עמודה לבדה אינה אחידה מעצם המבנה; המושא הטבעי למבחן התדירות הוא כל 1,968 המספרים יחד. שתי הרמות מוצגות.</li>
<li>סדר השורות נלקח כפי שהוא בקבצים (מניחים שזה סדר ההגרלות; הקבצים אינם כוללים תאריכים).</li>
</ul>

<h2>1. שכבה גרפית</h2>
<h3>גרף 1 — היסטוגרמת תדרים, כל 1,968 המספרים</h3><img src="{img('fig1_histogram_pooled.png')}" alt="histogram">
<p class="note">{len(out_band)} מתוך 37 עמודות חורגות מרצועת ±1.96σ (מספרים {', '.join(map(str, out_band))}). תחת אחידות, הצפי הוא כ-{0.05*37:.1f} חריגות.</p>
<h3>גרף 2 — היסטוגרמת תדרים לכל עמודה</h3><img src="{img('fig2_histogram_columns.png')}" alt="per column histograms">
<h3>גרף 3 — קורלוגרמה לכל עמודה (lags 1–20)</h3><img src="{img('fig3_correlogram_columns.png')}" alt="correlograms">
<p class="note">סה״כ {n_out} מתוך 120 מקדמים חורגים מרצועת 95% (הצפי תחת אי-תלות: כ-6). max|r| = {max(I[c]['max_abs_acf'] for c in COLS):.3f}.</p>
<h3>גרף 4 — קורלוגרמה של הרצף שורה-אחר-שורה</h3><img src="{img('fig4_correlogram_stream.png')}" alt="stream correlogram">
<p class="note">כל 30 ה-lags חורגים מהרצועה. זו תוצאה של המיון בתוך כל שורה (מחזוריות של 6), והיא תקפה גם לרצף הביטים שנבנה מאותו סדר.</p>

<h2>2. שכבה טבלאית</h2>
<h3>2.1 Scorecard — NIST SP 800-22 (15 מבחנים)</h3>
<p>קלט: {BE['n_bits']:,} ביטים. כלל הקידוד: <span class="ltr">{BE['rule']}</span>. מספר רצפים m=1, ולכן Proportion = 1/1 או 0/1. רווח הקבלה של NIST ל-m=1 הוא
<span class="ltr">0.99 ± 3·√(0.99·0.01/1) = [{R['nist_proportion_interval_m1'][0]:.4f}, {R['nist_proportion_interval_m1'][1]:.4f}]</span>.
מבחן אחידות ערכי ה-p (P-value<sub>T</sub>) של NIST דורש לפחות 55 רצפים ולכן אינו מחושב.</p>
<table><thead><tr><th>#</th><th>מבחן</th><th>פרמטרים</th><th>p-value</th><th>תתי-מבחנים ≥0.01</th><th>Proportion</th><th>סטטוס</th><th>הערה</th></tr></thead><tbody>{sc_rows}</tbody></table>
<p class="note">שבר האחדות שנצפה: {BE['ones_fraction']:.4f}. השבר הצפוי כשמספרים 1..37 אחידים לחלוטין ומקודדים ב-6 ביטים: {p1:.4f} (ולא 0.5). Monobit ומבחנים התלויים בו בודקים P(1)=0.5 ולכן נכשלים מעצם הקידוד.</p>

<h3>2.2 רשימת בדיקה לפי סעיף</h3>
<h3>{S1}</h3>{sec_tbl(S1)}
<h3>{S2}</h3>{sec_tbl(S2)}
<h3>{S3}</h3>{sec_tbl(S3)}
<h3>{S4}</h3>{sec_tbl(S4)}

<h2>3. שכבה מילולית — סיכום לבודק חיצוני</h2>
<p><b>מטרה ומושא.</b> נבדקה טבלה של {CH['rows']} שורות × 6 עמודות (6 קבצים שהתקבלו; הבקשה ציינה 7). כל שורה מכילה 6 מספרים שונים מתוך 1..37 בסדר ממוין. המדדים, המבחנים והספים הם אלה שהוגדרו בבקשה, והם הוחלו כפי שהם.</p>
<p><b>אחידות.</b> על כל 1,968 המספרים: χ²=({U['pooled']['chi2']['chi2']:.3f}, df=36, p={f4(U['pooled']['chi2']['p_value'])}) ו-KS בדיד (D={U['pooled']['ks']['D']:.4f}, p={f4(U['pooled']['ks']['p_value'])}) — שניהם עומדים בסף p&gt;0.05.
ברמת העמודה הבודדת כל 12 מבחני ה-χ²/KS נכשלים (p&lt;10⁻¹⁵); זה צפוי, כי כל עמודה היא מקום קבוע ברשימה ממוינת. מבחני NIST Frequency ו-Block Frequency על רצף הביטים נכשלים (p≈0).</p>
<p><b>אי-תלות.</b> בקורלוגרמות לפי עמודה {n_out}/120 מקדמים חורגים מרצועת 95% (max|r|={max(I[c]['max_abs_acf'] for c in COLS):.3f}).
מבחן הרצפים (עליות/ירידות) נכשל בעמודה Q (p={f4(I['Q']['runs_updown']['p_value'])}); מבחן הרצפים סביב החציון נכשל בעמודה R (p={f4(I['R']['runs_median']['p_value'])}); כל שאר מבחני הרצפים (כולל ערכים זהים עוקבים) עומדים בסף p&gt;0.05.
ברצף שורה-אחר-שורה כל ה-lags חורגים, בגלל המיון בתוך השורה. NIST Runs ו-Longest Run נכשלים. Random Excursions ו-Variant אינם ישימים (J={NI['Random Excursions']['stats']['J']} מחזורים, והמפרט דורש לפחות 500).</p>
<p><b>אנטרופיה.</b> אנטרופיית שאנון של כל המספרים היא {E['pooled']['ratio']*100:.3f}% מהמקסימום, מתחת לסף 99.9%. (לשם השוואה, התוחלת המקורבת למקור אחיד מושלם בגודל מדגם זה היא {mm_ratio*100:.3f}%, גם היא מתחת לסף.)
ברצף הביטים: {E['bits']['H_bits']:.5f}, מתחת לסף. Maurer אינו ישים (n קטן מ-387,840); Approximate Entropy נכשל.</p>
<p><b>NIST SP 800-22.</b> {n_pass}/15 PASS (Binary Matrix Rank, Linear Complexity — שניהם עם פחות בלוקים מהמומלץ), {n_fail} FAIL, {n_na} NOT APPLICABLE. אורך הקלט ({BE['n_bits']:,} ביטים) נמוך מהאורך המומלץ ב-NIST (10⁶) עבור חלק מהמבחנים.</p>
<p><b>מסקנה לפי הקריטריון שהוגדר.</b> הקריטריון "מעבר של 100% ממבחני האקראיות" <b>אינו מתקיים</b>, ולכן לפי הכלל שהוגדר הבדיקה אינה מוכיחה שהנתונים הם נתוני הגרלה.</p>
<p><b>גבולות ההיסק.</b> (א) מבחן סטטיסטי אינו יכול להוכיח אקראיות; מעבר משמעותו רק אי-דחייה. (ב) NIST SP 800-22 בודק את ההשערה "ביטים בלתי תלויים עם P(1)=½". מקור הגרלה הוגן של 6 מתוך 37, בקידוד 6 ביט ובמיון בתוך השורה, אינו מקיים השערה זו. שבר האחדות שנצפה ({BE['ones_fraction']:.4f}) קרוב לצפוי מהקידוד ({p1:.4f}). לכן כשלי NIST כאן אינם מבחינים בין הגרלה הוגנת להגרלה שאינה הוגנת. (ג) ריבוי מבחנים: ברמת 0.05, צפויות כ-5% דחיות גם תחת אקראיות מלאה.</p>

<h2>4. שיטות, תקינות המימוש ושחזור</h2>
<ul>
<li>התוכנה הרשמית של NIST (STS 2.1.2) לא הייתה זמינה, כי האתר csrc.nist.gov חסום ברשת של סביבת הריצה. ספריית <span class="ltr">nistrng</span> נבדקה ונפסלה, כי במבחן Random Excursions שלה יש שגיאת ספירה. לכן 15 המבחנים מומשו ישירות מהמפרט (<span class="ltr">scripts/sp800_22.py</span>) עם קבועי STS 2.1.2.</li>
<li>אימות המימוש (<span class="ltr">scripts/02_validate_sp800_22.py</span>): {val_ok}/{val_n} התאמות לדוגמאות המפרט ולערכי הייחוס של הרחבת e (10⁶ ביטים). אי-התאמות: {'; '.join(val_bad) if val_bad else 'אין'}. שתי הדוגמאות הקצרות של DFT המודפסות במפרט אינן משתחזרות, אך ערך הייחוס המלא של e עבור DFT משתחזר.</li>
<li>KS על נתונים בדידים: הסטטיסטי מחושב על נקודות התומך 1..37. <span class="ltr">scipy.stats.kstest</span> מניח התפלגות רציפה ומנפח את D כשיש ערכים חוזרים, ולכן לא שימש. ערך ה-p לפי קולמוגורוב שמרני במקרה בדיד.</li>
<li>מבחן עליות/ירידות: הפרשים אפסיים הושמטו. מבחן החציון: ערכים השווים לחציון הושמטו. מבחן ערכים זהים עוקבים: בינומי מדויק, כאשר ההסתברות לזוג זהה היא Σp̂² של העמודה.</li>
<li>פרמטרי NIST עבור n=11,808: Block Frequency M=128; Longest Run M=128; Non-overlapping m=9, N=8; Overlapping m=9, M=1032; Linear Complexity M=500; Serial m=10; ApEn m=7 (הערכים הגדולים ביותר שמותרים לפי תנאי המפרט).</li>
<li>שחזור: <span class="ltr">pip install numpy scipy pandas matplotlib openpyxl &amp;&amp; python scripts/01_merge.py &amp;&amp; python scripts/02_validate_sp800_22.py &amp;&amp; python scripts/03_analysis.py &amp;&amp; python scripts/04_report.py</span></li>
<li>טביעות SHA-256 של קבצי המקור: <span class="ltr">{'; '.join(f'{k}: {v[:16]}…' for k, v in CH['sha256_raw'].items())}</span></li>
</ul>
<h3>פלט מכונה מלא (step 3 stdout)</h3><pre>{(MO / '03_analysis_stdout.txt').read_text()}</pre>
</body></html>"""
(REP / "randomness_report.html").write_text(html, encoding="utf-8")
cl.to_csv(REP / "checklist.csv", index=False)
sc.to_csv(REP / "nist_scorecard.csv", index=False)
print("checklist:", n_cl_pass, "/", n_cl, "PASS ; NIST:", n_pass, n_fail, n_na)
print(sc[["test", "p_value", "n_ge_001", "n_subtests", "proportion", "status"]].to_string(index=False))
