## תשובה בקצרה

**הנתונים לא עברו 100% (15/15) ממבחני האקראיות.** בסוללת NIST המבחנים הסתיימו כך, והתוצאה זהה בסף 0.01 ובסף 0.05:
- **2 עברו** (Matrix Rank ו-Linear Complexity), אבל בשניהם אורך הנתונים קצר ממה שהתקן מחייב.
- **10 נכשלו.**
- **3 לא ניתנים לחישוב**, כי הנתונים קצרים מדי.

לכן לא ניתן להוכיח שאלה "נתוני הגרלה" לפי הקריטריון שהגדרת. גם לא ניתן לקבוע את ההפך. רוב הכישלונות נובעים מאופן הבנייה של הנתונים, ולא בהכרח ממקורם:
- הנתונים ממוינים בתוך כל שורה (Q<W<E<R<T<Y).
- המספרים 1–37 ממופים לביטים בקידוד בינארי של 6 ביטים, שאינו מאוזן מעצם הגדרתו.

במבחנים שנעשו ישירות על המספרים (ולא על הביטים) הנתונים עמדו ברוב הספים. היוצא מן הכלל הוא האנטרופיה.

הקובץ שוחזר וה-SHA-256 תואם: `b60bd7b7cd35bb708573189dfccb3f0c989890fea527f8b67b573fb89bb763b4` (328 שורות ועוד כותרת).

**לגבי "7 הקבצים":** קיבלתי טבלה אחת שכבר מאוחדת. אין לי דרך לוודא שהיא אכן מורכבת מ-7 קבצים.

---

## 1. צ׳קליסט לפי סעיפי הבקשה

| # | סעיף | בוצע | תוצאה |
|---|---|---|---|
| 0 | טבלה אחת מאוחדת ויצוא שלה | ✅ | `full_table_8_columns.csv`, ה-SHA תואם |
| 1.1 | מבחן תדירות χ² (בחרתי ב-χ² מבין "χ² או KS") | ✅ | מספרים 1–37: p=0.153282 → PASS@0.01, PASS@0.05. עמודת ABC: p=0.231916 → PASS@0.01, PASS@0.05 |
| 1.2 | NIST Frequency, כולל יחס מעבר | ✅ | p=0.000000, יחס מעבר 0/1 → FAIL בשני הספים |
| 1.3 | NIST Block Frequency | ✅ | p=0.000000, יחס מעבר 0/1 → FAIL בשני הספים |
| 1.4 | היסטוגרמת תדרים | ✅ | `fig1_frequency_histogram.png`. העמודות בין 39 ל-72 סביב תוחלת של 53.19 |
| 2.1 | אוטו-קורלציה וקורלוגרם | ✅ | `fig2_correlogram.png`. ברצף הממוין: כל 20 ה-lags מחוץ לרצועה. בעמודות הבודדות: 0–3 lags מחוץ לרצועת 95% |
| 2.2 | מבחן רצפים מורחב (עליות וירידות) | ✅ | ברצף הממוין: FAIL. בעמודה Q: FAIL בשני הספים. בעמודות W, E, R, T, Y, ABC: PASS בשני הספים |
| 2.3 | NIST Runs | ✅ | p=0.0, כי תנאי הקדם לא התקיים → FAIL |
| 2.4 | Longest Run of Ones | ✅ | p=0.000000 → FAIL |
| 2.5 | Random Excursions ו-Variant | ✅ (נוסה) | לא ישים: יש J=8 מחזורים, ו-NIST דורש לפחות 500 |
| 3.1 | אנטרופיית שאנון | ✅ | יחס למקסימום: 0.996930 במספרים, 0.993771 ב-ABC, 0.981636 בביטים → אף אחד מהם לא מגיע ל-0.999 |
| 3.2 | Maurer | ✅ (נוסה) | לא ניתן לחישוב: התקן דורש לפחות 387,840 ביטים, ויש 11,808 |
| 3.3 | Approximate Entropy | ✅ | p=0.000000 → FAIL |
| 4 | Binary Matrix Rank ושאר הסוללה (15 מבחנים) | ✅ | ראו את ה-Scorecard |
| 5 | שלוש שכבות: גרפית, טבלאית ומילולית, ודו"ח אקסל | ✅ | הגרפים, `nist_scorecard.csv`, `report.xlsx` ופרק 4 בהמשך |

**לא השתמשתי בערכים לדוגמה מהבקשה (0.4129, 0.9885 וכו׳) כערכי ייחוס.** הספים היחידים שהופעלו הם 0.01, 0.05 ו-0.999.

## 2. Scorecard (סוללת NIST SP 800-22)

לגבי הקידוד: כל מספר ראשי נכתב כמספר בינארי בן 6 ביטים, לפי סדר השורות ובכל שורה לפי הסדר Q, W, E, R, T, Y. כך מתקבלים n=11,808 ביטים, מתוכם 4,964 אחדות.

לגבי יחס המעבר: יש רק רצף אחד, ולכן היחס יכול להיות רק 0/1 או 1/1. לפי NIST, יחס מעבר בעל משמעות דורש הרבה רצפים.

| # | מבחן | p-value | יחס מעבר @0.01 | יחס מעבר @0.05 | סטטוס @0.01 | סטטוס @0.05 |
|---|---|---|---|---|---|---|
| 1 | Frequency (Monobit) | 0.000000 | 0/1 | 0/1 | FAIL | FAIL |
| 2 | Block Frequency | 0.000000 | 0/1 | 0/1 | FAIL | FAIL |
| 3 | Runs | 0.000000 (תנאי קדם לא התקיים) | 0/1 | 0/1 | FAIL | FAIL |
| 4 | Longest Run of Ones | 0.000000 | 0/1 | 0/1 | FAIL | FAIL |
| 5 | Binary Matrix Rank | 0.893303 (רק 11 מטריצות; NIST דורש לפחות 38) | 1/1 | 1/1 | PASS | PASS |
| 6 | DFT | 0.000000 | 0/1 | 0/1 | FAIL | FAIL |
| 7 | Non-overlapping Template (148 תבניות) | מינימום 2.59e-39, חציון 0.098666 | 108/148 | 90/148 | FAIL | FAIL |
| 8 | Overlapping Template | 0.001755 (רק 11 בלוקים) | 0/1 | 0/1 | FAIL | FAIL |
| 9 | Maurer's Universal | NA | NA | NA | NA | NA |
| 10 | Linear Complexity | 0.555031 (רק 23 בלוקים) | 1/1 | 1/1 | PASS | PASS |
| 11 | Serial (שני ערכי p) | 0.000000 / 0.000000 | 0/1 | 0/1 | FAIL | FAIL |
| 12 | Approximate Entropy | 0.000000 | 0/1 | 0/1 | FAIL | FAIL |
| 13 | Cumulative Sums (קדימה ואחורה) | 0.000000 / 0.000000 | 0/1 | 0/1 | FAIL | FAIL |
| 14 | Random Excursions | NA (J=8) | NA | NA | NA | NA |
| 15 | Random Excursions Variant | NA (J=8) | NA | NA | NA | NA |
| | **סיכום** | | | | **2 PASS, 10 FAIL, 3 NA** | **2 PASS, 10 FAIL, 3 NA** |

**כל ערכי ה-p שאינם מ-NIST:**
- χ² על המספרים 1–37: χ²=44.627, df=36, p=0.153282.
- χ² על ABC: χ²=8.085, df=6, p=0.231916.
- מבחן רצפים מורחב: ברצף הממוין z=−35.23, p≈0. בעמודות: Q 0.008635, W 0.087033, E 0.475464, R 0.349718, T 0.821850, Y 0.962530, ABC 0.706051.

## 3. שלוש הטבלאות הנדרשות

**3.1 טבלת חריגות**

| # | חריגה | ערך | סף שנחצה |
|---|---|---|---|
| 1 | רצף הממוין (Q..Y) מראה מחזוריות של 6: r=+0.7285 ב-lag 6, +0.7275 ב-lag 12, +0.7173 ב-lag 18 | כל 20 ה-lags מחוץ לרצועה | 95% ו-99% |
| 2 | ברצף הממוין יש 653 רצפי עליות וירידות, מול תוחלת של 1311.7 | z=−35.23 | 0.01 ו-0.05 |
| 3 | בזרם הביטים יש 4,964 אחדות מתוך 11,808 (0.4204) | S_n=−1880 | 0.01 ו-0.05 |
| 4 | אין אף רצף של 9 אחדות ברצף הביטים (Overlapping v=[11,0,0,0,0,0]), ובמבחן Longest Run אין אף בלוק עם רצף ארוך מ-6 | 15 תבניות קיבלו W=0 בכל 8 הבלוקים | 0.01 |
| 5 | בעמודה Q יש יותר רצפי עליות וירידות מהצפוי | z=+2.63, p=0.008635 | 0.01 ו-0.05 |
| 6 | אוטו-קורלציה בעמודות הבודדות: Y ב-lag 11 (r=−0.1610) מחוץ לרצועת 99%. מחוץ לרצועת 95% בלבד: Q lag 8, E lags 10/14/16, R lags 7/14, T lags 11/16, Y lag 6 | 11 חריגות מתוך 140 | 95% (אחת מהן גם 99%) |
| 7 | אנטרופיית המספרים 0.99693 ואנטרופיית ABC 0.99377 | מתחת ל-0.999 | 0.999 |
| 8 | תדירויות קיצון: המספר 21 הופיע 72 פעמים, המספר 12 הופיע 39 פעמים | תוחלת 53.19 | (χ² הכולל עבר) |
| 9 | בעמודת Day: שורה 153 היא שבת אחרי שבת (החזרה היחידה ברצף); הימים ראשון (7), רביעי (2) ושני (1) נדירים | | |
| 10 | Rank ו-Linear Complexity עברו, אבל עם 11 ו-23 בלוקים, הרבה מתחת לדרישת התקן | | |
| 11 | ב-Random Excursions יש רק J=8 מחזורים (הדרישה היא 500) | | |

**3.2 טבלת מעברים (יום ← היום של השורה הבאה, ספירות)**

| מ\אל | חמישי | ראשון | רביעי | שבת | שלישי | שני |
|---|---|---|---|---|---|---|
| חמישי | 0 | 1 | 0 | 4 | 40 | 0 |
| ראשון | 1 | 0 | 0 | 0 | 6 | 0 |
| רביעי | 0 | 0 | 0 | 2 | 0 | 0 |
| שבת | 44 | 0 | 2 | 1 | 89 | 1 |
| שלישי | 0 | 6 | 0 | 129 | 0 | 0 |
| שני | 0 | 0 | 0 | 1 | 0 | 0 |

מעבר בין הספים במבחן התבניות: 18 תבניות עוברות בסף 0.01 ונכשלות בסף 0.05 (108 מול 90). בכל שאר המבחנים הסטטוס זהה בשני הספים.

**3.3 תכונות, רכיבים ומידע לא מוסבר**

| רכיב | מה לא ידוע |
|---|---|
| עמודה ABC (1–7) | משמעותה. אולי "מספר חזק", אבל זה לא נמסר |
| עמודה Day | מה היא מייצגת, ולמה מופיעים ימים חריגים (ראשון, רביעי, שני) |
| שמות העמודות (אותיות מקלדת) וסדרן בכותרת: Y,T,W,R,ABC,Q,E | מה הן מייצגות, ולמה הסדר מעורבב |
| עמודת row | האם היא מייצגת סדר כרונולוגי, ובאיזה כיוון |
| "7 קבצים" | לא התקבלו כקבצים נפרדים |
| הערכים לדוגמה בבקשה (0.4129, 0.9885 וכו׳) | מקורם לא ידוע. לא השתמשתי בהם |
| סדר ההגרלה המקורי בתוך כל שורה | אבד במיון |

## 4. חמש הרשימות המפורשות

1. **דברים שאינני יכול להסביר:** את החריגות הבודדות באוטו-קורלציה (Q ב-lag 8, Y ב-lag 11) ואת עודף רצפי העליות והירידות בעמודה Q (p=0.0086). אלו יכולות להיות תנודות מקריות או אות אמיתי, והנתונים לא מאפשרים להבחין בין השניים. לא ביצעתי תיקון לריבוי השוואות, כפי שנדרש.
2. **דברים שאינני יודע מה משמעותם:** העמודות ABC, Day ו-row, שמות העמודות, וטענת "7 הקבצים".
3. **דברים שאינני מבין את הסיבתיות שלהם:** הימים החריגים, ושבת אחרי שבת בשורה 153. ייתכנו שינויי לוח זמנים או חגים, אבל אין לי מידע על כך.
   - לעומת זאת, את רוב כישלונות NIST אני כן יכול לייחס לסיבה מבנית. המספר הגדול ביותר האפשרי הוא 37 = 100101, ולכן הביט העליון הוא 1 רק עבור 32–37 ורצף של 7 אחדות או יותר לא יכול להיווצר. גם המיון בתוך השורה יוצר את המחזוריות של 6.
   - זה הסבר מבני, לא מבחן נוסף.
4. **דברים חריגים שמתרחשים:** כל מה שבטבלת החריגות. הבולטים ביותר הם המחזוריות של 6, היעדר מוחלט של רצפי אחדות ארוכים, ו-J=8.
5. **תכונות או רכיבים שלא אמורים להופיע:**
   - אם זו סדרה של מספרים בלתי תלויים, לא אמורים להופיע הסדר הקבוע Q<W<E<R<T<Y ושלוש התוצאות שלו: המחזוריות של 6, ההבדלים השיטתיים בין העמודות (Q בין 1 ל-23, Y בין 11 ל-37) ו-50 שורות שבהן Y=37.
   - זה מתיישב עם הנחת העבודה שנמסרה, אבל זה הופך את הרצף הממוין לבלתי מתאים למבחנים שמניחים אי-תלות.

## 5. סיכום מילולי לבודק חיצוני

- **מה נבדק:** 328 שורות, כל אחת עם 6 מספרים ראשיים בטווח 1–37 ועמודה נוספת (ABC) בטווח 1–7. כל המבחנים שהתבקשו בוצעו מול ספים של 0.01 ושל 0.05.
- **ברמת המספרים:**
  - התפלגות התדירויות מתיישבת עם אחידות (χ² p=0.153, וב-ABC p=0.232).
  - הסדרות לאורך השורות בעמודות W, E, R, T, Y ו-ABC עברו את מבחן הרצפים.
  - בעמודה Q מבחן הרצפים נכשל בשני הספים.
  - יחס האנטרופיה (0.9969 במספרים) לא הגיע לסף 0.999.
- **ברמת הביטים, לפי NIST:** 2 מבחנים עברו, 10 נכשלו ו-3 לא ישימים. הכישלונות נקבעים במידה מכרעת על ידי בחירת הקידוד (6 ביטים לטווח 1–37 יוצר חוסר איזון) ועל ידי המיון בתוך השורה. לכן הם אינם ראיה ישירה נגד היותה של סדרת ההגרלות אקראית.
- **הקריטריון "מעבר של 100%" לא התקיים.**
- **הערה אפיסטמולוגית:** גם אילו כל המבחנים היו עוברים, הדבר לא היה מוכיח שמקור הנתונים הוא הגרלה. מבחני אקראיות יכולים רק לא לדחות מודל, לא לאשר מקור.

## 6. גרסאות התוכנה

Python 3.13.18, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, matplotlib 3.11.2, openpyxl 3.1.5.

את מבחני NIST מימשתי מחדש בפייתון לפי המפרט SP 800-22 rev1a. זה לא תוכנת ה-STS הרשמית, שאינה זמינה בסביבה. מבחינת בודק חיצוני, זו הסתייגות מהותית.

## 7. הקוד המלא

**`build_csv.py`:** המשתנה `RAW` מכיל את מחרוזת הנתונים בדיוק כפי שנמסרה, והושמט כאן רק בגלל אורכו.
```python
import hashlib
RAW = """row,Y,T,W,R,ABC,Q,E,Day 1,37,29,7,25,7,2,15,שלישי 2,... (verbatim, through) 328,28,22,9,21,1,5,12,שלישי"""
data = ("\n".join(RAW.split(" ")) + "\n").encode("utf-8")
open("full_table_8_columns.csv","wb").write(data)
h = hashlib.sha256(data).hexdigest()
print("lines:", data.count(b"\n"), "sha256:", h)
print("MATCH" if h=="b60bd7b7cd35bb708573189dfccb3f0c989890fea527f8b67b573fb89bb763b4" else "MISMATCH -> STOP")
```

**`analysis.py`:**
```python
# Lotto-table audit: frequency/chi2, autocorrelation, runs up/down, Shannon entropy, NIST SP 800-22 (15 tests).
# Deterministic. No simulation, no resampling, no seeds, no multiplicity correction, no p-value adjustment.
import hashlib, json, math, itertools, sys
import numpy as np, pandas as pd, scipy, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import gammaincc as igamc, erfc
from scipy.stats import norm, chi2 as chi2d

A1, A2 = 0.01, 0.05
out = {}
def P(*a):
    print(*a); sys.stdout.flush()

raw = open("full_table_8_columns.csv","rb").read()
P("SHA256", hashlib.sha256(raw).hexdigest())
df = pd.read_csv("full_table_8_columns.csv", encoding="utf-8")
P("shape", df.shape, "columns", list(df.columns))
MAIN = ["Q","W","E","R","T","Y"]
ok = ((df["Q"]<df["W"])&(df["W"]<df["E"])&(df["E"]<df["R"])&(df["R"]<df["T"])&(df["T"]<df["Y"]))
P("rows with Q<W<E<R<T<Y:", int(ok.sum()), "/", len(df))
P("ranges:", {c:(int(df[c].min()),int(df[c].max())) for c in MAIN+["ABC"]})
P("row ids consecutive 1..328:", list(df.row)==list(range(1,329)))
P("duplicate full main-number rows:", int(df.duplicated(subset=MAIN).sum()))

def pf(p):  # pass flags vs both thresholds
    return ("PASS" if p>A1 else "FAIL", "PASS" if p>A2 else "FAIL")

# ---------------- 1. Uniformity: chi-square frequency test ----------------
P("\n=== 1. FREQUENCY TEST (chi-square vs uniform) ===")
main = df[MAIN].to_numpy().ravel()            # row-major stream Q,W,E,R,T,Y
cnt = np.bincount(main, minlength=38)[1:]
exp = len(main)/37
x2 = float(((cnt-exp)**2/exp).sum()); p_main = float(chi2d.sf(x2, 36))
P(f"main numbers 1..37 pooled: N={len(main)} expected/cell={exp:.4f} chi2={x2:.6f} df=36 p={p_main:.6f} {pf(p_main)}")
P("counts 1..37:", cnt.tolist())
abc = df.ABC.to_numpy(); ca = np.bincount(abc, minlength=8)[1:]; ea = len(abc)/7
x2a = float(((ca-ea)**2/ea).sum()); p_abc = float(chi2d.sf(x2a, 6))
P(f"ABC 1..7: N={len(abc)} expected/cell={ea:.4f} chi2={x2a:.6f} df=6 p={p_abc:.6f} {pf(p_abc)}")
P("counts 1..7:", ca.tolist())
out["chi2_main"]=(x2,p_main); out["chi2_abc"]=(x2a,p_abc)

fig,ax = plt.subplots(1,2,figsize=(13,4.2))
ax[0].bar(range(1,38),cnt,color="#4a78b5"); ax[0].axhline(exp,color="k",ls="--",lw=1,label=f"uniform expectation {exp:.2f}")
ax[0].set_title(f"Main numbers 1-37 (N={len(main)}), chi2 p={p_main:.4f}"); ax[0].set_xlabel("number"); ax[0].set_ylabel("count"); ax[0].legend()
ax[1].bar(range(1,8),ca,color="#c07a32"); ax[1].axhline(ea,color="k",ls="--",lw=1,label=f"uniform expectation {ea:.2f}")
ax[1].set_title(f"ABC 1-7 (N={len(abc)}), chi2 p={p_abc:.4f}"); ax[1].set_xlabel("value"); ax[1].legend()
plt.tight_layout(); plt.savefig("fig1_frequency_histogram.png",dpi=130); plt.close()

# ---------------- 2. Independence: autocorrelation + runs up/down ----------------
P("\n=== 2a. AUTOCORRELATION (lags 1..20; bands 1.96/sqrt(n) [0.05] and 2.5758/sqrt(n) [0.01]) ===")
def acf(x, K=20):
    x = np.asarray(x,float); x = x-x.mean(); d = (x*x).sum()
    return np.array([(x[:-k]*x[k:]).sum()/d for k in range(1,K+1)])
series = {"STREAM(row-major Q..Y)": main}
for c in MAIN+["ABC"]: series[c] = df[c].to_numpy()
acf_res = {}
for name,x in series.items():
    r = acf(x); n=len(x); b5=norm.ppf(1-A2/2)/math.sqrt(n); b1=norm.ppf(1-A1/2)/math.sqrt(n)
    acf_res[name]=(r,b5,b1)
    P(f"{name}: n={n} band05=±{b5:.4f} band01=±{b1:.4f} | lags outside band05: {[k+1 for k in range(20) if abs(r[k])>b5]} | outside band01: {[k+1 for k in range(20) if abs(r[k])>b1]}")
    P("   r =", " ".join(f"{v:+.4f}" for v in r))
fig,axs = plt.subplots(4,2,figsize=(13,12)); axs=axs.ravel()
for i,(name,(r,b5,b1)) in enumerate(acf_res.items()):
    a=axs[i]; a.bar(range(0,21),np.r_[1.0,r],color=["#999"]+["#4a78b5"]*20,width=0.6)
    for b,ls in ((b5,"--"),(b1,":")): a.axhline(b,color="r",ls=ls,lw=1); a.axhline(-b,color="r",ls=ls,lw=1)
    a.set_title(name,fontsize=10); a.set_ylim(-0.4,1.05)
plt.suptitle("Correlograms (lag 0 shown grey; dashed=95% band, dotted=99% band)"); plt.tight_layout(); plt.savefig("fig2_correlogram.png",dpi=130); plt.close()

P("\n=== 2b. RUNS UP-AND-DOWN TEST (zero differences removed; Wallis-Moore normal approx; two-sided) ===")
def runs_updown(x):
    d = np.diff(np.asarray(x,float)); s = np.sign(d); s = s[s!=0]
    n = len(s)+1; R = 1+int((s[1:]!=s[:-1]).sum())
    mu=(2*n-1)/3; var=(16*n-29)/90; z=(R-mu)/math.sqrt(var); p=2*norm.sf(abs(z))
    return n,R,mu,var,z,p
runs_res={}
for name,x in series.items():
    n,R,mu,var,z,p = runs_updown(x); runs_res[name]=(n,R,z,p)
    P(f"{name}: n_eff={n} runs={R} E={mu:.3f} Var={var:.3f} z={z:+.4f} p={p:.6f} {pf(p)}")

# ---------------- 3. Shannon entropy ----------------
P("\n=== 3. SHANNON ENTROPY (base 2; ratio to maximum; threshold ratio >= 0.999) ===")
def H(c):
    p=np.asarray(c,float); p=p[p>0]/p.sum(); return float(-(p*np.log2(p)).sum())
# bit stream: each main number written as 6-bit big-endian binary, row-major order Q,W,E,R,T,Y
bits = np.array([int(b) for v in main for b in format(int(v),"06b")],dtype=np.int8)
ent = {"main numbers (max log2 37)":(H(cnt),math.log2(37)), "ABC (max log2 7)":(H(ca),math.log2(7)),
       "bit stream (max 1)":(H(np.bincount(bits,minlength=2)),1.0)}
for k,(h,m) in ent.items():
    P(f"{k}: H={h:.6f} Hmax={m:.6f} ratio={h/m:.6f} -> {'MEETS' if h/m>=0.999 else 'DOES NOT MEET'} 0.999")
out["entropy"]=ent

# ---------------- 4. NIST SP 800-22 rev1a ----------------
n = len(bits); eps = bits.astype(int)
P(f"\n=== 4. NIST SP 800-22 (re-implemented in Python from the rev1a specification) ===")
P(f"bit conversion: each main number -> 6-bit binary; n={n} bits; ones={int(eps.sum())} zeros={n-int(eps.sum())}")
P("first 72 bits:", "".join(map(str,eps[:72])))
NI = []   # (test, p, note)
def add(t,p,note=""):
    NI.append((t,p,note)); P(f"{t}: p={'NA' if p is None else f'{p:.6f}'} {'' if p is None else pf(p)} {note}")

# 1 Frequency (monobit)
S = int((2*eps-1).sum()); add("01 Frequency (Monobit)", float(erfc(abs(S)/math.sqrt(n)/math.sqrt(2))), f"S_n={S}")
# 2 Block frequency
M=128; N=n//M; pis=eps[:N*M].reshape(N,M).mean(1); x=4*M*((pis-0.5)**2).sum()
add("02 Block Frequency", float(igamc(N/2,x/2)), f"M={M} N={N} chi2={x:.4f}")
# 3 Runs
pi=eps.mean(); tau=2/math.sqrt(n)
if abs(pi-0.5)>=tau: add("03 Runs", 0.0, f"prerequisite failed: |pi-0.5|={abs(pi-0.5):.5f} >= tau={tau:.5f}; NIST assigns p=0.0")
else:
    V=1+int((eps[1:]!=eps[:-1]).sum()); add("03 Runs", float(erfc(abs(V-2*n*pi*(1-pi))/(2*math.sqrt(2*n)*pi*(1-pi)))), f"V={V}")
# 4 Longest run of ones (n>=6272 -> M=128)
M=128; N=n//M; K=5; piL=[0.1174,0.2430,0.2493,0.1752,0.1027,0.1124]; v=np.zeros(6)
for blk in eps[:N*M].reshape(N,M):
    L=c=0
    for b in blk:
        c=c+1 if b else 0; L=max(L,c)
    v[min(max(L,4),9)-4]+=1
x=float(((v-N*np.array(piL))**2/(N*np.array(piL))).sum())
add("04 Longest Run of Ones", float(igamc(K/2,x/2)), f"M={M} N={N} v={v.astype(int).tolist()} chi2={x:.4f}")
# 5 Binary matrix rank
def rank2(mat):
    m=mat.copy(); r=0; rows,cols=m.shape
    for c in range(cols):
        piv=[i for i in range(r,rows) if m[i,c]]
        if not piv: continue
        m[[r,piv[0]]]=m[[piv[0],r]]
        for i in range(rows):
            if i!=r and m[i,c]: m[i]^=m[r]
        r+=1
        if r==rows: break
    return r
N=n//1024; ranks=[rank2(eps[i*1024:(i+1)*1024].reshape(32,32).astype(np.uint8)) for i in range(N)]
F32=ranks.count(32); F31=ranks.count(31); Fr=N-F32-F31
x=(F32-0.2888*N)**2/(0.2888*N)+(F31-0.5776*N)**2/(0.5776*N)+(Fr-0.1336*N)**2/(0.1336*N)
add("05 Binary Matrix Rank", float(math.exp(-x/2)), f"N={N} matrices (NIST requires N>=38, i.e. n>=38912) F32={F32} F31={F31} rest={Fr} chi2={x:.4f}")
# 6 DFT
X=2*eps-1; Sf=np.abs(np.fft.fft(X))[:n//2]; T=math.sqrt(math.log(1/0.05)*n); N0=0.95*n/2; N1=int((Sf<T).sum())
d=(N1-N0)/math.sqrt(n*0.95*0.05/4); add("06 Discrete Fourier Transform", float(erfc(abs(d)/math.sqrt(2))), f"N0={N0:.1f} N1={N1} d={d:.4f}")
# 7 Non-overlapping template (m=9, all aperiodic templates, N=8)
m=9; N=8; M=n//N; mu=(M-m+1)/2**m; s2=M*(1/2**m-(2*m-1)/2**(2*m))
temps=[t for t in itertools.product([0,1],repeat=m) if all(t[k:]!=t[:m-k] for k in range(1,m))]
blocks=eps[:N*M].reshape(N,M); nt_p=[]
for t in temps:
    t=np.array(t); W=[]
    for blk in blocks:
        j=w=0
        while j<=M-m:
            if np.array_equal(blk[j:j+m],t): w+=1; j+=m
            else: j+=1
        W.append(w)
    x=float(((np.array(W)-mu)**2/s2).sum()); nt_p.append((("".join(map(str,t))),float(igamc(N/2,x/2)),W))
ntv=np.array([q for _,q,_ in nt_p])
P(f"07 Non-overlapping Template: {len(temps)} aperiodic templates m=9, N={N}, M={M}")
for s,q,W in nt_p: P(f"   template {s}: W={W} p={q:.6f} {pf(q)}")
P(f"   templates p>0.01: {int((ntv>A1).sum())}/{len(ntv)} ; p>0.05: {int((ntv>A2).sum())}/{len(ntv)} ; min p={ntv.min():.6g} ; median p={np.median(ntv):.6f}")
out["nonover"]=nt_p
# 8 Overlapping template (m=9 all ones, M=1032, K=5)
m=9; M=1032; N=n//M; K=5; piO=np.array([0.364091,0.185659,0.139381,0.100571,0.070432,0.139865]); v=np.zeros(6)
for blk in eps[:N*M].reshape(N,M):
    w=sum(1 for j in range(M-m+1) if blk[j:j+m].all()); v[min(w,5)]+=1
x=float(((v-N*piO)**2/(N*piO)).sum())
add("08 Overlapping Template", float(igamc(K/2,x/2)), f"N={N} blocks (NIST spec designed for n=10^6, N=968) v={v.astype(int).tolist()} chi2={x:.4f}")
# 9 Maurer universal
add("09 Maurer's Universal Statistical", None, f"NOT COMPUTABLE: NIST minimum n=387840 (L=6, Q=640); available n={n}")
# 10 Linear complexity (M=500)
def bm(s):
    N_=len(s); c=np.zeros(N_,int); b=np.zeros(N_,int); c[0]=b[0]=1; L=0; mm=-1
    for i in range(N_):
        dd=(s[i]+int(np.dot(c[1:L+1],s[i-1::-1][:L])))%2 if L>0 else s[i]
        if dd:
            t=c.copy(); sh=i-mm; c[sh:]=(c[sh:]+b[:N_-sh])%2
            if L<=i//2: L=i+1-L; mm=i; b=t
    return L
M=500; N=n//M; K=6; muL=M/2+(9+(-1)**(M+1))/36-(M/3+2/9)/2**M
piC=np.array([0.010417,0.03125,0.125,0.5,0.25,0.0625,0.020833]); v=np.zeros(7); Ls=[]
for blk in eps[:N*M].reshape(N,M):
    L=bm(blk); Ls.append(L); Tt=(-1)**M*(L-muL)+2/9
    idx=0 if Tt<=-2.5 else 1 if Tt<=-1.5 else 2 if Tt<=-0.5 else 3 if Tt<=0.5 else 4 if Tt<=1.5 else 5 if Tt<=2.5 else 6; v[idx]+=1
x=float(((v-N*piC)**2/(N*piC)).sum())
add("10 Linear Complexity", float(igamc(K/2,x/2)), f"M={M} N={N} blocks (NIST recommends n>=10^6) L={Ls} v={v.astype(int).tolist()} chi2={x:.4f}")
# 11 Serial (m = floor(log2 n)-3 = 10, max allowed by m < floor(log2 n)-2)
def psi2(m):
    if m<=0: return 0.0
    e=np.r_[eps,eps[:m-1]]; idx=np.zeros(n,int)
    for k in range(m): idx=idx*2+e[k:k+n]
    c=np.bincount(idx,minlength=2**m); return float(2**m/n*(c.astype(float)**2).sum()-n)
m=int(math.floor(math.log2(n)))-3
p0,p1,p2=psi2(m),psi2(m-1),psi2(m-2); d1=p0-p1; d2=p0-2*p1+p2
add("11a Serial (p-value 1)", float(igamc(2**(m-2),d1/2)), f"m={m} del_psi2={d1:.4f}")
add("11b Serial (p-value 2)", float(igamc(2**(m-3),d2/2)), f"m={m} del2_psi2={d2:.4f}")
# 12 Approximate entropy (m = floor(log2 n)-6 = 7, max allowed by m < floor(log2 n)-5)
def phi(m):
    e=np.r_[eps,eps[:m-1]]; idx=np.zeros(n,int)
    for k in range(m): idx=idx*2+e[k:k+n]
    c=np.bincount(idx,minlength=2**m)/n; c=c[c>0]; return float((c*np.log(c)).sum())
m=int(math.floor(math.log2(n)))-6; ap=phi(m)-phi(m+1); x=2*n*(math.log(2)-ap)
add("12 Approximate Entropy", float(igamc(2**(m-1),x/2)), f"m={m} ApEn={ap:.6f} chi2={x:.4f}")
# 13 Cumulative sums
def cusum(X):
    z=int(np.max(np.abs(np.cumsum(X)))); s1=sum(norm.cdf((4*k+1)*z/math.sqrt(n))-norm.cdf((4*k-1)*z/math.sqrt(n)) for k in range(int((-n/z+1)//4),int((n/z-1)//4)+1))
    s2=sum(norm.cdf((4*k+3)*z/math.sqrt(n))-norm.cdf((4*k+1)*z/math.sqrt(n)) for k in range(int((-n/z-3)//4),int((n/z-1)//4)+1))
    return z,1-s1+s2
z,p=cusum(X); add("13a Cumulative Sums (forward)", float(p), f"z={z}")
z,p=cusum(X[::-1]); add("13b Cumulative Sums (backward)", float(p), f"z={z}")
# 14/15 Random excursions (+variant)
Sw=np.r_[0,np.cumsum(X),0]; zeros=np.where(Sw==0)[0]; J=len(zeros)-1
P(f"Random-excursion walk: J={J} cycles; final partial sum S_n={int(np.cumsum(X)[-1])}; NIST rule: if J<500 the test is not applicable")
add("14 Random Excursions", None, f"NOT APPLICABLE: J={J} < 500 (NIST)")
add("15 Random Excursions Variant", None, f"NOT APPLICABLE: J={J} < 500 (NIST)")

# ---------------- Scorecard ----------------
rows=[]
for t,p,note in NI:
    if p is None: rows.append([t,"NA","NA","NA","NA","NA",note]); continue
    rows.append([t,f"{p:.6f}","1/1" if p>A1 else "0/1","1/1" if p>A2 else "0/1",pf(p)[0],pf(p)[1],note])
nt_pass1=int((ntv>A1).sum()); nt_pass5=int((ntv>A2).sum())
rows.insert(6,["07 Non-overlapping Template (148 templates)",f"min={ntv.min():.6g}; median={np.median(ntv):.6f}",
   f"{nt_pass1}/148",f"{nt_pass5}/148","PASS" if nt_pass1==148 else "FAIL","PASS" if nt_pass5==148 else "FAIL","per-template p-values in machine output"])
sc=pd.DataFrame(rows,columns=["test","p_value","proportion@0.01","proportion@0.05","status@0.01","status@0.05","note"])
sc.to_csv("nist_scorecard.csv",index=False,encoding="utf-8")
P("\n=== SCORECARD ==="); P(sc.to_string(index=False))
fam = {}
for t,p,_ in NI:
    k=t[:2]; fam.setdefault(k,[]).append(p)
fam["07"]=list(ntv)
def famstat(ps,a):
    if any(p is None for p in ps): return "NA"
    return "PASS" if all(p>a for p in ps) else "FAIL"
s1=[famstat(v,A1) for v in fam.values()]; s5=[famstat(v,A2) for v in fam.values()]
P(f"\n15-test summary @0.01: PASS={s1.count('PASS')} FAIL={s1.count('FAIL')} NA={s1.count('NA')} of 15")
P(f"15-test summary @0.05: PASS={s5.count('PASS')} FAIL={s5.count('FAIL')} NA={s5.count('NA')} of 15")

# ---------------- mandated descriptive tables ----------------
P("\n=== DAY column ===")
P(df.Day.value_counts().to_string())
tr=pd.crosstab(df.Day[:-1].reset_index(drop=True),df.Day[1:].reset_index(drop=True),rownames=["from"],colnames=["to"])
tr.to_csv("table_transitions_day.csv",encoding="utf-8"); P("Day->next Day transition counts:"); P(tr.to_string())
same=[int(df.row[i+1]) for i in range(len(df)-1) if df.Day[i]==df.Day[i+1]]
P("rows whose Day equals previous row's Day:", same)
P("rows with Y=37:", int((df.Y==37).sum()), " rows with Q=1:", int((df.Q==1).sum()))
P("consecutive-number pairs inside a row (rows count):", int(sum(any(b-a==1 for a,b in zip(r[:-1],r[1:])) for r in df[MAIN].to_numpy())))
P("\nversions: python", sys.version.split()[0], "numpy", np.__version__, "scipy", scipy.__version__, "pandas", pd.__version__, "matplotlib", matplotlib.__version__)
with pd.ExcelWriter("report.xlsx") as w:
    df.to_excel(w,sheet_name="full_table",index=False); sc.to_excel(w,sheet_name="NIST_scorecard",index=False)
    pd.DataFrame([[t,q,str(W)] for t,q,W in nt_p],columns=["template","p","W"]).to_excel(w,sheet_name="nonoverlap_templates",index=False)
    tr.to_excel(w,sheet_name="day_transitions")
```

## 8. פלט המכונה המלא (`machine_output.txt`)

הפלט מוצג כפי שהוא. היוצא מן הכלל הוא 148 שורות התבניות, שמוצגות בפורמט דחוס: `template | W לכל 8 הבלוקים | p | @0.01 | @0.05`. הערכים זהים לפלט המקורי, והשורות המלאות נמצאות בקובץ.

```
SHA256 b60bd7b7cd35bb708573189dfccb3f0c989890fea527f8b67b573fb89bb763b4
shape (328, 9) columns ['row', 'Y', 'T', 'W', 'R', 'ABC', 'Q', 'E', 'Day']
rows with Q<W<E<R<T<Y: 328 / 328
ranges: {'Q': (1, 23), 'W': (2, 28), 'E': (3, 33), 'R': (6, 35), 'T': (7, 36), 'Y': (11, 37), 'ABC': (1, 7)}
row ids consecutive 1..328: True
duplicate full main-number rows: 0

=== 1. FREQUENCY TEST (chi-square vs uniform) ===
main numbers 1..37 pooled: N=1968 expected/cell=53.1892 chi2=44.627033 df=36 p=0.153282 ('PASS', 'PASS')
counts 1..37: [67, 57, 50, 51, 47, 46, 58, 49, 48, 54, 59, 39, 46, 47, 42, 45, 46, 64, 50, 48, 72, 52, 53, 66, 61, 57, 67, 55, 52, 54, 45, 48, 47, 63, 67, 46, 50]
ABC 1..7: N=328 expected/cell=46.8571 chi2=8.085366 df=6 p=0.231916 ('PASS', 'PASS')
counts 1..7: [38, 53, 39, 45, 60, 50, 43]

=== 2a. AUTOCORRELATION (lags 1..20; bands 1.96/sqrt(n) [0.05] and 2.5758/sqrt(n) [0.01]) ===
STREAM(row-major Q..Y): n=1968 band05=±0.0442 band01=±0.0581 | lags outside band05: [1..20 all] | outside band01: [1..20 all]
   r = +0.2764 -0.1695 -0.3292 -0.2477 +0.1118 +0.7285 +0.0993 -0.2773 -0.3891 -0.2741 +0.1015 +0.7275 +0.1046 -0.2719 -0.3879 -0.2776 +0.0912 +0.7173 +0.1019 -0.2656
Q: n=328 band05=±0.1082 band01=±0.1422 | lags outside band05: [8] | outside band01: []
   r = -0.0726 +0.0389 +0.0223 -0.0107 +0.0665 +0.0575 +0.0470 +0.1148 -0.0141 -0.0781 +0.0471 +0.0258 -0.0687 +0.0704 -0.0257 +0.0331 -0.0176 +0.0438 -0.0799 +0.0418
W: n=328 band05=±0.1082 band01=±0.1422 | lags outside band05: [] | outside band01: []
   r = -0.0894 +0.0505 -0.0119 +0.0076 -0.0293 -0.0337 -0.0609 +0.0363 +0.0412 +0.0848 -0.0305 +0.0205 +0.0583 +0.0949 -0.0642 +0.0171 +0.0489 +0.0033 -0.0470 +0.0221
E: n=328 band05=±0.1082 band01=±0.1422 | lags outside band05: [10, 14, 16] | outside band01: []
   r = +0.0311 -0.0583 -0.0119 -0.0025 -0.0453 -0.0226 -0.0818 +0.0321 +0.0595 +0.1411 -0.0093 +0.0636 -0.0360 +0.1176 -0.0432 -0.1179 -0.0298 +0.0456 +0.0263 +0.0138
R: n=328 band05=±0.1082 band01=±0.1422 | lags outside band05: [7, 14] | outside band01: []
   r = +0.0615 +0.0062 -0.0383 +0.0620 -0.0207 -0.0132 -0.1209 -0.0096 +0.0585 +0.0659 -0.0588 +0.0158 -0.0113 +0.1375 -0.0208 -0.0637 -0.0697 -0.0711 -0.0102 +0.0107
T: n=328 band05=±0.1082 band01=±0.1422 | lags outside band05: [11, 16] | outside band01: []
   r = +0.0556 +0.0282 -0.0500 +0.0745 -0.0012 +0.0849 -0.0816 +0.0256 -0.0106 -0.0017 -0.1098 -0.0064 +0.0265 -0.0181 -0.0697 -0.1101 -0.0215 -0.0327 +0.0221 -0.0031
Y: n=328 band05=±0.1082 band01=±0.1422 | lags outside band05: [6, 11] | outside band01: [11]
   r = +0.0768 -0.0330 -0.0791 +0.0411 +0.0055 +0.1315 +0.0292 -0.0134 +0.0934 +0.0493 -0.1610 -0.0401 -0.0098 -0.0091 -0.0281 -0.0049 -0.0537 -0.0679 +0.0454 -0.0163
ABC: n=328 band05=±0.1082 band01=±0.1422 | lags outside band05: [] | outside band01: []
   r = +0.0859 -0.0192 +0.0008 -0.0058 -0.0380 -0.0005 -0.1034 -0.0299 +0.0209 +0.0469 -0.0526 -0.0243 -0.0556 -0.0437 +0.0217 +0.0663 +0.0290 -0.0300 +0.0796 -0.0175

=== 2b. RUNS UP-AND-DOWN TEST (zero differences removed; Wallis-Moore normal approx; two-sided) ===
STREAM(row-major Q..Y): n_eff=1968 runs=653 E=1311.667 Var=349.544 z=-35.2301 p=0.000000 ('FAIL', 'FAIL')
Q: n_eff=286 runs=209 E=190.333 Var=50.522 z=+2.6262 p=0.008635 ('FAIL', 'FAIL')
W: n_eff=310 runs=219 E=206.333 Var=54.789 z=+1.7113 p=0.087033 ('PASS', 'PASS')
E: n_eff=316 runs=205 E=210.333 Var=55.856 z=-0.7136 p=0.475464 ('PASS', 'PASS')
R: n_eff=317 runs=218 E=211.000 Var=56.033 z=+0.9351 p=0.349718 ('PASS', 'PASS')
T: n_eff=310 runs=208 E=206.333 Var=54.789 z=+0.2252 p=0.821850 ('PASS', 'PASS')
Y: n_eff=285 runs=190 E=189.667 Var=50.344 z=+0.0470 p=0.962530 ('PASS', 'PASS')
ABC: n_eff=283 runs=191 E=188.333 Var=49.989 z=+0.3772 p=0.706051 ('PASS', 'PASS')

=== 3. SHANNON ENTROPY (base 2; ratio to maximum; threshold ratio >= 0.999) ===
main numbers (max log2 37): H=5.193458 Hmax=5.209453 ratio=0.996930 -> DOES NOT MEET 0.999
ABC (max log2 7): H=2.789867 Hmax=2.807355 ratio=0.993771 -> DOES NOT MEET 0.999
bit stream (max 1): H=0.981636 Hmax=1.000000 ratio=0.981636 -> DOES NOT MEET 0.999

=== 4. NIST SP 800-22 (re-implemented in Python from the rev1a specification) ===
bit conversion: each main number -> 6-bit binary; n=11808 bits; ones=4964 zeros=6844
first 72 bits: 000010000111001111011001011101100101010010010110011000011001011011011100
01 Frequency (Monobit): p=0.000000 ('FAIL', 'FAIL') S_n=-1880
02 Block Frequency: p=0.000000 ('FAIL', 'FAIL') M=128 N=92 chi2=359.3750
03 Runs: p=0.000000 ('FAIL', 'FAIL') prerequisite failed: |pi-0.5|=0.07961 >= tau=0.01841; NIST assigns p=0.0
04 Longest Run of Ones: p=0.000000 ('FAIL', 'FAIL') M=128 N=92 v=[51, 15, 26, 0, 0, 0] chi2=188.3537
05 Binary Matrix Rank: p=0.893303 ('PASS', 'PASS') N=11 matrices (NIST requires N>=38, i.e. n>=38912) F32=3 F31=7 rest=1 chi2=0.2257
06 Discrete Fourier Transform: p=0.000000 ('FAIL', 'FAIL') N0=5608.8 N1=5716 d=9.0529
07 Non-overlapping Template: 148 aperiodic templates m=9, N=8, M=1476
000000001|2,1,2,5,0,0,2,2|0.274883|PASS|PASS
000000011|1,1,1,3,0,0,1,1|0.144206|PASS|PASS
000000101|0,0,2,1,1,2,1,2|0.234126|PASS|PASS
000000111|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
000001001|11,5,4,1,6,4,4,10|0.000000|FAIL|FAIL
000001011|2,1,2,0,1,3,0,2|0.324341|PASS|PASS
000001101|2,1,0,2,3,0,0,2|0.206537|PASS|PASS
000001111|0,1,0,0,1,0,0,1|0.017771|PASS|FAIL
000010001|10,10,14,11,13,9,12,5|0.000000|FAIL|FAIL
000010011|1,1,2,0,1,0,1,1|0.133285|PASS|PASS
000010101|0,2,1,1,3,3,0,1|0.269663|PASS|PASS
000010111|1,0,1,0,0,1,1,0|0.032239|PASS|FAIL
000011001|3,6,4,7,5,5,6,8|0.000913|FAIL|FAIL
000011011|1,3,0,3,2,2,0,4|0.417731|PASS|PASS
000011101|3,1,1,1,0,0,1,3|0.206537|PASS|PASS
000011111|2,1,2,0,0,2,0,1|0.144206|PASS|PASS
000100011|5,3,4,3,4,4,5,1|0.657413|PASS|PASS
000100101|7,3,4,2,3,4,3,6|0.209795|PASS|PASS
000100111|0,0,1,0,1,0,1,1|0.032239|PASS|FAIL
000101001|4,10,11,6,9,6,10,6|0.000000|FAIL|FAIL
000101011|2,0,1,3,1,2,0,0|0.155881|PASS|PASS
000101101|4,4,4,5,6,0,3,2|0.282235|PASS|PASS
000101111|1,0,0,0,0,1,0,0|0.009608|FAIL|FAIL
000110011|3,2,2,4,1,1,1,2|0.754479|PASS|PASS
000110101|3,3,3,6,4,2,4,4|0.736924|PASS|PASS
000110111|1,0,0,1,0,3,0,3|0.073968|PASS|PASS
000111001|12,12,13,11,7,13,13,5|0.000000|FAIL|FAIL
000111011|3,0,0,1,1,3,2,0|0.168343|PASS|PASS
000111101|3,3,2,1,5,6,2,1|0.413454|PASS|PASS
000111111|3,1,0,0,1,1,0,1|0.085524|PASS|PASS
001000011|6,9,5,4,4,8,10,4|0.000000|FAIL|FAIL
001000101|2,9,14,9,11,3,8,5|0.000000|FAIL|FAIL
001000111|4,3,3,2,2,1,1,2|0.876036|PASS|PASS
001001011|3,2,7,3,1,5,6,6|0.037808|PASS|FAIL
001001101|4,4,4,3,2,1,8,4|0.118139|PASS|PASS
001001111|0,1,0,0,2,2,0,2|0.085524|PASS|PASS
001010011|2,1,4,3,3,4,3,4|0.939217|PASS|PASS
001010101|3,3,6,3,8,8,4,9|0.000015|FAIL|FAIL
001010111|0,0,0,1,0,1,0,0|0.009608|FAIL|FAIL
001011011|2,3,3,1,1,2,3,2|0.911957|PASS|PASS
001011101|3,5,5,4,3,1,6,5|0.255215|PASS|PASS
001011111|0,0,1,2,0,0,0,0|0.013735|PASS|FAIL
001100101|2,2,1,6,4,3,4,8|0.047008|PASS|FAIL
001100111|1,0,1,3,1,1,0,0|0.085524|PASS|PASS
001101011|2,2,3,3,2,3,0,4|0.834993|PASS|PASS
001101101|2,3,2,3,3,1,7,4|0.395565|PASS|PASS
001101111|0,0,1,0,0,0,2,1|0.025132|PASS|FAIL
001110101|3,2,4,1,1,0,0,1|0.238812|PASS|PASS
001110111|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
001111011|2,1,1,1,2,1,2,1|0.529745|PASS|PASS
001111101|0,0,1,0,0,0,0,0|0.005108|FAIL|FAIL
001111111|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
010000011|8,7,3,10,10,3,2,8|0.000000|FAIL|FAIL
010000111|3,3,2,1,0,3,3,0|0.489446|PASS|PASS
010001011|4,8,8,8,7,2,7,3|0.000002|FAIL|FAIL
010001111|3,0,0,0,2,0,0,0|0.021425|PASS|FAIL
010010011|6,6,5,7,5,5,5,5|0.006316|FAIL|FAIL
010010111|0,2,1,1,0,0,3,1|0.116131|PASS|PASS
010011011|3,3,5,4,2,4,11,4|0.000699|FAIL|FAIL
010011111|0,1,2,0,2,1,0,0|0.062337|PASS|PASS
010100011|16,9,8,9,9,12,8,7|0.000000|FAIL|FAIL
010100111|1,0,2,1,2,0,1,2|0.234126|PASS|PASS
010101011|6,5,9,5,12,10,5,11|0.000000|FAIL|FAIL
010101111|1,0,1,0,2,0,2,2|0.144206|PASS|PASS
010110011|3,6,8,6,11,9,10,9|0.000000|FAIL|FAIL
010110111|0,1,0,0,0,0,1,2|0.025132|PASS|FAIL
010111011|3,7,10,6,7,5,5,6|0.000002|FAIL|FAIL
010111111|1,1,1,1,0,0,1,1|0.098666|PASS|PASS
011000111|2,4,2,3,2,7,3,2|0.465269|PASS|PASS
011001111|1,1,2,5,3,3,2,0|0.469828|PASS|PASS
011010111|1,0,1,2,2,2,0,1|0.234126|PASS|PASS
011011111|1,1,1,0,0,1,1,1|0.098666|PASS|PASS
011101111|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
011111111|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
100000000|2,1,2,5,0,0,2,2|0.274883|PASS|PASS
100010000|11,6,3,9,9,11,10,11|0.000000|FAIL|FAIL
100100000|10,12,14,7,5,6,8,9|0.000000|FAIL|FAIL
100101000|14,7,8,8,12,9,6,6|0.000000|FAIL|FAIL
100110000|4,4,3,3,3,5,3,1|0.872195|PASS|PASS
100111000|3,1,6,1,4,3,2,2|0.532759|PASS|PASS
101000000|1,0,1,0,0,1,2,1|0.078636|PASS|PASS
101000100|2,2,0,1,2,4,3,3|0.704919|PASS|PASS
101001000|4,7,5,6,3,2,5,4|0.078918|PASS|PASS
101001100|0,0,1,0,1,1,1,0|0.032239|PASS|FAIL
101010000|1,3,2,0,1,1,2,0|0.251413|PASS|PASS
101010100|1,0,2,0,0,1,0,2|0.062337|PASS|PASS
101011000|3,5,2,1,4,0,1,7|0.083120|PASS|PASS
101011100|2,1,0,1,1,0,0,0|0.045017|PASS|FAIL
101100000|1,6,0,1,2,0,1,1|0.065375|PASS|PASS
101100100|1,3,1,0,1,3,0,1|0.206537|PASS|PASS
101101000|4,1,1,7,4,2,2,4|0.228458|PASS|PASS
101101100|1,1,2,0,1,1,1,0|0.133285|PASS|PASS
101110000|5,0,1,2,1,3,4,5|0.305593|PASS|PASS
101110100|2,0,1,0,0,0,1,1|0.045017|PASS|FAIL
101111000|2,1,0,1,1,3,5,1|0.274883|PASS|PASS
101111100|2,0,1,1,1,0,0,0|0.045017|PASS|FAIL
110000000|1,3,0,4,1,2,0,2|0.309104|PASS|PASS
110000010|9,5,1,4,5,4,4,5|0.007085|FAIL|FAIL
110000100|0,2,3,5,5,2,3,1|0.431808|PASS|PASS
110001000|8,2,7,7,10,5,7,8|0.000000|FAIL|FAIL
110001010|5,1,2,1,3,3,4,4|0.720994|PASS|PASS
110010000|4,13,12,4,5,9,7,5|0.000000|FAIL|FAIL
110010010|1,6,4,4,4,10,2,2|0.001587|FAIL|FAIL
110010100|6,2,3,3,4,6,8,5|0.015567|PASS|FAIL
110011000|5,4,2,4,4,4,4,2|0.812037|PASS|PASS
110011010|1,4,3,1,0,1,4,5|0.320710|PASS|PASS
110100000|3,9,3,3,4,3,3,3|0.081942|PASS|PASS
110100010|1,0,1,5,4,6,3,2|0.183089|PASS|PASS
110100100|3,2,3,4,1,0,2,3|0.733874|PASS|PASS
110101000|6,7,4,2,3,3,5,5|0.091159|PASS|PASS
110101010|4,2,5,2,4,3,1,4|0.777647|PASS|PASS
110101100|2,2,3,1,1,2,0,2|0.587306|PASS|PASS
110110000|2,3,3,4,4,4,6,4|0.686859|PASS|PASS
110110010|4,4,6,2,2,2,4,2|0.649250|PASS|PASS
110110100|0,1,0,3,2,2,1,2|0.324341|PASS|PASS
110111000|5,2,2,2,3,4,5,5|0.627905|PASS|PASS
110111010|2,0,1,1,1,0,0,0|0.045017|PASS|FAIL
110111100|0,1,0,0,0,1,1,0|0.017771|PASS|FAIL
111000000|1,0,1,4,1,2,3,2|0.462430|PASS|PASS
111000010|0,1,2,1,1,0,3,1|0.191765|PASS|PASS
111000100|3,1,1,2,2,3,1,1|0.696837|PASS|PASS
111000110|3,3,1,2,0,2,2,1|0.616624|PASS|PASS
111001000|2,7,12,6,5,5,5,8|0.000000|FAIL|FAIL
111001010|7,1,1,2,5,5,4,2|0.115534|PASS|PASS
111001100|2,3,1,2,2,0,1,0|0.324341|PASS|PASS
111010000|4,8,4,2,3,3,5,4|0.120795|PASS|PASS
111010010|2,2,4,5,0,1,3,3|0.553407|PASS|PASS
111010100|1,1,1,0,0,1,0,1|0.057168|PASS|PASS
111010110|2,0,0,0,0,0,0,0|0.007371|FAIL|FAIL
111011000|3,2,4,4,6,5,6,4|0.242467|PASS|PASS
111011010|1,1,1,1,3,3,2,1|0.587306|PASS|PASS
111011100|0,0,1,1,0,0,1,0|0.017771|PASS|FAIL
111100000|2,2,0,3,1,2,2,2|0.696837|PASS|PASS
111100010|3,0,1,2,3,3,2,1|0.646110|PASS|PASS
111100100|0,1,2,0,0,3,0,1|0.067927|PASS|PASS
111100110|0,0,0,0,1,0,0,0|0.005108|FAIL|FAIL
111101000|2,6,4,3,6,6,7,4|0.022070|PASS|FAIL
111101010|2,0,2,0,2,1,0,1|0.144206|PASS|PASS
111101100|1,0,1,1,0,1,1,1|0.098666|PASS|PASS
111101110|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
111110000|2,3,2,4,2,2,2,1|0.930098|PASS|PASS
111110010|2,1,2,1,1,2,0,0|0.234126|PASS|PASS
111110100|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
111110110|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
111111000|4,5,3,4,1,3,3,3|0.872195|PASS|PASS
111111010|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
111111100|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
111111110|0,0,0,0,0,0,0,0|0.002677|FAIL|FAIL
   templates p>0.01: 108/148 ; p>0.05: 90/148 ; min p=2.58825e-39 ; median p=0.098666
08 Overlapping Template: p=0.001755 ('FAIL', 'FAIL') N=11 blocks (NIST spec designed for n=10^6, N=968) v=[11, 0, 0, 0, 0, 0] chi2=19.2122
09 Maurer's Universal Statistical: p=NA  NOT COMPUTABLE: NIST minimum n=387840 (L=6, Q=640); available n=11808
10 Linear Complexity: p=0.555031 ('PASS', 'PASS') M=500 N=23 blocks (NIST recommends n>=10^6) L=[250, 249, 250, 251, 251, 250, 252, 250, 251, 250, 249, 250, 250, 251, 251, 250, 251, 250, 250, 250, 251, 251, 247] v=[1, 0, 2, 11, 8, 1, 0] chi2=4.9129
11a Serial (p-value 1): p=0.000000 ('FAIL', 'FAIL') m=10 del_psi2=4464.4770
11b Serial (p-value 2): p=0.000000 ('FAIL', 'FAIL') m=10 del2_psi2=907.0569
12 Approximate Entropy: p=0.000000 ('FAIL', 'FAIL') m=7 ApEn=0.597706 chi2=2253.9479
13a Cumulative Sums (forward): p=0.000000 ('FAIL', 'FAIL') z=1881
13b Cumulative Sums (backward): p=0.000000 ('FAIL', 'FAIL') z=1881
Random-excursion walk: J=8 cycles; final partial sum S_n=-1880; NIST rule: if J<500 the test is not applicable
14 Random Excursions: p=NA  NOT APPLICABLE: J=8 < 500 (NIST)
15 Random Excursions Variant: p=NA  NOT APPLICABLE: J=8 < 500 (NIST)

=== SCORECARD ===  (same content as table in section 2; wide notes column in nist_scorecard.csv)

15-test summary @0.01: PASS=2 FAIL=10 NA=3 of 15
15-test summary @0.05: PASS=2 FAIL=10 NA=3 of 15

=== DAY column ===
שבת 137 | שלישי 136 | חמישי 45 | ראשון 7 | רביעי 2 | שני 1
Day->next Day transition counts: (table 3.2)
rows whose Day equals previous row's Day: [153]
rows with Y=37: 50  rows with Q=1: 67
consecutive-number pairs inside a row (rows count): 205

versions: python 3.13.18 numpy 2.5.3 scipy 1.18.1 pandas 3.0.5 matplotlib 3.11.2
```

בקטע הזה יש שלושה מקומות שקיצרתי, וכל השאר זהה לקובץ:
- **השורה של ה-STREAM:** הרשימות `[1, 2, …, 20]` נכתבו כאן כ-`[1..20 all]`.
- **ה-SCORECARD וטבלת המעברים:** הופנו לטבלאות שלמעלה, ששורותיהן זהות.
- **ספירת הימים:** הודפסה כאן בשורה אחת במקום שורה לכל יום.

## קבצים

כל הקבצים בתיקייה `/home/user/lotto_audit/`:
- **נתונים:** `full_table_8_columns.csv`
- **קוד:** `build_csv.py`, `analysis.py`
- **פלט מלא:** `machine_output.txt`
- **טבלאות:** `nist_scorecard.csv`, `table_transitions_day.csv`
- **גרפים:** `fig1_frequency_histogram.png`, `fig2_correlogram.png`
- **דו"ח אקסל:** `report.xlsx` (גיליונות: הטבלה המלאה, Scorecard, 148 התבניות ומעברי הימים)