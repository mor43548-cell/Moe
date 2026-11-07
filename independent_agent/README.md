# Blind independent replication

- Separate cloud session: `session_01BnjUbMJzsah5Ls2o1vvy5q`. Fresh container, **no git repository** (tag `config:no-git-repo`),
  no access to this repository, to my code or to my results.
- What it received (and only that):
  1. The full merged table (`full_table_7_files.csv`, SHA-256 `4000ca5c068b37957c41e28ef48dc4d1e219b564408e3a575eae2d232d6ad697`),
     with the original column names Y,T,W,R,ABC,Q,E and no interpretation.
  2. The working assumption: "in all 328 rows Q < W < E < R < T < Y".
  3. A ban on any simulation, random or dummy data, shuffling/resampling, and seeds.
  4. The requester's message, word for word.
  5. A technical note: rebuild the CSV from the text and verify the hash; put all outputs in the final message.
- Not given: the lotto interpretation, the column relabelling, my tests or my results.

## Result
The independent session finished (about 7 minutes, no errors). Its full final message (report, machine output, code) was read
from its transcript with list_events and is summarized in the comparison in the chat. Its own files stayed in its container
(/home/user/lotto_report/), and I have no access to them. The verbatim report remains in the transcript of session_01BnjUbMJzsah5Ls2o1vvy5q and has not been copied here.
