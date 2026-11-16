"""Round 9: put the report together from its parts in runs/round9_logs/ (the parts are written by hand from the logs beside them).

    python tools/round9_report.py [--commit=<hash>] [--out=<file>]
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "runs", "round9_logs")
OUT = os.path.join(HERE, "..", "..", "..", "Vol_III", "round9_reruns_2026-11-16.md")


def rd(n):
    return open(os.path.join(R, n), encoding="utf-8").read().strip()


def demote(t):
    return re.sub(r"^(#+) ", lambda m: "#" * (len(m.group(1)) + 1) + " ", t, flags=re.M)


def main():
    commit = next((a[9:] for a in sys.argv if a.startswith("--commit=")), "(not committed yet)")
    out = next((a[6:] for a in sys.argv if a.startswith("--out=")), OUT)
    A, B, C, D, E = (rd(f"report_part_{x}.md") for x in "ABCDE")
    b1, b34 = B.split("## 3. The dynamics cases"); b34 = "## 3. The dynamics cases" + b34
    head, item2 = A.split("@@PART1@@")
    c5, a3 = C.split("## A3."); a3 = "## A3." + a3
    d_app, d_a1 = D.split("## A1."); d_a1 = "## A1." + d_a1
    e_a2, e_rest = E.split("## What was not run"); e_rest = "## What was not run" + e_rest
    f5b = demote(rd("item5b_report.md")).replace("## Item 5b", "## 5b. The stator pulse, the stray field, the coil ends", 1)
    f5c = demote(rd("item5c_report.md")).replace("### Item 5c.", "## 5c.", 1)
    f5c = re.sub(r"`[cC]:[^`]*round9_lines_array\.py`", "`tools/round9_lines_array.py`", f5c)
    parts = [head, b1, item2, b34, c5, f5b, d_app, f5c, d_a1, e_a2, a3, e_rest]
    t = "\n\n".join(p.strip() for p in parts) + "\n"
    t = t.replace("@@COMMIT@@", commit)
    open(out, "w", encoding="utf-8", newline="\n").write(t)
    print(f"wrote {out}: {len(t.splitlines())} lines; m-dashes {t.count(chr(8212))}, arrows {t.count(chr(8594))}")
    for ln in t.splitlines():
        if ln.startswith("## ") or ln.startswith("# "):
            print("  ", ln)


if __name__ == "__main__":
    main()
