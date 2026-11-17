"""Round 13: put the report together from its parts in runs/round13_logs/.

    python tools/round13_report.py [--commit=<hash>] [--tests="72 passed, 11 skipped"]
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "runs", "round13_logs")
OUT = os.path.join(HERE, "..", "..", "..", "Vol_III", "round13_2026-11-17.md")
NL = chr(10)


def rd(n):
    return open(os.path.join(R, n), encoding="utf-8").read().strip()


def main():
    arg = lambda k, d: next((a[len(k) + 3:] for a in sys.argv if a.startswith(f"--{k}=")), d)
    parts = [rd("report_head.md")]
    for f, title in (("items1and2_report.md", "Items 1 and 2: the ruled stator as the default, and the code appendix"), ("items3to5_report.md", "Items 3 to 5: the figures and the spreader arms")):
        a = rd(f)
        a = "## " + title + " (the agent's report)" + a[a.index(NL):]
        a = a.replace(NL + "## Look at first", NL + "### Look at first, in this part", 1).replace(NL + "## What was not run", NL + "### What was not run, in this part", 1)
        parts.append(a)
    t = (NL + NL).join(parts + [rd("report_items6and7.md"), rd("report_tail.md")]) + NL
    t = t.replace("@@COMMIT@@", arg("commit", "(not committed yet)")).replace("@@TESTS@@", arg("tests", "not run"))
    open(OUT, "w", encoding="utf-8", newline=NL).write(t)
    print(f"wrote {OUT}: {len(t.splitlines())} lines; m-dashes {t.count(chr(8212))}, arrows {t.count(chr(8594))}, unfilled {t.count('@@')}")
    for ln in t.splitlines():
        if ln.startswith("## "):
            print("  ", ln)


if __name__ == "__main__":
    main()
