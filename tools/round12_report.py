"""Round 12: put the report together from its parts in runs/round12_logs/.

    python tools/round12_report.py [--commit=<hash>] [--tests="70 passed, 11 skipped"]
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "runs", "round12_logs")
OUT = os.path.join(HERE, "..", "..", "..", "Vol_III", "round12_2026-11-17.md")


def rd(n):
    return open(os.path.join(R, n), encoding="utf-8").read().strip()


def main():
    arg = lambda k, d: next((a[len(k) + 3:] for a in sys.argv if a.startswith(f"--{k}=")), d)
    a = rd("items1to4_report.md")
    a = re.sub(r"^# .*\n", "## Items 1 to 4: the 6 m stator at 2,500 Hz (the agent's report)\n", a, count=1)
    a = a.replace("\n## Look at first", "\n### Look at first, items 1 to 4", 1).replace("\n## What was not run", "\n## What was not run, items 1 to 4", 1)
    t = "\n\n".join([rd("report_head.md"), a, rd("report_item5.md"), rd("report_tail.md")]) + "\n"
    t = t.replace("@@COMMIT@@", arg("commit", "(not committed yet)")).replace("@@TESTS@@", arg("tests", "not run"))
    open(OUT, "w", encoding="utf-8", newline="\n").write(t)
    print(f"wrote {OUT}: {len(t.splitlines())} lines; m-dashes {t.count(chr(8212))}, arrows {t.count(chr(8594))}, unfilled {t.count('@@')}")
    for ln in t.splitlines():
        if ln.startswith("## "):
            print("  ", ln)


if __name__ == "__main__":
    main()
