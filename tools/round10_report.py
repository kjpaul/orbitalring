"""Round 10: put the report together from its parts in runs/round10_logs/ (the parts are written from the logs beside them).

    python tools/round10_report.py [--commit=<hash>] [--tests="70 passed, 11 skipped"] [--out=<file>]
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "runs", "round10_logs")
OUT = os.path.join(HERE, "..", "..", "..", "Vol_III", "round10_2026-11-17.md")


def rd(n):
    return open(os.path.join(R, n), encoding="utf-8").read().strip()


def demote(t):
    return re.sub(r"^(#+) ", lambda m: "#" * (len(m.group(1)) + 1) + " ", t, flags=re.M)


def main():
    arg = lambda k, d: next((a[len(k) + 3:] for a in sys.argv if a.startswith(f"--{k}=")), d)
    commit = arg("commit", "(not committed yet)"); tests = arg("tests", "not run"); out = arg("out", OUT)
    i2 = demote(rd("item2_report.md")).replace("## Round 10, item 2: the tailored winding for a local load", "## 2. The tailored winding for a local load", 1)
    i3 = demote(rd("item3_report.md")).replace("## Item 3: the mass driver's stator", "## 3. The mass driver's stator", 1)
    t = "\n\n".join([rd("report_item1.md"), i2, i3, rd("report_tail.md")]) + "\n"
    t = t.replace("@@COMMIT@@", commit).replace("@@TESTS@@", tests)
    open(out, "w", encoding="utf-8", newline="\n").write(t)
    print(f"wrote {out}: {len(t.splitlines())} lines; m-dashes {t.count(chr(8212))}, arrows {t.count(chr(8594))}")
    for ln in t.splitlines():
        if ln.startswith("## ") or ln.startswith("# "):
            print("  ", ln)


if __name__ == "__main__":
    main()
