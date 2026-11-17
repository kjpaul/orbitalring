"""Round 11: put the report together from its parts in runs/round11_logs/ (the parts are written from the logs beside them).

    python tools/round11_report.py [--commit=<hash>] [--tests="70 passed, 11 skipped"] [--out=<file>]
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "runs", "round11_logs")
OUT = os.path.join(HERE, "..", "..", "..", "Vol_III", "round11_2026-11-17.md")


def rd(n):
    return open(os.path.join(R, n), encoding="utf-8").read().strip()


def demote(t):
    return re.sub(r"^(#+) ", lambda m: "#" * (len(m.group(1)) + 1) + " ", t, flags=re.M)


def retitle(t, new):
    return re.sub(r"^## .*$", new, demote(t), count=1, flags=re.M)


def main():
    arg = lambda k, d: next((a[len(k) + 3:] for a in sys.argv if a.startswith(f"--{k}=")), d)
    commit = arg("commit", "(not committed yet)"); tests = arg("tests", "not run"); out = arg("out", OUT)
    i1b = retitle(rd("item1b_report.md"), "## 1b. One plate under the cable, and a wider face")
    i2 = retitle(rd("item2_report.md"), "## 2. A stator for the mass driver that carries a stated launch")
    t = "\n\n".join([rd("report_head.md"), rd("report_item1.md"), i1b, i2, rd("report_tail.md")]) + "\n"
    for k in ("cases114", "look114", "script1b", "changed1b"):
        t = t.replace("@@" + k.upper() + "@@", rd("fill_" + k + ".md"))
    t = t.replace("@@COMMIT@@", commit).replace("@@TESTS@@", tests)
    open(out, "w", encoding="utf-8", newline="\n").write(t)
    print(f"wrote {out}: {len(t.splitlines())} lines; m-dashes {t.count(chr(8212))}, arrows {t.count(chr(8594))}, unfilled {t.count('@@')}")
    for ln in t.splitlines():
        if ln.startswith("## ") or ln.startswith("# "):
            print("  ", ln)


if __name__ == "__main__":
    main()
