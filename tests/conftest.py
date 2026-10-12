"""pytest setup for the reference checks: the repository root on sys.path, matplotlib without a display, and two helpers.

Run from the repository root:  python -m pytest tests
"""
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE = os.path.join(ROOT, "tests", "reference")
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.environ.setdefault("MPLBACKEND", "Agg")


def run_script(args, env=None, timeout=600):
    """Run one script from the repository root, as run_all.py does, with its files sent to a scratch folder.
    Returns the combined stdout and stderr as text."""
    scratch = os.path.join(ROOT, "runs", "_pytest_scratch")
    os.makedirs(scratch, exist_ok=True)
    full_env = dict(os.environ, MPLBACKEND="Agg", ORBITAL_RING_OUT=scratch)
    full_env.update(env or {})
    r = subprocess.run([sys.executable] + list(args), cwd=ROOT, env=full_env, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", timeout=timeout)
    return r.stdout + r.stderr


def reference_lines(name):
    """The headline lines of the book's log <name>.txt, copied verbatim from the reference run."""
    with open(os.path.join(REFERENCE, name + ".txt"), encoding="utf-8") as f:
        return [l.rstrip("\n") for l in f if l.strip()]


def saved_log(name):
    """The log <name>.txt of the reader's own full run (python run_all.py), from ORBITAL_RING_OUT (default runs/)
    under the repository root; the test is skipped when it is not there."""
    out = os.environ.get("ORBITAL_RING_OUT", "runs")
    path = out if os.path.isabs(out) else os.path.join(ROOT, out)
    path = os.path.join(path, name + ".txt")
    if not os.path.isfile(path):
        pytest.skip(f"{name}.txt is not in the output folder; run  python run_all.py --only={name}  first "
                    "(this one is too slow to run inside the test)")
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read()


def assert_lines_present(text, name):
    """Every reference line of <name> must appear, character for character, in the text."""
    have = {l.rstrip() for l in text.split("\n")}
    missing = [l for l in reference_lines(name) if l not in have]
    assert not missing, f"{name}: {len(missing)} reference line(s) not found in the output:\n  " + "\n  ".join(missing)
