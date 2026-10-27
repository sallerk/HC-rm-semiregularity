"""
compare_outputs.py -- compare the result files of a rerun with the reference outputs.

Usage:  python -I compare_outputs.py REFERENCE_DIR RERUN_DIR
Every .json, .out, .log and .txt file in REFERENCE_DIR is compared with the file of the same name in RERUN_DIR:
.json files as parsed JSON values, without their top-level entry "seconds" (the running times); the other files
line by line, with trailing white space and the line-end convention ignored, and with running times printed in
the form "[1.2s]" or "(0.1s)" (a number followed by s, in brackets or parentheses) replaced by a placeholder.
The random tests of the scripts use fixed seeds, so nothing else is set aside.  Each file is reported as "equal"
or "DIFFERENT" (with the first line that differs).  Expected last line: "ALL EQUAL".
"""

import json
import os
import re
import sys

TIME = re.compile(r"([\[(])\d+(?:\.\d+)?s([\])])")
KINDS = (".json", ".out", ".log", ".txt")


def text_lines(path):
    with open(path, encoding="utf-8") as f:
        return [TIME.sub(r"\1<time>\2", line.rstrip()) for line in f]


def json_value(path):
    with open(path, encoding="utf-8") as f:
        value = json.load(f)
    if isinstance(value, dict):
        value.pop("seconds", None)
    return value


def first_difference(A, B):
    for i, (a, b) in enumerate(zip(A, B)):
        if a != b:
            return "line %d differs" % (i + 1)
    return "lengths %d and %d" % (len(A), len(B))


def main():
    ref, new = sys.argv[1], sys.argv[2]
    allequal = True
    names = sorted(n for n in os.listdir(ref) if n.endswith(KINDS))
    for name in names:
        a = os.path.join(ref, name)
        b = os.path.join(new, name)
        if not os.path.exists(b):
            print("%-24s missing in rerun" % name)
            allequal = False
            continue
        if name.endswith(".json"):
            A, B = json_value(a), json_value(b)
            note = "" if A == B else "  (JSON values differ)"
        else:
            A, B = text_lines(a), text_lines(b)
            note = "" if A == B else "  (%s)" % first_difference(A, B)
        eq = A == B
        print("%-24s %s%s" % (name, "equal" if eq else "DIFFERENT", note))
        allequal &= eq
    extra = sorted(n for n in os.listdir(new) if n.endswith(KINDS) and n not in names)
    for name in extra:
        print("%-24s only in rerun" % name)
    print("ALL EQUAL" if allequal and not extra else "DIFFERENCES FOUND")


if __name__ == "__main__":
    main()
