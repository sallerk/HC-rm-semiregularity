#!/bin/sh
# run_all.sh -- rerun the checks of this repository; the results go into the folder given as argument
# (default: rerun).
# Usage (from the top folder of the repository):   sh run_all.sh rerun
# Then:   python3 -I compare_outputs.py rm-semiregularity/outputs rerun/rm-semiregularity/outputs
#         python3 -I compare_outputs.py rm-twisted/outputs rerun/rm-twisted/outputs
# The interpreter is python3, or the value of the environment variable PYTHON (for example PYTHON=python).
# The scripts run one after the other: about 6 minutes in all, plus about 28 minutes for cm_fields.py, whose
# search (arguments 4 6 16) is the longest (see README.txt, Sections 2 and 5). A script that fails does not stop
# the others, and the last lines list the scripts that report that all their checks pass.
OUT=${1:-rerun}
PY=${PYTHON:-python3}
T=$(dirname "$0")
S="$T/rm-semiregularity/code"
R="$T/rm-twisted/code"
OS="$OUT/rm-semiregularity/outputs"
OR="$OUT/rm-twisted/outputs"
mkdir -p "$OS" "$OR" || exit 1
for s in semireg shape_search cone6 split; do
  "$PY" -I -B "$S/$s.py" "$OS/$s.json" > "$OS/$s.out"
done
for s in check_algebra check_sym2 check_planes check_product examples nonproduct bmb_m bmb_m_all \
         principal_chi14 cm_ranks; do
  "$PY" -I -B "$R/$s.py" > "$OR/$s.log"
done
"$PY" -I -B "$R/multiplier.py" 3 > "$OR/multiplier.log"
"$PY" -I -B "$R/cm_fields.py" 4 6 16 > "$OR/cm_fields.log"
"$PY" -I -B "$R/ex68_lambda.py" "$R" > "$OR/ex68_lambda.log"
"$PY" -I -B "$R/p88_shape.py" "$OR/p88_shape.txt" > /dev/null
"$PY" -I -B "$R/p88_check.py" "$OR/p88_check.txt" > /dev/null
grep -H "ALL CHECKS PASS\|ALL CONTROLS PASS" "$OS"/*.out "$OR"/*.log "$OR"/*.txt
