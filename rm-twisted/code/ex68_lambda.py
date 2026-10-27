"""Print the normalizing factor lambda of Example 6.8 (rm-semiregularity note) for the four rows.

The rows are built in rm-twisted/code/examples.py: E = delta(g1 e12 + g1 g2 e34), Theta = Tr E divided by its
content cth, so that Theta is primitive. Hence lambda = 1/cth. Also re-runs the checks of examples.py.
usage: python -I -B ex68_lambda.py RM_TWISTED_CODE_DIR"""
import importlib.util
import os
import sys
from fractions import Fraction as Fr

code_dir = sys.argv[1]
spec = importlib.util.spec_from_file_location("examples", os.path.join(code_dir, "examples.py"))
EX = importlib.util.module_from_spec(spec)
spec.loader.exec_module(EX)
SS = EX.SS
for t, ((a1, b1), (a2, b2)), q, h, u, chi in EX.EXAMPLES:
    g1, g2 = SS.Fe(t, a1, b1), SS.Fe(t, a2, b2)
    dl = SS.delta(t)
    P, s0 = SS.setup(t, dl * g1, dl * g1 * g2)
    cth = SS.content(P[(1, 0)])
    print("t=%d: content of Tr E = %s, lambda = %s" % (t, cth, Fr(1) / Fr(cth)))
EX.main()
