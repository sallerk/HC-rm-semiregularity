Supporting code and data for the note
"The semiregularity step in Markman's approach to Weil classes for quartic CM fields"
================================================================================================================

Status: accompanies the note (K. Saller, 2026), DOI 10.5281/zenodo.23682740; not peer reviewed.

Setting. X is a complex abelian variety of dimension n = 2m with real multiplication by a real quadratic field
F = Q(sqrt t) of discriminant d_F, Theta = Theta_1 + Theta_2 a compatible polarization, split along the two real
places, q > 0 rational or, where the note says so, q in F totally positive, and K = F(sqrt -q). Markman's secant
space B in H^ev(X, Q) has dimension 4 and, for rational q, the basis alpha, beta, alpha~, beta~; P = <alpha, beta>
and P~ = <alpha~, beta~> are the secant planes of the two imaginary quadratic subfields of K. With
c = x alpha + y beta + z alpha~ + w sqrt t beta~, r = x + z, v = z - x and g = y + w sqrt t, the rank-one cone of B
is q r v = N(g). HT^k(X) is the space of polyvector fields of degree k, K(c) in HT^2(X) the annihilator of a class
c under the Clifford action, r_k(c) the rank of the contraction with c on HT^k, and chi the Euler form (Mukai
pairing). For a lattice Lambda = H_1(X, Z) the integral classes are those of H^*(X, Z) = wedge^* Hom(Lambda, Z). The
twisted secant spaces are the e^{-L_h} B, h in F, and Pi is the quartic invariant of Section 6 of the note.

Role of the computations. The note's Appendix B ("Computations") lists six items with lettered sub-checks; Section 3
below says which script makes each check. In the note's description, where the text cites the appendix, the check
confirms a computation given there or carries out one that the text describes; the other checks test proved
statements on examples, and the searches and enumerations are complete only within the stated bounds. The
statements that the note marks as computed, or whose support it calls computational, are:
  - the sharpness of Theorem 5.6 on the 12 free lattices (Introduction, and after Theorem 5.6): check (3e);
  - the remark after Theorem 5.6 that on random lattices Lemma 5.5 is often an equality with nonzero glue: (3c);
  - the remark after Corollary 5.8: for t = 5 and q = 1, on the principal lattice, the integral classes of B with
    chi <= 14 are +-alpha, +-beta and the four classes +-(beta +- sqrt5 beta~)/2: (3d);
  - Example 6.8, the four classes at the bound of Theorem 6.7, with the factors lambda = 1/3, 1/2, 1, 1: (4d);
  - the evidence for Conjecture 6.10, the content of Pi on 112 lattices: (4e);
  - after Theorem 6.13, the classes with chi = 10 and Pi = 20 for d_F = 5, and the finite search that found no
    class at the bound for d_F = 8, 12, 13: (4f);
  - Example 7.10, and the integral points of the cone with c_1 not in Q Theta: (5c);
  - after Theorem 8.4, the stabilizer dimensions 70 at five rational points of B off the cone and 143 at alpha,
    beta, alpha~, beta~: (6a).
The other checks confirm computations that the note carries out by hand (for instance the identities in Section 6
that the note says are "also checked by computer", the identities of Remark 6.6, and Proposition 6.15 for
n = 4, 6, 8), or test proved statements on examples.

Not checked by any script: the proofs (the scripts work in explicit models, at finitely many points, lattices and
parameters); the cited results (among them those of Markman, Perry, Orlov, Sato-Kimura and Bass, and the
statements quoted from [BMB], [MS26] and [OAWP]); anything about Ext groups, obstruction maps or the
semiregularity of actual sheaves, since the scripts compute only classes, annihilators, ranks of linear maps on
cohomology, Euler forms and integrality; and the arguments of Section 3 (genus four) beyond the computations of
Example 3.3, of Section 4 (the count) beyond its numerical inputs, and of Appendix A.

Labels. The headers and output lines of the scripts in rm-semiregularity use the numbering of the
working document NOTE.txt (not part of this repository), for example "NOTE Theorem 9.8" or "Proposition 9.11"; the
scripts in rm-twisted use the numbering of a second working document, NOTE.md (not part of this repository), or
earlier numberings of the note itself. Section 4 below translates these labels into the numbering of the note.
The letters of the checks, such as (G) or (S6), are those of the scripts; check (kx) of the note is sub-check (x)
of item (k) of its Appendix B.

Licence: the code (*.py, *.sh) is under the MIT licence; the outputs and this text are under CC BY 4.0 (see the
LICENSE files of the repository). No third-party files are included.


1. Requirements
---------------
Python 3 with python-flint (for shape_search.py, cone6.py, split.py and, through twist_lib.py or twist_scan.py,
which load shape_search.py, for every script in rm-twisted except bmb_m.py and bmb_m_all.py), SymPy (semireg.py,
cone6.py, split.py, cm_ranks.py, bmb_m.py, bmb_m_all.py) and NumPy (imported by semireg.py, which cone6.py,
split.py and cm_ranks.py load). compare_outputs.py uses the standard library only. Tested versions: Python 3.12.9,
python-flint 0.9.0, NumPy 2.1.3, SymPy 1.13.1, on Windows 10; these are the versions that Appendix B of the note
names. The arithmetic is exact (integers, fractions, real quadratic fields), except that floating point bounds the
ranges of the enumerations and selects the test values of the finite searches (every class found is checked
exactly), and that cm_ranks.py computes the coordinates c_ef of its classes in floating point, so its count of
nonzero c_ef is a sanity check, as its docstring says.


2. Layout and how to run
------------------------
rm-semiregularity/code/      semireg.py, shape_search.py, cone6.py, split.py
rm-semiregularity/outputs/   their reference outputs: NAME.json (the results) and NAME.out (the screen output)
rm-twisted/code/      twist_lib.py, twist_scan.py (modules, see below), check_algebra.py, check_sym2.py,
                      check_planes.py, check_product.py, examples.py, ex68_lambda.py, nonproduct.py, multiplier.py,
                      principal_chi14.py, cm_fields.py, cm_ranks.py, bmb_m.py, bmb_m_all.py, p88_shape.py,
                      p88_check.py
rm-twisted/outputs/   their reference outputs: NAME.log (the screen output), and p88_shape.txt, p88_check.txt (the
                      output files of these two scripts, which print the same lines to the screen)
run_all.sh, compare_outputs.py

twist_lib.py and cm_ranks.py load shape_search.py and semireg.py by the relative path ../../rm-semiregularity/code,
so the folder layout must be kept. Imports: cone6.py imports semireg.py and shape_search.py; split.py imports
semireg.py, shape_search.py and cone6.py; twist_lib.py loads shape_search.py; check_sym2.py, check_algebra.py and
principal_chi14.py load twist_lib.py; check_planes.py, check_product.py, examples.py, nonproduct.py, multiplier.py,
cm_fields.py and p88_shape.py load check_sym2.py; cm_fields.py also loads twist_scan.py (for its lattice
enumeration routine), p88_shape.py loads nonproduct.py, p88_check.py loads p88_shape.py, cm_ranks.py loads
semireg.py and twist_lib.py, and ex68_lambda.py loads examples.py from the folder given as its argument.
twist_scan.py is included as a module of cm_fields.py only; its own script mode, a heuristic search that the note
does not use, is not run here.

From the top folder of this repository:
    sh run_all.sh rerun
    python3 -I compare_outputs.py rm-semiregularity/outputs rerun/rm-semiregularity/outputs
    python3 -I compare_outputs.py rm-twisted/outputs rerun/rm-twisted/outputs
run_all.sh runs the scripts one after the other, each as python3 -I -B SCRIPT [ARGUMENTS], and writes the results
into rerun/rm-semiregularity/outputs and rerun/rm-twisted/outputs. The arguments are: the output file
NAME.json for the four scripts of rm-semiregularity (their screen output goes to NAME.out); none for
check_algebra.py, check_sym2.py, check_planes.py, check_product.py, examples.py, nonproduct.py, bmb_m.py,
bmb_m_all.py, principal_chi14.py and cm_ranks.py; 3 for multiplier.py (three random lattices for each field and
each N, the default); 4 6 16 for cm_fields.py (MMAX = 4, QD = 6, EMAX = 16, the search of check (4f)); the folder
rm-twisted/code for ex68_lambda.py; and the output file NAME.txt for p88_shape.py and p88_check.py. The shell
script calls python3; set the environment variable PYTHON to use another interpreter name (for example
PYTHON=python). Every script except nonproduct.py and multiplier.py, which print data only, ends with the line
"ALL CHECKS PASS: True" (p88_shape.py: "ALL CONTROLS PASS: True") when all its checks pass; run_all.sh lists these
lines at the end. compare_outputs.py compares each .json, .out, .log and .txt file with the reference; it sets
aside the running times, which the scripts of rm-semiregularity write into their .json files (entry
"seconds") and into some screen lines (for example "[1.2s]"). The random tests use fixed seeds. Its expected last
line is "ALL EQUAL"; it also ignores the line ends, which are LF in rm-semiregularity/outputs and
CRLF (Windows) in rm-twisted/outputs. Running times in the rerun of Section 5: about 6 minutes in all without
cm_fields.py (shape_search.py about 4 minutes, every other script at most a minute); cm_fields.py with the
arguments 4 6 16 about 28 minutes.

Reference outputs. The files in rm-semiregularity/outputs and six files in rm-twisted/outputs (cm_ranks.log,
bmb_m.log, bmb_m_all.log, principal_chi14.log, p88_shape.txt, p88_check.txt) were made in the working tree between
20 and 26 October 2026; the platform and the versions of those runs are not recorded in them. cm_fields.log is the
output of the rerun of Section 5; the output made in the working tree on 24 October 2026 differs from it only in
the printing of two values of q (Section 5). For check_algebra.py, check_sym2.py, check_planes.py,
check_product.py, examples.py, nonproduct.py, multiplier.py and ex68_lambda.py no output had been kept; their
reference outputs (NAME.log) were made from the files of this repository when it was assembled, on 27 October 2026
(Section 5), so for these eight scripts the comparison of Section 5 shows reproducibility only.

Changes against the working tree: in the working tree the folder rm-semiregularity lies one level deeper, and the
relative paths to it in twist_lib.py (line 25), cm_ranks.py (line 46) and twist_scan.py (lines 3 and 23) are
changed to the layout of this repository. In principal_chi14.py (line 17) and twist_scan.py (line 23) an absolute
path of the author's machine, to twist_lib.py and to shape_search.py, is replaced by the same file named relative to
the script (TLP = __file__.replace("principal_chi14.py", "twist_lib.py") and
SRC = __file__.replace("twist_scan.py", "../../rm-semiregularity/code/shape_search.py")).
ex68_lambda.py was written in a scratch folder during a revision of the note, not in the code folders of the
working tree; it is copied unchanged. run_all.sh and compare_outputs.py are new. All other files are unchanged.


3. The scripts and the items of the note's Appendix B
-----------------------------------------------------
Statement numbers are those of the note. "S/" stands for rm-semiregularity/code/, "T/" for
rm-twisted/code/; the output lines are those of the reference output of the script named, which begin with the
letter of the check.

Item 1 (Section 2) -- S/semireg.py, for (n, t, q) = (4, 3, 1), (4, 5, 1), (4, 5, 3), (6, 57, 3), (8, 11, 1):
  K(B) = R + G of dimension n^2/2 and dim K(P) = n^2, by an explicit basis and ranks modulo a prime: check (A).
  dim K(c) = n^2/2 at two random rational points of B in each case and at Markman's beta' (n = 4), with the exact
  identity that puts beta' in B: (B). The ranks r_k(c) on every HT^k: on P and at a random point of B for
  n = 4, 6, 8, by the block formula, checked against direct ranks modulo a prime for n = 4, 6 in every degree and
  for n = 8 in degree 2, and on P for n = 3 by direct ranks: (G); on the rank-one cone for n = 4, 6, 8, at four
  pure tensors each, by the block formula against the closed formula of Theorem 2.3(iv): (T), lines "rank profile
  at four pure tensors"; direct ranks modulo a prime at cone points in every degree for n = 6 are in S/cone6.py
  (U), and for n = 4 only dim K(c) = 10 at two cone points, in (T). For n = 4 with c symbolic, the determinants of
  the mixed and same-block matrices, which give dim K(c) = 16, 10, 8: (T), first line. The injectivity on HT^1
  (Lemma 4.4): (E), for t = 3, 5 and q = 1, 2, 3. The Gram matrices of Proposition 2.1 for n = 2, 4, 6, 8, q
  symbolic: (G), line "Mukai Gram matrices on B", and (D) for n = 4; chi(c) = 32 q^2 (q r v - N(g)) for n = 6:
  (T), line "n=6: Mukai form on B".
Item 2 (Example 3.3) -- S/semireg.py (C), first line: genus 4, the differentials w_1, ..., w_4, the single
  relation Q = -(w_1 w_4 - w_2 w_3) among the ten products, of rank 4, and its components in Sym^2 V_1, V_1 V_2
  (none) and Sym^2 V_2. The other lines of (C) (a control curve and the other exponent patterns) are not used by
  the note.
Item 3 (Section 5) -- S/shape_search.py, except (3d) in part and (3e):
  (a) (H): free O_F-lattices, 25 polarization types and 552 values of q for each of t = 5, 2, 3, 13, least chi of
      a class off P u P~ with chi <= 14 (12, for t = 5); (I): cyclic quotients, (J): quotients by elementary
      abelian subgroups of X[p], p = 2, 3, 5 (Smith normal forms).
  (b) (L): the criterion for odd classes and chi = chi(c_ev) + 4 Tr(a), 2200 tests on 88 lattices; (M): the
      inequality that excludes +-E_1 (Proposition 5.3), 105 lattices; (N): 4 d_F | chi(E_2), and no nonzero
      integral class of B with chi <= 6 for five values of q, 430 lattices.
  (c) (O): the steps of the proof of Theorem 5.4 for 103 pairs (prime, lattice); the counts "Lemma 9.9 tight with
      nonzero glue" add up to 35.
  (d) (R): the classes with c_1 in Q Theta and chi <= 14 on principal free lattices (t = 5, 2, 3, 13, q = 1);
      T/principal_chi14.py: all classes with chi <= 14 on the same lattices, with the four further classes for
      t = 5 and their r = 0 and indefinite g.
  (e) T/p88_check.py, check (A) (output p88_check.txt, one line per lattice, 171 lattices): the content of
      varpi(4N(m) - rho^2) on the lattice of integral classes Q_m + rho Theta_1 Theta_2, divided by 4 d_F; on the
      12 lines "free(...)" (two polarization types for each t) it is 1.
Item 4 (Section 6) -- rm-twisted:
  (a) T/check_algebra.py (A1)-(A4): the matrices of alpha, beta, alpha~, sqrt t beta~ (against shape_search.py),
      the twist X -> E_h X sigma(E_h)^T, the formula for chi, and degreewise integrality against integrality of
      the form. (A5), (A6) check a normal form and the cone of the working document and are not cited.
  (b) T/check_sym2.py (S1)-(S3): tr eps = 0, (varpi/4) tr eps^2 = chi, det, and the characteristic polynomial of
      eps^2 on random twisted secant classes; T/check_planes.py (K1)-(K3): the planes H_+, H_- with eps^2 = m_+-
      there, and eps = 0 on H_+ for c in P and on H_- for c in P~.
  (c) T/check_product.py (P1)-(P3): varpi = 4 W_A W_B, the integrality of S_AB and the formula for Pi, on 9 free
      polarization types for each of t = 5, 2, 3, 13 with 40 random classes each.
  (d) T/examples.py: the four classes of Example 6.8 (integral, off P u P~ and the cone, chi(c_P), chi(c_P~),
      Pi = 4 d_F, chi = s(d_F)); T/ex68_lambda.py prints the factors lambda = 1/3, 1/2, 1, 1 of Example 6.8 (the
      inverse of the content of Tr E) and runs examples.py again.
  (e) The content of Pi on L_R(Lambda): T/check_sym2.py (S6): 36 free lattices (content/(4 d_F) = 1) and 16
      random lattices; T/multiplier.py: 36 random lattices O_F^4 + sum Z gamma_i/N, N = 2, 3, 4; T/nonproduct.py:
      24 lattices, namely the 8 lattices "O^4, f=2" and "O^4, f=3" of Lemma 6.9 (content/(4 d_F) = f^2) and 16
      variants ("+e24") and superlattices ("+ Z gamma/f"). Together 112 lattices; on the 52 random ones the values
      of content/(4 d_F) are 1, 2, 3, 4, 8, 9, 16, 27, 32, 64.
  (f) T/cm_fields.py 4 6 16: (M1) the matrices of B(Theta, q); (M2) the identities of Proposition 6.12 on random
      twisted classes; (M3) the finite search over q and h on two free lattices for each of t = 5, 2, 3, 13, with
      the search bound s_cm(d_F) + 8 and the least chi found (10, 14, 22, 18 for d_F = 5, 8, 12, 13).
  (g) T/cm_ranks.py (i)-(iii): Proposition 6.15 for the five values of (t, q) of the note, and the biquadratic
      control q = 3 + sqrt 5.
  (h) T/bmb_m.py (on twisted secant spaces) and T/bmb_m_all.py (with all nine coordinates of R symbolic): the
      identities of Remark 6.6.
  (i) T/p88_check.py, check (B): for 11 totally positive q in F on each of the 171 lattices of (3e), the odd
      integral classes (a in O(Lambda), chi = 4 Tr(a)) and the least chi <= 14 of a nonzero integral class, by t
      and kind of K (last lines of p88_check.txt: 8 for biquadratic K, no entry for the other fields);
      T/p88_shape.py: the search for the shape +-(E_1^q + rho Theta_1 Theta_2) on 430 lattices with
      t = 2, 7, 14, 17 (line TOTAL: 0 solutions), with control 1 (the degree-4 part of alpha for q = 1), control 2
      (chi = 6 and Pi = 1 for the shape) and control 3 (a direct search on six lattices).
Item 5 (Section 7) -- S/cone6.py, except (5d):
  (a) (U): ranks (1, 12, 48, 74, 48, 12, 1) at cone points (block formula and direct ranks modulo a prime), the
      count forcing e_1 = 12, e_3 = 74, the ranks (1, 6, 6, 1) on one block, the inclusion of the annihilators on
      3584 products, and dim K(c) = 18 at rational cone points for t = 57, q = 3.
  (b) (V): e^Theta f(Theta) = alpha_0 + beta_0, the ranks and Mukai square of ch(E'), and the operator A.
  (c) (W): the identities of Proposition 7.9(a)-(c) (first line), and the enumeration of integral cone points
      below the height bound H on principal lattices and on two quotient lattices for (t, q) = (57, 3); among the
      classes listed are those of Example 7.10 (t = 5, q = 1: r = 1, g = 2, chi(F) = 4; t = 57, q = 3: r = 1,
      g = 15, chi(F) = 2025).
  (d) S/shape_search.py (P), the four lines "n=6": the identities of Lemma 7.2 and 32 | chi for t = 5, 57, 13 and
      q = 1, 1/3, 3/2. (The lines "n=8" of (P) check a statement of the working document that the note does not
      contain.)
Item 6 (Section 8) -- S/split.py: (a) (Y); (b) (Z), six configurations; (c) (AA), with eleven Z-bases of orders
  and N = 1, ..., 9, that is d = N - 1 = 0, ..., 8.

Checks of the included scripts that Appendix B does not use: semireg.py (F), (K), (Q), (S) and the last line of
(T) (classes built from line bundles); the lines "n=8" of shape_search.py (P); cone6.py (X); check_algebra.py (A5),
(A6); check_sym2.py (S4), (S5). They are rerun with the rest.


4. Labels of the working documents and numbers of the note
-----------------------------------------------------------
rm-semiregularity (NOTE.txt; "NOTE" in the scripts):
  Theorem 3.1(i)-(iii)                          Theorem 2.3(i)-(iii)
  Theorem 3.1(iv) (the rank-one cone)           Theorem 2.3(ii) (the cone q r v = N(g), Section 2)
  Proposition 3.4 (weight blocks)               the proof of Theorem 2.3, and Theorem 2.3(iv)
  Section 5, Example 5.3                        Section 3, Example 3.3
  Propositions 5.6, 5.7, Example 5.8,           not in the note (mixed Jacobians: Remark 3.6)
    Remark 5.9 (symmetric points)
  Section 7                                     Section 4
  Proposition 7.1                               Proposition 2.1
  Lemma 7.2, Proposition 7.3                    Lemma 4.4, Lemma 4.3
  Theorem 7.4, 7.4(b)                           Theorem 4.5, 4.5(ii)
  Proposition 8.1(e)                            Corollary 5.8
  Section 9                                     Sections 5, 7 and 8
  Theorem 9.2 (free lattices)                   the free case of Corollary 5.7
  Theorem 9.4                                   Theorem 5.2
  Proposition 9.5                               Proposition 5.3
  Proposition 9.6 (the line Q Theta_1 Theta_2)  Theorem 5.4; for dim K(Theta_1 Theta_2) = 10, Theorem 2.3(ii)
  Theorem 9.8                                   Theorem 5.4
  Lemma 9.9                                     Lemma 5.5
  Corollary 9.10                                Corollary 5.7
  Proposition 9.11                              Lemma 7.2 (n = 6; the case n = 8 is not in the note)
  Corollary 9.12                                Remark 4.9 and Proposition 7.3 (the all-degree form)
  Theorem 9.14                                  not in the note
  Theorem 9.15                                  Theorem 7.4
  Theorem 9.16                                  Theorem 7.5
  Proposition 9.17                              Proposition 7.9
  Example 9.18                                  Example 7.10
  Proposition 9.19                              the members X_0 (x) O of Remarks 8.6 and 8.7
  Lemma 9.20                                    Lemma 8.1
  Theorem 9.21                                  Theorem 8.4
  Corollary 9.22                                Remark 8.6
  Proposition 9.23                              Remark 8.7
  Propositions 10.4, 10.5, 10.7                 not in the note
  Section 11 (OPEN), items 11(f), 11(l)         Section 9 (open cases)
rm-twisted (NOTE.md):
  Section 2 (the Hermitian model)               Section 6, Lemma 6.2
  Section 3 (the RM-surface model)              Section 6, Lemma 6.3 and Proposition 6.4
  Proposition 3.1                               Proposition 6.4 (the planes H_+, H_-)
  Theorem 3.4                                   Theorem 6.5
  Section 3.3; the remark on [BMB]              Theorem 6.7; Remark 6.6
  Theorem 3.6                                   Theorem 6.7
  Section 3.5                                   Lemma 6.9, Conjecture 6.10
  Section 3.6, Theorem 3.9 (the bound s_cm)     Propositions 6.11-6.12, Theorem 6.13
rm-twisted, earlier numberings of the note:
  Corollary 5.7 (principal_chi14.py)            Corollary 5.8
  "Cor 5.6 for every quartic CM field"          Proposition 6.11
    (p88_shape.py)
  Lemma 6.7 (p88_shape.py, p88_check.py)        Lemma 6.9
  Theorem 2.3 (cm_ranks.py), Lemma 6.2(a), Lemma 6.3, Proposition 6.4 (bmb_m_all.py) and Example 6.8
  (ex68_lambda.py) are already numbered as in the note.
In bmb_m.py, "[BMB] Proposition 19.36" is the numbering of [BMB]. In semireg.py, split.py and twist_scan.py,
"[M2]" is Markman's paper [MarW] of the note (arXiv:2502.03415), and "arXiv:2509.23079" in semireg.py is [MarRM].


5. Rerun
--------
When this repository was assembled, on 27 October 2026, it was run on the author's machine (Windows 10, with the
versions listed in Section 1). First the eight scripts without kept outputs (Section 2) were run from the
repository, and their screen outputs became their reference outputs. Then run_all.sh was run into a folder outside
the repository (about 33 minutes, of which about 28 for cm_fields.py), and compare_outputs.py compared the two
output folders with the reference outputs:
  - rm-semiregularity/outputs: all eight files equal, "ALL EQUAL". The files differ from the reference outputs
    only in the running times and in the line ends, which compare_outputs.py sets aside.
  - rm-twisted/outputs: 14 of the 15 files equal, all 14 byte for byte. cm_fields.log differed from the output
    made in the working tree on 24 October 2026 in two lines, and there only in the sign "+": that output predates
    a fix of the printing of q (Q2.__repr__ in twist_lib.py), which dropped the "+", so that it printed "11/41*r3"
    for 11/4+1*r3 and "5/61/6*r13" for 5/6+1/6*r13. All other characters agree, the computed values among them.
    The reference cm_fields.log is now the output of this rerun.
Every script ran to the end and printed its pass line ("ALL CHECKS PASS: True"; p88_shape.py: "ALL CONTROLS PASS:
True"), except nonproduct.py and multiplier.py, which print data only. The outputs of the working tree in
rm-twisted/outputs other than cm_fields.log (cm_ranks.log, bmb_m.log, bmb_m_all.log, principal_chi14.log,
p88_shape.txt, p88_check.txt) were reproduced byte for byte. For the eight scripts whose reference outputs were
made at assembly, the comparison shows only that two runs on the same machine agree.


6. What the computations do not show
------------------------------------
- The scripts work in explicit models: in dimensions n = 3, 4, 6, 8 for the annihilators and ranks, on finitely
  many lattices, polarization types, values of q and h, and random classes for the integrality, and below stated
  bounds for the enumerations and searches. The statements in general are proved in the note.
- They compute classes in cohomology and linear algebra on it. That a sheaf with a given Chern character exists,
  what its Ext groups and obstruction map are, and whether it is semiregular, is not computed; the count of
  Section 4 and the arguments of Section 3 and Appendix A are not checked.
- The cited theorems are not checked.
- Conjecture 6.10 is tested on 112 lattices; the computations do not prove it.
- All scripts were written with the assistance of an AI system (Claude, by Anthropic).
