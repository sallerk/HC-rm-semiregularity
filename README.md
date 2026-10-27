# HC-rm-semiregularity

Code and data supporting the note *The semiregularity step in Markman's approach to Weil classes for quartic CM
fields* (K. Saller, 2026), DOI: [10.5281/zenodo.23682740](https://doi.org/10.5281/zenodo.23682740).

**Subject.** Markman reduces the algebraicity of the Weil classes on abelian varieties of Weil type for a CM field
K = F(sqrt(-q)) to the semiregularity of an object built from two sheaves whose Chern characters lie in a secant
space B, on an abelian variety X with real multiplication by F. For imaginary quadratic K he proves the needed
semiregularity on Jacobians of genus three, through two properties of a sheaf that is equivariant for a finite
group of translations. The note studies this step when F is real quadratic, so that K is a quartic CM field and B
has dimension four. It shows that the mechanism does not transfer in two cases: on abelian fourfolds with real
multiplication, for finite groups of translations and any lattice, and for the actions conjugated by a line-bundle
twist (graph-type actions) when the lattice of the quotient is locally split. The tools are a count of invariant
Ext groups and the integrality of the classes of B, in particular a divisibility of the integral classes of degree
four on every lattice and a quartic invariant Pi on locally split lattices. The note also treats abelian sixfolds
(a rank-one cone in B and a Galois obstruction) and Markman's sheaves on Jacobians of genus four. The semiregularity
question itself stays open, and none of the results concerns the Hodge conjecture.

**Role of the computations.** The note's Appendix B ("Computations") lists six items. Most of them confirm
computations that the note carries out by hand, or test proved statements on examples; the searches and
enumerations are complete only within the bounds stated there. Some statements rest on a computation, and the note
says so where they are stated: the sharpness of Theorem 5.6 on twelve free lattices, the list of small integral
classes in the remark after Corollary 5.8, Example 6.8, the evidence for Conjecture 6.10 (112 lattices), the class
attaining the bound of Theorem 6.13 for d_F = 5 and the search that found none at the bound for d_F = 8, 12, 13,
Example 7.10, and the dimensions of stabilizers stated after Theorem 8.4.
[`README.txt`](README.txt) gives the requirements, the commands, what each script checks, how the scripts' labels
translate into the numbering of the note, and the result of the rerun made when this repository was assembled.

**Status:** not peer reviewed.

**Disclosure:** the arguments, the computations and the code were produced with the assistance of an AI system
(Claude, by Anthropic).

## Licences

- **Code** (`*.py`, `*.sh`): MIT License, see [`LICENSE`](LICENSE).
- **Everything else** (README files and computed outputs): Creative Commons Attribution 4.0 International
  (CC BY 4.0), see [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
- No third-party files are included.
