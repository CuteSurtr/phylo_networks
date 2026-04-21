# phylomath

A from scratch library for the mathematics of phylogenetic trees and networks. The library implements the combinatorial, geometric, and algebraic structures that underpin modern phylogenetics: split systems and the Buneman compatibility theorem, enumeration of rooted binary trees, Robinson Foulds and quartet distances, the BHV tree space of Billera, Holmes, and Vogtmann with its CAT(0) geometry, the exact Owen Provan polynomial time geodesic algorithm, Frechet means via Bacak proximal iteration, tanglegram crossing minimization, the Neighbor Net circular split algorithm, rearrangement graphs on NNI moves, and the tropical metric on tree space with Fermat Weber points. All core algorithms are verified on exact and sampled data and include golden tests pinned to published values.

## Foundations

Phylogenetic trees on a set of taxa form a space whose mathematical structure is unusually rich. Billera, Holmes, and Vogtmann showed in 2001 that the set of rooted binary phylogenetic trees with positive edge lengths on n leaves forms a metric space made of (2n minus 3) double factorial orthants of dimension n minus 2, glued along faces corresponding to shared splits. This space is CAT(0), meaning it has unique geodesics and well defined means. Owen and Provan gave the first polynomial time algorithm for computing BHV geodesics in 2011. Tropical geometry provides an alternative metric on the same space via the Hilbert projective distance on dissimilarity vectors. Rearrangement graphs, tanglegrams, and Neighbor Net together form the combinatorial and graph theoretic half of the subject, while phylogenetic networks extend the tree model to directed acyclic graphs that accommodate hybridization and horizontal gene transfer.

## Core algorithms

Trees and splits are implemented in `tree.py` and `split.py` with a minimal Newick parser, postorder traversal, the Buneman compatibility predicate, and the canonical split set of a tree. Robinson Foulds, weighted Robinson Foulds, and the quartet distance are in `rf.py`. Enumeration of all rooted binary semi labeled trees on n leaves and the NNI neighborhood of a tree are in `enumerate.py`. The BHV tree space is built in `bhv.py` with the cone path distance and the Vogtmann common split decomposition. The full Owen Provan GTP algorithm is in `gtp.py`, including the closure lattice construction, enumeration of maximal chains, the path space length formula, and a geodesic point walker that evaluates the unique geodesic between two trees at any fraction t in zero to one, used by `frechet.py` to compute the Frechet mean via Bacak proximal iteration in a CAT(0) space.

## Networks, tanglegrams, and rearrangement graphs

Phylogenetic networks are modeled in `network.py` as rooted directed acyclic graphs with explicit reticulation nodes. The module provides predicates for tree based networks, counts of reticulations, and the Bryant and Moulton Neighbor Net agglomerative algorithm that builds a circular split system from a distance matrix. Tanglegrams between two leaf labeled binary trees are handled in `tanglegram.py` with a merge sort based inversion count for crossings, an O(n log n) dynamic program for one tree crossing minimization, and an exponential brute force for two tree crossing minimization suitable for small trees. Rearrangement graphs on binary trees are built in `rearrange.py` using NNI moves and Robinson Foulds bounded adjacency, together with shortest NNI distance via networkx.

## Tropical geometry

The tropical metric on tree space is implemented in `tropical.py`. The module computes the dissimilarity vector of a tree (pairwise path distances across all leaf pairs), the tropical d infinity distance modulo constant shift, and the tropical Fermat Weber point as a linear program via scipy linprog. The tropical line between two points is also provided as a sampling function for visualization. Together with the BHV geodesic and the Frechet mean, this provides two genuinely different CAT(0) style average tree notions for the same input data.

## Layout

```
literature/   foundational papers (Billera Holmes Vogtmann 2001, Owen Provan 2011 GTP geodesic, Brown Owen mean and variance, Miller Owen Provan polyhedral averaging, tropical geometry of tree space, Moon counting labeled trees, Matrix tree theorem and Cayley, Allman Rhodes phylogenetics lectures, tree rearrangement uniqueness intractability, connecting BHV spaces with non identical leaves)
data/         example Newick trees (primate mitochondrial Frechet sample with five bootstrap trees)
docs/         lit_review.md (literature synthesis on BHV geometry, Owen Provan geodesics, tropical tree space, networks, rearrangements) and PLAN.md (layered roadmap)
src/phylomath/ the library (splits, trees, BHV, GTP geodesic, Frechet, networks and Neighbor Net, tanglegrams, rearrangement graphs, tropical metric, Fermat Weber, visualization)
tests/        pytest suite (forty eight tests covering splits compatibility, RF and quartet distances, enumeration counts, BHV cone path vs GTP geodesic agreement on shared split cases, Owen Provan path space formula, triangle inequality, symmetry, Pythagorean identity across orthogonal orthants, Frechet mean recovery on identical and averaged samples, NNI graph structure, tanglegram crossing count, Neighbor Net on a circular matrix, tropical distance invariance under constant shift, Fermat Weber reduction on colinear points, and golden tests pinned to published Owen Provan examples)
results/      figures and benchmark output
```

## Quick start

```bash
cd src
python3 -m phylomath.demo
```

The demo builds random rooted binary trees and validates the Schroder double factorial enumeration count, computes the RF and quartet distances between two example trees, constructs the NNI rearrangement graph for n equals four and five, demonstrates the Vogtmann common split decomposition, computes the BHV cone path distance and the exact Owen Provan GTP geodesic side by side on matched and totally incompatible pairs, verifies the Pythagorean identity across orthogonal orthants, computes the Frechet mean of a sample of bootstrap trees via Bacak proximal iteration, shows tanglegram crossings under default and optimal layouts, builds a phylogenetic network and classifies it as tree based, runs Neighbor Net on a circular distance matrix, and finally computes the tropical Fermat Weber point of a set of dissimilarity vectors. Figures are written to `results/`.

## Testing

```bash
python3 -m pytest tests/ -q
```

Forty eight tests, all passing. Coverage includes split compatibility via Buneman, RF identity symmetry and positivity, quartet distance zero on identical trees, the Schroder double factorial count for n in three to six, NNI neighbor sanity and graph connectivity, BHV cone path equals GTP geodesic when sharing splits, path space formula on single pair and ascending and descending ratio cases, triangle inequality on five taxa triples, symmetry, Pythagorean identity across orthogonal orthants, Frechet mean zero on identical samples and correct averaging on weighted branches, Neighbor Net recovery of a circular ordering, tropical distance symmetry and shift invariance, Fermat Weber reduction on colinear points, and golden numerical tests pinned to the published Owen Provan worked examples.

## Dependencies

The core library requires `numpy`, `scipy` (for `scipy.optimize.linprog` in the tropical Fermat Weber LP), `networkx` (for the NNI rearrangement graph and shortest path), and `matplotlib` (for plotting). Install the full stack with

```bash
pip install numpy scipy networkx matplotlib pytest
```

## References

The literature review in `docs/lit_review.md` synthesizes Billera, Holmes, and Vogtmann 2001 (the geometry of phylogenetic tree space and CAT(0) property), Owen 2011 and Owen and Provan 2011 (the polynomial time geodesic algorithm with the path space enumeration via closure lattices), Miller, Owen, and Provan polyhedral averaging and Brown and Owen mean and variance papers (statistical summaries on CAT(0) spaces), Lin, Sturmfels, Tang, and Yoshida tropical geometry of phylogenetic tree space, Moon counting labeled trees (classical combinatorics), Huson and Bryant reviews of split networks and Neighbor Net, and a recent 2025 paper on tree rearrangement graphs admitting paths of decreasing Robinson Foulds distance. PDFs are included in the `literature/` folder.
