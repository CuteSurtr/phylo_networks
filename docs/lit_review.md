# Literature Synthesis -- Project 4: Phylogenetic Networks & Tree Space

## 1. Core foundations

### 1.1 Trees, splits, and compatibility
A phylogenetic tree on label set X = {1,...,n} (with root 0) is a rooted
semi-labeled tree whose leaves bijective with X and whose interior edges
each induce a **split** A|B, a bipartition of X ? {0} with |A|, |B| >= 2.

**Compatibility (Buneman):** two splits e = X|X', e' = Y|Y' are compatible
iff one of the four intersections X?Y, X?Y', X'?Y, X'?Y' is empty.
A set of pairwise-compatible splits encodes a unique tree topology
(Semple-Steel thm 3.1.4; Buneman).

### 1.2 Counting
- Binary rooted semi-labeled trees on n leaves: (2n?3)!! = (2n?2)!/(2^(n?1)(n?1)!)
- Unlabeled binary rooted trees: Catalan C_{n?1}
- Labeled trees (Cayley): n^(n?2) spanning trees on K_n.
- Link of origin L_n: simplicial complex of (n?3)-dim with one k-simplex
  per tree having k+1 interior edges; has the homotopy type of a wedge of
  (n?1)! spheres of dim (n?3) (Vogtmann 1990).

### 1.3 Tree rearrangements
- **NNI** (nearest-neighbor interchange): collapse an interior edge to
  degree-4 vertex, re-expand in one of the other two ways. Max distance
  O(n log n) (Sleator et al.), 2n?6 within a single associahedron.
- **SPR** (subtree prune + regraft): prune a subtree, regraft on any
  other edge. SPR distance is NP-hard.
- **TBR** (tree bisection + reconnection): generalizes SPR. Also NP-hard.
- **Robinson-Foulds (RF)** distance: |Sigma(T_1) ? Sigma(T_2)|, the symmetric
  difference between split sets. Computable in O(n) via Day's algorithm.

## 2. BHV tree space (Billera-Holmes-Vogtmann 2001)

### 2.1 Construction
T_n is built by taking (2n?3)!! Euclidean orthants of dimension (n?2),
one per binary topology. A point in orthant O(tau) with coords (l_1,...,l_{n?2})
is the tree with topology tau and interior edge lengths l_i >= 0.
Orthants are glued along boundary faces corresponding to shared splits
(setting an edge length to 0 collapses that split).

### 2.2 Metric and geometry
- Distance = length of shortest segmented path.
- **Theorem 4.1 (BHV):** T_n is a **CAT(0)** space. Proof: subdivide each
  orthant into unit cubes; then T_n is a cubical complex; the link of
  every vertex is a flag simplicial complex (pairwise-compatibility =>
  full compatibility). Gromov's theorem => CAT(0).
- Consequence: unique geodesic between any two points; well-defined
  Frechet means, centroids, convex hulls.
- Cone path: T -> 0 -> T' is a geodesic iff the angle at origin >= pi.

### 2.3 Tropical / algebraic connections
T_n is homeomorphic (not isometric) to the tropical Grassmannian
Gr(2,n) and to the Bergman fan of the graphic matroid of K_n
(Speyer-Sturmfels). Edge-length coordinates ? distance coordinates
via the four-point condition.

## 3. Owen-Provan geodesic algorithm (2011)

### 3.1 Key reductions
1. **Vogtmann decomposition:** if T_1, T_2 share split e,
   d(T_1,T_2)^2 = d(T_1^X,T_2^X)^2 + d(T_1^Y,T_2^Y)^2 + (|e|_{T_1} ? |e|_{T_2})^2.
   Recurse until no shared splits (the **essential problem**).
2. For the essential problem: enumerate candidate orthant sequences via
   the **path poset** K(Sigma_1, Sigma_2).

### 3.2 Combinatorial structure
- **Incompatibility poset** P(Sigma_1, Sigma_2): equivalence classes of Sigma_2
  under identical crossing set X_{Sigma_1}(*), ordered by inclusion.
- **Path poset** K(Sigma_1, Sigma_2): closed sets of Sigma_2 under the closure
  operator A ? \bar A = {f : X(f) ? X(A)}, ordered by inclusion.
- Maximal chains in K ? maximal path spaces (Theorem 3.7).
- Remark 3.8: |K| can be exponential.

### 3.3 Euclidean reduction
For a fixed candidate path space, the path length minimization reduces
to a **touring problem with obstacles** in R^k:  shortest path from
positive orthant to negative orthant through a sequence of orthants
with i coords <= 0 and rest >= 0 (Theorem 4.4). Solvable in O(k) per
path space (Theorems 4.10, 4.11).

### 3.4 Dynamic programming + divide/conquer
Combine via Theorem 5.2 -- the optimal geodesic between T_1 and T_2
is related to geodesics between simpler subtrees, giving practical
algorithms (O(n^4) expected; polynomial-time GTP by Owen 2011).

## 4. Statistics on BHV tree space

- **Frechet mean** mu = argmin_p Sigma d(p, T_i)^2 exists uniquely on any
  CAT(0) space (Sturm 2003).
- Iterative algorithms: Bacak (2014) proximal point; Sturm's law of
  large numbers for "random walks" on CAT(0) spaces.
- **Miller-Owen-Provan (2015):** polyhedral computational geometry gives
  exact Frechet mean via locally-Euclidean approximations.
- **Brown-Owen (2020):** Frechet variance via geodesic summation;
  applied to bootstrap-like confidence regions on trees.
- **Barden-Le CLT (2018):** asymptotic normality of sample means on
  T_n (stickiness at lower strata complicates limit).
- **Nye et al. PCA:** principal paths / principal geodesics in BHV space.

## 5. Phylogenetic networks (reticulation)

### 5.1 Model
A rooted phylogenetic network on X is a rooted DAG where:
- leaves are bijective with X,
- every node has in-degree + out-degree >= 2 except leaves/root,
- **reticulation nodes** have in-degree >= 2 (hybridization, HGT,
  recombination).

### 5.2 Classes (progressively restrictive)
- **General networks** -> hard to infer, NP-hard reticulation min.
- **Tree-child**: every non-leaf has at least one tree-child descendant.
- **Tree-based**: contains a spanning base tree.
- **Normal / level-k**: bounded reticulation per biconnected component.
- **Galled trees / galled networks**: reticulations isolated in cycles.

### 5.3 Encodings
- **Split networks** (Huson-Bryant): planar representation of
  incompatible splits; **Neighbor-Net** builds circular split systems
  via agglomerative clustering on pairwise distances -> always planar.
- **Reticulation networks**: DAG-based explicit models.
- Metrics on networks: tripartition distance, mu-distance, tree-based
  distances; most NP-hard to compute.

### 5.4 Rearrangement moves
- **Tail moves** and **head moves** generalize SPR to networks
  (Bordewich-Semple-Tokes).
- Connect rooted network space under these moves; combinatorial
  diameter / connectedness questions active.

## 6. Tropical geometry of tree space

### 6.1 Tropical metric
"Palm tree space" uses the tropical (Hilbert projective) metric
d_inf(u,v) = max_ij((u_i ? v_i) ? (u_j ? v_j)) on the quotient
R^{(n choose 2)} / R*(1,...,1). Equivalent to the ?_inf metric up to
translation.

### 6.2 Why it matters
- Unlike BHV, tropical statistics admit tractable **tropical Fermat-
  Weber points**, tropical PCA, tropical convexity (via max-plus
  polytopes).
- Lin, Sturmfels, Tang, Yoshida: tropical regression, tropical
  logistic regression, tropical SVM on tree space.
- Computational advantage: tropical convex hull in polynomial time
  (vs. BHV exact convex hull which is open).

## 7. Open / advanced topics touched by the downloads

- **BHV with non-identical leaves:** Ren et al. 2017 (connecting BHV
  spaces of different taxa); 2024 extension-space distances.
- **Stickiness / finite-sample CLT:** statistics on positively-curved
  strata.
- **Robinson-Foulds networks:** generalizations of RF to reticulation.
- **Information geometry:** Riemannian structure on pdfs of tree
  models (Garba et al., JMB 2021).

---

## Reading-list summary (files in `literature/`)

| Paper | Role |
|-------|------|
| BHV_2001_geometry_tree_space | foundational construction of T_n |
| Owen_Provan_arxiv_v2 | polynomial geodesic algorithm |
| Owen_Provan_2011_geodesic_algorithm | fast variant |
| Brown_Owen_mean_variance_trees | Frechet mean/variance |
| Miller_et_al_polyhedral_averaging | exact means via polyhedra |
| Lin_tropical_geometry_treespace | tropical statistics framework |
| uniqueness_intractability_networks | network reconstruction complexity |
| connecting_BHV_spaces_different_taxa | cross-taxa BHV geometry |
| tree_rearrangement_RF_2025 | rearrangement graphs + RF |
| Allman_Rhodes_phylogenetics_lectures | textbook background |
| Moon_counting_labelled_trees | classical enumeration |
| matrix_tree_theorem_Cayley | Cayley/Kirchhoff proofs |
