# Project Plan -- Phylogenetic Networks & Tree Space

## Goal
Build a research-quality library that exercises **algebra + analysis +
combinatorics + graph theory** in concert, using phylogenetic trees and
networks as the domain. End state: compute distances, means, and
rearrangement paths between trees, and explore a restricted class of
phylogenetic networks -- all from mathematical first principles.

## Layers (each layer = increasing math depth)

### Layer 0 -- Foundations (combinatorics + graph theory)
- Tree data structure with parent/child pointers.
- **Split** representation as frozenset of taxa + bitmask.
- Newick I/O (ASCII parser + writer).
- Enumerate all (2n?3)!! rooted binary trees (small n).
- Verify Cayley / Schroder counts.

### Layer 1 -- Tree metrics (combinatorics)
- **Robinson-Foulds** distance via split symmetric difference.
- Weighted RF (?_1 and ?_2 on split-length vectors).
- **Quartet distance** (combinatorial, O(n^2)).
- **NNI neighborhood** enumeration; build an NNI-adjacency graph
  (small n) and compute shortest-path distances on it.

### Layer 2 -- BHV tree space (analysis + topology)
- Represent points as (topology, length-vector); the topology is a
  compatible split set; length-vector is in the orthant R_{>=0}^{n?2}.
- Implement the **link of origin L_n**: simplicial complex of pairwise-
  compatible split collections. Verify flag-complex property.
- Implement **cone path** distance d_cone(T_1, T_2) = |T_1| + |T_2|.
- Implement **common-split decomposition** (Vogtmann thm) reducing
  geodesic to the "no-common-splits" problem.
- Implement the **GTP algorithm** (Owen-Provan) on small instances:
  - build incompatibility poset P(Sigma_1, Sigma_2),
  - enumerate maximal chains in the closure lattice K,
  - for each candidate path space, solve the Euclidean touring problem
    in linear time,
  - take the minimum.
- Unit-tests: CAT(0) triangle inequality, cone path upper-bounds geodesic,
  monotonicity under common-split decomposition.

### Layer 3 -- Statistics on BHV space (analysis)
- **Frechet mean**: Bacak proximal iteration
  p_{k+1} = Geodesic(p_k, T_{(k mod N)+1}, 1/(k+1))
  -- exploits unique BHV geodesics.
- **Sample Frechet variance**; bootstrap confidence.
- **Principal geodesics** (Nye): line-of-best-fit minimizing sum of
  squared projection distances; solved via alternating optimization.

### Layer 4 -- Combinatorial enumerations (graph theory)
- Count tanglegram crossings (OTCM via sorting-by-LCAs; O(n log n)).
- Generate NNI / SPR rearrangement graphs for small n; confirm
  diameters against published bounds.
- Compute SPR distance via ILP / BFS for small trees.

### Layer 5 -- Phylogenetic networks (graph theory + algebra)
- Rooted DAG data structure with reticulation bookkeeping.
- Implement **Neighbor-Net** to build a circular split system from a
  distance matrix (this is pure combinatorial linear algebra: the
  non-negative LSQ for edge weights is a convex QP).
- Draw planar split network for circular systems.
- **Tree-child network** construction from multiple gene trees
  (basic FPT algorithm).
- Distinguish tree-based vs. non-tree-based networks on small inputs.

### Layer 6 -- Tropical geometry (algebra + analysis)
- Implement **tropical distance** d_inf on ?_inf quotient space.
- Map BHV tree -> dissimilarity vector (via four-point embedding).
- Compute **tropical Fermat-Weber point** via LP.
- Tropical convex hull; tropical PCA (best-fit tropical line).

## Milestones

| # | Deliverable | Layers | Time budget |
|---|-------------|--------|-------------|
| M1 | Tree class + Newick + split enumeration + RF + quartet | 0-1 | 1 session |
| M2 | Link-of-origin + cone path + Vogtmann decomposition | 2 | 1 session |
| M3 | GTP geodesic algorithm (complete Owen-Provan) | 2 | 2 sessions |
| M4 | Frechet mean + variance on real bootstrap data | 3 | 1 session |
| M5 | Tanglegram crossing + rearrangement graphs | 4 | 1 session |
| M6 | Neighbor-Net + network DAG + tree-child networks | 5 | 2 sessions |
| M7 | Tropical distance, FW point, tropical PCA | 6 | 1 session |
| M8 | Case study: vertebrate mitochondrial or HIV dataset | -- | 1 session |

## Data sources
- `dendropy` / `ete3` libs for Newick tests (not computation).
- Small real dataset: primate mitochondrial 12s rRNA (12 taxa, 1 kb),
  available via the TreeBASE study S2018.
- Simulated: generate gene trees under coalescent (msprime) for
  benchmarking network inference.

## Stack
- Python 3.10+ (dataclasses, typing).
- NumPy / SciPy (matrix tree theorem, Laplacians, LP).
- NetworkX (graph utilities for rearrangement graphs, network DAGs).
- matplotlib (visualization).
- pytest (unit tests).
- No heavy phylo libs in the core path -- we implement algorithms
  ourselves. Import dendropy **only** to validate parsing.

## Success criteria
- BHV geodesic results agree with published examples (Owen's paper
  Figures 4-5).
- Frechet mean converges to the true tree on simulated bootstrap trees.
- NNI graph for n=5 has the right vertex/edge counts (15 vertices,
  15 edges per published combinatorics).
- Neighbor-Net on a circular distance matrix recovers the exact
  network.
