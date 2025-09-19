

# cflg: an algorithmic toolkit for continuous set covering on networks

Facility location is an important application in operations research.  Continuous set covering on networks generalizes the classical dicrete set covering problems in graphs, which allows both the demands and facilities continuously locating in edges.
cflg contains several algorithms for solving continuous set covering on networks. cflg is written in Julia based on [JuMP](https://jump.dev/JuMP.jl/stable/installation/).


## Installation
cflg has been developped and tested under *Linux* System.

To use cflg, it requires:
- Julia,
- JuMP,

and one of the following MILP solvers:
- CPLEX,
- GLPK,
- Gurobi.


## Benchmarks
There are six benchmarks `city`, `Kgroup_A`, `Kgroup_B`, `random_A`, `random_B`, `tree_A`, `tree_B` and `test`. Each benchmark contains instances in `.txt` file format.

Each instance file contains data of a graph. Its first line contains the number of nodes and edges of the graph; Edges are represented as a list from the second line: each line indicates an edge's end nodes and weight.

## Usage

You can run `cflg.jl` via the following command as the following example:
```
julia src/main.jl data_folder instance result_folder solver time_limit algorithm raidus
```


The algorithm option names and short descriptions (features are inferred from the bitmask flags):


| Algorithm | Model | Edge vs Edge-Vertex | Delimitation (P0) | Bound tighten (P1) | Valid Ineq. (V) | Disjunctive (D) | Long edge (L) | Indicator (I) | Attached pruning (A) | Comment |
|---|:---:|:-------------------:|:-----------------:|:------------------:|:---------------:|:-------:|:-----------:|:-------------:|:--------------------:|---------|
| `EF`     | Yes | Edge        | No  | No  | No  | No  | No  | No  | No  | Edge model (Fröhlich et al.) |
| `EFP0`   | Yes | Edge        | Yes | No  | No  | No  | No  | No  | No  | Edge model + delimitation (simple processing) |
| `EFP`    | Yes | Edge        | Yes | Yes | No  | No  | No  | No  | No  | Edge big‑M + delimitation + bound tightening |
| `EVF`    | Yes | Edge‑Vertex | No  | No  | No  | No  | No  | No  | No  | Edge‑vertex big‑M formulation |
| `EVFP`   | Yes | Edge‑Vertex | Yes | Yes | No  | No  | No  | No  | No  | EVF + processing (bound tightening + delimitation) |
| `EVFPV`  | Yes | Edge‑Vertex | Yes | Yes | Yes | No  | No  | No  | No  | EVFP + simple valid inequalities |
| `LEVFP`  | Yes | Edge‑Vertex | Yes | Yes | No  | No  | Yes | No  | No  | Long edge‑vertex big‑M + processing |
| `LEFPI`  | Yes | Edge        | Yes | Yes | No  | No  | Yes | Yes | No  | Edge model with long‑edge + indicator |
| `LEFP`   | Yes | Edge        | Yes | Yes | No  | No  | Yes | No  | No  | Long‑edge edge model + processing |
| `LEFPD`  | Yes | Edge        | Yes | Yes | No  | Yes | Yes | No  | No  | LEFP + disjunctive formulation |
| `LEFPV`  | Yes | Edge        | Yes | Yes | Yes | No  | Yes | No  | No  | LEFP + valid inequalities |
| `None`   | No  | —           | No  | No  | No  | No  | No  | No  | No  | Not a model — compute statistics (original, degree‑2‑free, subdivided) |

Use the exact algorithm name (e.g., `LEFPB`, `None`) when invoking the program.




## References

If you find cflg useful in your work, we kindly request that you cite the following paper draft ([arXiv preprint](https://arxiv.org/abs/2203.00284)), which is recommended reading for advanced users:

    @misc{pelegrín2022continuous,
        title={Continuous Covering on Networks: Strong Mixed Integer Programming Formulations},
        author={Mercedes Pelegrín and Liding Xu},
        year={2022},
        eprint={2203.00284},
        archivePrefix={arXiv},
        primaryClass={math.OC}
    }

