

# CFLG: an algorithmic toolkit for continuous set covering on networks

Facility location is an important appication in operations research.  Continuous set covering on networks generalizes the classical dicrete set covering problems in graphs, which allows both the demands and facilities continuously locating in edges.
CFLG contains several algorithms for solving continuous set covering on networks. CFLG is written in Julia based on [JuMP](https://jump.dev/JuMP.jl/stable/installation/).

## Updates
Newer version supports disjunctive programming and indicator formulations, new valid inequalities, and SCIP as solver option.
The package now is generated and managed by `pkg`.


## Prerequisite
CFLG has been developped and tested in *Linux* System. 

To use CFLG, it requires the following installation:
- Julia,
- JuMP,

and the following packages of MILP solvers:
- CPLEX.jl,
- GLPK.jl,
- Gurobi.jl,
- SCIP.jl.

It is necessary to have CPLEX.jl to build CFLG, the installation of CPLEX.jl requires a cplex binary executable. To use Gurobi, the system should have a licence installed.

## Installation
Clone the repository: 
```
git clone https://github.com/lidingxu/CFLG.git
cd CFLG
```
We recommend to build the package locally:
```
julia
julia> ]
pkg> activate .
(CFLG) pkg> build
```
Run a simple test:
```
(CFLG) pkg> test
```


## Benchmarks
There are six benchmarks `city`, `Kgroup_A`, `Kgroup_B`, `random_A`, `random_B` and `test`. Each benchmark contains instances in `.txt` file format.

Each instance file contains data of a graph. Its first line contains the number of nodes and edges of the graph; Edges are represented as a list from the second line: each line indicates an edge's end nodes and weight.  


## Usage

You can run `CFLG`  via the following command:
```
julia src/main.jl instance_dir instance_name  output_dir solver_name time_limit formulation cover
```

The settings of the above arguments are as follows.
  * `instance_dir`: directory of the instance.
  * `instance_name`:  the instance's name.
  * `output_dir`: the directory for the output result.
  * `solver_name`: one of `CPLEX`, `GLPK`, `Gurobi`, and `SCIP`.
  * `time_limit`:  time limit in seconds.
  * `formulation`: one of the following options:
  ```
    EF: edge model formulation, from "Covering edges in networks", Fröhlich et al.
    EFP: edge model big-M formulation with processing (bound tightenning and delimited cover)
    EFPV: edge model big-M formulation with processing (bound tightenning and delimited cover) and rank-1 valid inequalities
    EFPV2: edge model big-M formulation with processing (bound tightenning and delimited cover) and rank-2 valid inequalities   
    EFPD: edge model disjunctive programming formulation with processing (delimited cover)
    EFPI: edge model indicator constraint formulation with processing (delimited cover)
    EVF: edge-vertex model big-M formulation 
    EVFP0: edge-vertex model big-M formulation with simple processing (delimited cover)
    EVFP: edge-vertex model big-M formulation with processing (bound tightenning and delimited cover)
    LEVFP: long edge-vertex model big-M formulation with processing (bound tightenning and delimited cover)
    None: not a model, just record statistics of original graph, degree-2-free graph, subdivided (splitted) graph
  ```
  * `cover`: the covering raidus, `Small` or `Large`.

To reproduce the computational results in the accompanied paper, clear the `results` directory.
Edit the options in the `run.sh` file:
* `formulations`: an array of test formulations, for example, `("None" "EF" "EFP" "EFPV" "EFPV2" "EFPI" "EFPD" "EVFP" "LEVFP")`
* `covers`: an array of test coverage, for example, `("Small" "Large")`
* `solver`: one of the solver "Gurobi", "CPLEX", "GLPK", "SCIP"
* `gnuparalleltest=1`: enable the GNU parallel test option `gnuparalleltest` according to their availability in your system. # enable: 1, disable: 0
* `juliabin`: the path of julia binary executable

Execute the following command to run all tests:
```
/bin/bash run.sh 
```

Run the following codes to produce the figures and tables in the `experiment` directory:
```
/bin/bash parse_result.sh
```


## References

If you find CFLG useful in your work, we kindly request that you cite the following papers, which are recommended reading for advanced users:

    @misc{xu2024formulations,
        title={Formulations of the continuous set-covering problem on networks: a comparative study},
        author={Liding Xu and Claudia D'Ambrosio},
        year={2024},
        eprint={2404.10467},
        archivePrefix={arXiv}
    }

    @article{PELEGRIN2023102835,
      title = {Continuous covering on networks: Improved mixed integer programming formulations},
      journal = {Omega},
      volume = {117},
      pages = {102835},
      year = {2023},
      issn = {0305-0483},
      doi = {https://doi.org/10.1016/j.omega.2023.102835},
      url = {https://www.sciencedirect.com/science/article/pii/S0305048323000014},
      author = {Mercedes Pelegrín and Liding Xu},
    }


