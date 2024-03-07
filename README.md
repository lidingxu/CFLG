

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

and at least one package of the following MILP solvers:
- CPLEX.jl,
- GLPK.jl,
- Gurobi.jl,
- SCIP.jl.


## Installation
The installation is manual. 
```
git clone https://github.com/lidingxu/CFLG.git
cd CFLG
```


## Benchmarks
There are six benchmarks `city`, `Kgroup_A`, `Kgroup_B`, `random_A`, `random_B` and `test`. Each benchmark contains instances in `.txt` file format.

Each instance file contains data of a graph. Its first line contains the number of nodes and edges of the graph; Edges are represented as a list from the second line: each line indicates an edge's end nodes and weight.  


## Usage

You can run `CFLG.jl`  via the following command:
```
julia ./CFLG.jl solver time_limit instance algorithm raidus
```

The settings of the above arguments are as follows.
  * `solver`: one of the solvers in `CPLEX`, `GLPK`, `Gurobi`， and `SCIP`.
  * `time_limit`:  time limit in seconds.
  * `instance`:  the instance path.
  * `algorithm`: one of the algorithms in ``.
  * `raidus`: the coverage raidus.

To reproduce the computational results in the accompanied paper, clear the `results` directory, go to the `test` directory, and execute the following command in the terminal (in *Linux*)
```
/bin/bash run.sh
```
Note that you should change the solver option `solver` and the GNU parallel test option `gnuparalleltest` according to their availability in your system.



## References

If you find CFLG useful in your work, we kindly request that you cite the following paper draft ([arXiv preprint](https://arxiv.org/abs/2203.00284)), which is recommended reading for advanced users:

    @misc{pelegrín2022continuous,
        title={Continuous Covering on Networks: Strong Mixed Integer Programming Formulations}, 
        author={Mercedes Pelegrín and Liding Xu},
        year={2022},
        eprint={2203.00284},
        archivePrefix={arXiv},
        primaryClass={math.OC}
    }


