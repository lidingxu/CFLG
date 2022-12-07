

# CFLG: an algorithmic toolkit for continuous set covering on networks

Facility location is an important appication in operations research.  Continuous set covering on networks generalizes the classical dicrete set covering problems in graphs, which allows both the demands and facilities continuously locating in edges.
CFLG contains several algorithms for solving continuous set covering on networks. CFLG is written in Julia based on [JuMP](https://jump.dev/JuMP.jl/stable/installation/).

## Updates
The version 2.0.0 supports disjunctive programming formulation, new valid inequalities, and an implementation in SCIP.
The package now is generated and managed by `pkg`.


## Installation
CFLG has been developped and tested in *Linux* System. 

To use CFLG, it requires:
- Julia,
- JuMP,

and one of the following MILP solvers:
- CPLEX,
- GLPK,
- Gurobi,
- SCIP.




## Benchmarks
There are six benchmarks `city`, `Kgroup_A`, `Kgroup_B`, `random_A`, `random_B` and `test`. Each benchmark contains instances in `.txt` file format.

Each instance file contains data of a graph. Its first line contains the number of nodes and edges of the graph; Edges are represented as a list from the second line: each line indicates an edge's end nodes and weight.  

## Usage

You can run `CFLG.jl`  via the following command:
```
julia ./CFLG.jl solver time_limit instance algorithm raidus
```

The settings of the above arguments are as follows.
  * `solver`: one of the solvers in `CPLEX`, `GLPK` and `Gurobi`.
  * `time_limit`:  time limit in seconds.
  * `instance`:  the instance path.
  * `algorithm`: one of the algorithms in `EF`, `F0`, `F`, `SF`, `RF`, `SFD` and `None`.
  * `raidus`: the coverage raidus.

To reproduce the computational results in the accompanied paper, clear the `results` directory, go to the `test` directory, and execute the following command in the terminal (in *Linux*)
```
/bin/bash run.sh
```
Note that you should change the solver option `solver` and the GNU parallel test option `gnuparalleltest` according to their availability in your system.


### Algorithm Option
The algorithm option is summarized in the following table.


|     | Continuous demand |  Continuous facility |   Set delimitation  | Strengthenning| Long edge| Model size | Input graph     | Comment|  
|-----|-------------------|----------------------|---------------------|---------------|----------|------------|-----------------|--------|
| `EF`|      Yes          | Yes                  | No                  |   No          | No       | Very large | Subdivided graph|From [Covering edges in networks](https://onlinelibrary.wiley.com/doi/full/10.1002/net.21924) 
| `F0`|      Yes          | Yes                  | No                  |   No          | No       | Large     | Subdivided graph| Naive model
| `F` |      Yes          | Yes                  | Yes                 |   No          | No       | Medium     | Subdivided graph|Complete model
| `SF`|      Yes          | Yes                  | Yes                 |   Yes         | No       | Meidum     | Subdivided graph|Strenghtenned model
| `RF`|      Yes          | Yes                  | Yes                 |   Yes         | Yes      | Small      | Degree-2-free graph| Reduced model
| `SFD`|      No         | Yes                  | Yes                 |   Yes         | No       | Very small     | Subdivided graph|Discrete model

The `None` option is used for computing the parameters of subdivided and degree-2-free graphs.



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


