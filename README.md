# Generation Expansion Planning (GEP) models considering uncertainties on renewable energy resources (RES) using python/pyomo

The following files solve the GEP problem for three scenarios of wind and solar production using different approaches:

* **Stochastic-GEP.py**: Two-Stage Stochastic Generation Expansion Planning

The models are developed in [GAMS](https://www.gams.com/) and solved with [CPLEX](https://www.ibm.com/analytics/cplex-optimizer), but you could use any other solver (e.g., [GUROBI](https://www.gurobi.com/), [Cbc](https://github.com/coin-or/Cbc)).

The main references to model the optimization problems are:

[1] [Optimization Techniques by Andrés Ramos Galán](https://pascua.iit.comillas.edu/aramos/OT.htm)
[2] [Good Optimization Modeling Practices using Pyomo](https://pascua.iit.comillas.edu/aramos/simio/transpa/s_GoodOptimizationModelingPracticesPyomo.pdf)
