# ROADMAP Simulations

Simulation implementation and simulation results for the trial operating characteristics.

To build this site, you need to have a recent version of quarto install plus there is a dependency on the `roadmap.data` R package, see [https://github.com/maj-biostat/roadmap.data].


## Running simulations

Simulations are configured via Makefile. 
Run `make list` to see all targets. 
At the time of writing `make sim01` and `make sim02` are the primary simulations.
The former runs a single analysis and summarises the power at different treatment effect sizes.
The range of effects are from log(1/2) to log(2) in all arms.

