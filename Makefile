#
# The general format of a rule is:
# target: pre_req_i
# tab->recipe
#

# Definitions
DAT := ./data
OUT := ./site
LOG := ./logs

# R scripts
R_SCRIPTS = $(wildcard ./R/*.R)
R_OPTS=--vanilla

# PASSWORD ?= $(shell bash -c 'read -s -p "Password: " pwd')
# PASSWORD ?= $(shell stty -echo; read -p "Password: " pwd; stty echo; echo $$pwd)
# PASSWORD ?= $(shell bash -c 'read -s -p "Password: " pwd; echo $$pwd')


.PHONY: list
list:
	@LC_ALL=C $(MAKE) -pRrq -f $(firstword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/(^|\n)# Files(\n|$$)/,/(^|\n)# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'

# Note: the following PHONY simply allows me to say "make sim_dat" rather than
# having to spell out the full target name. It is a way to create abbreviations.
# Note: the $@ is an auto variable that is substituted with the target name so
# in this case $@ is set to $(DAT_SIM)/dat_sim.qs, which is the output we want.

# Removes all derived files
.PHONY: clean
clean:
	rm -fR logs/*
	rm -f data/*.sqlite
	rm -fr _site/*

# @echo "\nPassword is" $(PASSWORD)
.PHONY: preview
preview:
	quarto preview

.PHONY: render
render:
	quarto render

.PHONY: publish
publish:
	quarto publish gh-pages

.PHONY: data
data:
	Rscript $(R_OPTS) ./R/data.R run_data

.PHONY: sim01
sim01:
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-01.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-02.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-03.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-04.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-05.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-06.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-07.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-08.yml


.PHONY: sim02
sim02:
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc02-01.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc02-02.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc02-03.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc02-04.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc02-05.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc02-06.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc02-07.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc02-08.yml


.PHONY: sim03
sim03:
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc03-01.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc03-02.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc03-03.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc03-04.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc03-05.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc03-06.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc03-07.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc03-08.yml


