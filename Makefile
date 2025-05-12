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
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-v01.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-v02.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-v03.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-v04.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-v05.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-v06.yml
	Rscript $(R_OPTS) ./R/sim-01.R run_sim_01 cfg-sim01-sc01-v07.yml

# sequential design sim-02.R 
.PHONY: sim02
sim02:
	Rscript $(R_OPTS) ./R/sim-02.R run_sim_02 cfg-sim02-sc01-v01.yml
	Rscript $(R_OPTS) ./R/sim-02.R run_sim_02 cfg-sim02-sc01-v02.yml
	Rscript $(R_OPTS) ./R/sim-02.R run_sim_02 cfg-sim02-sc01-v03.yml
	Rscript $(R_OPTS) ./R/sim-02.R run_sim_02 cfg-sim02-sc01-v04.yml
	Rscript $(R_OPTS) ./R/sim-02.R run_sim_02 cfg-sim02-sc01-v05.yml

# sequential design sim-02.R 
.PHONY: sim03
sim03:
	Rscript $(R_OPTS) ./R/sim-03.R run_sim_03 cfg-sim03-sc01-v01.yml
	Rscript $(R_OPTS) ./R/sim-03.R run_sim_03 cfg-sim03-sc01-v02.yml
	Rscript $(R_OPTS) ./R/sim-03.R run_sim_03 cfg-sim03-sc01-v03.yml
	Rscript $(R_OPTS) ./R/sim-03.R run_sim_03 cfg-sim03-sc01-v04.yml
	Rscript $(R_OPTS) ./R/sim-03.R run_sim_03 cfg-sim03-sc01-v05.yml
	Rscript $(R_OPTS) ./R/sim-03.R run_sim_03 cfg-sim03-sc01-v06.yml
	Rscript $(R_OPTS) ./R/sim-03.R run_sim_03 cfg-sim03-sc01-v07.yml

# sequential design sim-02.R 
.PHONY: sim04
sim04:
	Rscript $(R_OPTS) ./R/sim-04.R run_sim_04 cfg-sim04-sc01-v01.yml
	Rscript $(R_OPTS) ./R/sim-04.R run_sim_04 cfg-sim04-sc01-v02.yml
	Rscript $(R_OPTS) ./R/sim-04.R run_sim_04 cfg-sim04-sc01-v03.yml
	Rscript $(R_OPTS) ./R/sim-04.R run_sim_04 cfg-sim04-sc01-v04.yml
	Rscript $(R_OPTS) ./R/sim-04.R run_sim_04 cfg-sim04-sc01-v05.yml
	Rscript $(R_OPTS) ./R/sim-04.R run_sim_04 cfg-sim04-sc01-v06.yml
	Rscript $(R_OPTS) ./R/sim-04.R run_sim_04 cfg-sim04-sc01-v07.yml

# 30/Jan/2025
# sequential design sim-05.R 
.PHONY: sim05
sim05:
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v01.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v02.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v03.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v04.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v05.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v06.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v07.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v08.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v09.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v10.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v11.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v12.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v13.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v14.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v15.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-05.R run_sim_05 ./sim05/cfg-sim05-sc01-v16.yml
	sleep 10
	rm tmp/*.csv

# 9/May/2025
# sequential design sim-06.R 
.PHONY: sim06
sim06:
	Rscript $(R_OPTS) ./R/sim-06.R run_sim_06 ./sim06/cfg-sim06-sc01-v01.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-06.R run_sim_06 ./sim06/cfg-sim06-sc01-v02.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-06.R run_sim_06 ./sim06/cfg-sim06-sc01-v03.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-06.R run_sim_06 ./sim06/cfg-sim06-sc01-v04.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-06.R run_sim_06 ./sim06/cfg-sim06-sc01-v05.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-06.R run_sim_06 ./sim06/cfg-sim06-sc01-v06.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-06.R run_sim_06 ./sim06/cfg-sim06-sc01-v07.yml
	sleep 10
	rm tmp/*.csv
	Rscript $(R_OPTS) ./R/sim-06.R run_sim_06 ./sim06/cfg-sim06-sc01-v08.yml
	sleep 10
	rm tmp/*.csv
	



