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

