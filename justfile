
cleantmp:
  find tmp -delete
  mkdir tmp

sim07 cfg:
  Rscript --vanilla ./R/data-sim07-all.R sim07_sim_loop {{cfg}} 
  just cleantmp
  
runsim07:
  just sim07 ../etc/sim07/cfg-sim07-sc01-v01.yml
  just sim07 ../etc/sim07/cfg-sim07-sc01-v02.yml

sim09 cfg:
  Rscript --vanilla ./R/data-sim09.R sim09_sim_loop {{cfg}} 
  just cleantmp
  
runsim09:
  just sim09 ../etc/sim09/cfg-sim09-sc01-v01.yml
  #just sim09 ../etc/sim09/cfg-sim09-sc01-v02.yml
  just sim09 ../etc/sim09/cfg-sim09-sc01-v03.yml
  #just sim09 ../etc/sim09/cfg-sim09-sc01-v04.yml
  #just sim09 ../etc/sim09/cfg-sim09-sc01-v05.yml
  just sim09 ../etc/sim09/cfg-sim09-sc01-v06.yml
  #just sim09 ../etc/sim09/cfg-sim09-sc01-v07.yml
  just sim09 ../etc/sim09/cfg-sim09-sc01-v08.yml
  #just sim09 ../etc/sim09/cfg-sim09-sc01-v09.yml
  just sim09 ../etc/sim09/cfg-sim09-sc01-v10.yml
  just sim09 ../etc/sim09/cfg-sim09-sc01-v11.yml
  #just sim09 ../etc/sim09/cfg-sim09-sc01-v12.yml
  #just sim09 ../etc/sim09/cfg-sim09-sc01-v13.yml
  #just sim09 ../etc/sim09/cfg-sim09-sc01-v14.yml
  