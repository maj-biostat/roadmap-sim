
cleantmp:
  find tmp -delete
  mkdir tmp


sim09 cfg:
  Rscript --vanilla ./R/data-sim09.R sim09_sim_loop {{cfg}} 
  just cleantmp

runsim09:
  just sim09 ../etc/sim09/cfg-sim09-sc01-v01.yml
  
  