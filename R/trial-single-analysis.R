source("./R/init.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_single_analysis"
} else {
  log_info("Run method ", args[1])
}


idx_sim = 1
trial_single_analysis <- function(idx_sim = 1){
  
  # Get all data associated with this trial index
  db <- dbConnect(RSQLite::SQLite(), g_cfg$db_name)
  d_qry <- dbSendQuery(db, "SELECT * FROM trial_data WHERE sim = ?")
  dbBind(d_qry, list(idx_sim))
  d_all <- data.table(dbFetch(d_qry))
  dbClearResult(d_qry)
  dbDisconnect(db)
  

  
  
  
  
  
  
  
}