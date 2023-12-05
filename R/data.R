source("./R/init.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_data"
} else {
  log_info("Run method ", args[1])
}


g_silo = c("early_acute", "late_acute" , "chronic") 
g_pr_silo <- c(0.3, 0.5, 0.2)
names(g_pr_silo) <- g_silo

g_joint <- c("knee", "hip")
g_pr_joint <- c(0.6, 0.4)
names(g_pr_joint) <- g_joint

g_fct_abty <- c("rif", "soc")
g_fct_la_surg <- c("dair", "one", "two")
g_fct_ch_surg <- c("one", "two")
g_fct_abdu <- c("shrt", "long")


# Trial config specifies the setup for the sim
# Hard boiled list for now. 
# Will make more flexible down the line.
get_trial_cfg <- function(){
  
  # All combinations
  
  X <- CJ(
    silo = c("early_acute", "late_acute", "chronic"),
    surg = c("dair", "one", "two"),
    abdu = c("shrt", "long"),
    abty = c("rif", "soc")
  )
  setkey(X, silo, surg, abdu, abty)
  # Allowable set of combinations
  X <- rbind(
    X[silo == "early_acute" & surg == "dair" & abdu == "long"],
    X[silo == "late_acute" & surg == "dair" & abdu == "long"],
    X[silo == "late_acute" & surg %in% c("one", "two"), ],
    X[silo == "chronic" & surg %in% c("one", "two"), ]
  )
  X <- cbind(id_x = 1:nrow(X), X)

  # Todo - 
  # make up probability of response for each group.
  # function to compute risk diffs,  risk ratios, odds ratios
  
  l <- list(
    N_sim = 1000,
    N = 2500,
    X = X,
    pr_y = rep(0.7, nrow(X)),
    pr_la_surg = c(dair = 0.5, one = 0.25, two = 0.25),
    # 50pct chance of being gram positive
    pr_gpos = 0.5,
    # enrolment times - poisson process pt per day
    nhpoislambda = 2500/(365*4),
    # ramp up over 6 months (180 days)
    nhomogintens = function(t) pmin(t/180, 1)
  )
  
  l
  
}

g_trial_cfg <- get_trial_cfg()


#



get_data_demo_1 <- function(
    N = 10000,
    a0 = qlogis(0.55),
    # late acute does worse
    a_la = -0.5,
    # rif does nothing.
    c_abty = 0,
    # hip joint infection decreases your chances of success
    z_joint = -0.4,
    pr_ea_joint = 0.6,
    pr_la_joint = 0.35
    ){
  
  d <- data.table()
  d[, pt := 1:N]
  
  # Stratification
  silo_sub <- g_silo[1:2]
  pr_silo_sub <- g_pr_silo[1:2]/sum(g_pr_silo[1:2])
  
  silo <- factor(sample(silo_sub, size = N, replace = T, prob = pr_silo_sub), levels = silo_sub)
  d[, silo := silo]
  
  # Balanced rand for abty across all silo
  d[, abty := sample(g_fct_abty, .N, replace = T)]
  
  # Just for illustration - binary 1 => hip
  # early acute tend to have more hip infections
  d[silo == "early_acute", joint := rbinom(.N, 1, pr_ea_joint)]
  d[silo == "late_acute", joint := rbinom(.N, 1, pr_la_joint)]
  
  setcolorder(d, c("pt", "silo", "joint"))
  
  # Outcome model
  
  # For EA outcome model
  # eta = a_ea + c_abty I(abty == "rif") + z_joint I(joint == 1)
  # eta = a_la + c_abty I(abty == "rif") + z_joint I(joint == 1)
  
  # Treatment success at 12 months
  
  d[, pr_y := plogis( 
    a0 + a_la  * (silo == "late_acute") + c_abty * (abty == "rif") + z_joint * joint)]

  d[, y := rbinom(.N, 1, prob = pr_y)]
  
  d
}

N = 10000
# baseline log-odds of trt success by silo
b_silo = c(early_acute=qlogis(0.65), late_acute=qlogis(0.55), chronic=qlogis(0.45))
# hip does worse
b_joint = c(knee = 0, hip = -0.5)
# surgery

# could do this with ragged, but this is quicker.    
v_early_surg = c(early_na=0)
v_late_surg = c(late_na=0, late_dair = 0, late_one=0.1, late_two=-0.4)
v_chronic_surg = c(chronic_na=0, chronic_one=0, chronic_two=0.2)

v_early_abdu = c(early_na=0)
v_late_abdu = c(late_na=0, late_one_short=0, late_one_long=0.3, late_two_short=0, late_two_long=0.1)
v_chronic_abdu = c(chronic_na=0, chronic_one_short=0, chronic_one_long=-0.2, chronic_two_short=0,chronic_two_long=0.5)

v_early_abty = c(early_na=0, soc = 0, rif = 0.3)
v_late_abty = c(late_na=0, soc = 0, rif = 0.3)
v_chronic_abty = c(chronic_na=0, soc = 0, rif = 0.3)

get_data_demo_2 <- function(
    N = 10000,
    # baseline log-odds of trt success by silo
    b_silo = c(early_acute=qlogis(0.65), late_acute=qlogis(0.55), chronic=qlogis(0.45)),
    # hip does worse
    b_joint = c(knee = 0, hip = -0.5),
    # surgery

    # could do this with ragged, but this is quicker.    
    v_early_surg = c(early_na=0),
    v_late_surg = c(late_na=0, late_dair = 0, late_one=0.1, late_two=-0.4),
    v_chronic_surg = c(chronic_na=0, chronic_one=0, chronic_two=0.2),
    
    v_early_abdu = c(early_na=0),
    v_late_abdu = c(late_na=0, late_one_short=0, late_one_long=0.3, late_two_short=0, late_two_long=0.1),
    v_chronic_abdu = c(chronic_na=0, chronic_one_short=0, chronic_one_long=-0.2, chronic_two_short=0,chronic_two_long=0.5),
    
    v_early_abty = c(early_na=0, soc = 0, rif = 0.3),
    v_late_abty = c(late_na=0, soc = 0, rif = 0.3),
    v_chronic_abty = c(chronic_na=0, soc = 0, rif = 0.3)
    
  ){
  
  # Design options
  fct_early_surg <- c("early_na")
  fct_early_abdu <- c("early_na")
  fct_early_abty <- c("early_na", "soc", "rif")
  
  fct_late_surg <- c("late_na", "late_dair", "late_one", "late_two")
  fct_late_abdu <- c("late_na", "late_one_short", "late_one_long", "late_two_short", "late_two_long")
  fct_late_abty <- c("late_na", "soc", "rif")
  
  fct_chronic_surg <- c("chronic_na", "chronic_one", "chronic_two")
  fct_chronic_abdu <- c("chronic_na", "chronic_one_short", "chronic_one_long", "chronic_two_short", "chronic_two_long")
  fct_chronic_abty <- c("chronic_na", "soc", "rif")

  X_early <- CJ(surg = fct_early_surg, abdu = fct_early_abdu, abty = fct_early_abty)
  X_late <- CJ(surg = fct_late_surg, abdu = fct_late_abdu, abty = fct_late_abty)
  X_chronic <- CJ(surg = fct_chronic_surg, abdu = fct_chronic_abdu, abty = fct_chronic_abty)
  
  # Data
  d <- data.table()
  d[, pt := 1:N]
  
  # Stratification
  silo <- factor(sample(g_silo, size = N, replace = T, prob = g_pr_silo), levels = g_silo)
  d[, silo := silo]
  joint <- factor(sample(g_joint, size = N, replace = T, prob = g_pr_joint), levels = g_joint)
  d[, joint := joint]
  
  # Early
  d[silo == "early_acute", surg := "early_na"]
  d[silo == "early_acute", abdu := "early_na"]
  d[silo == "early_acute", abty := sample(fct_early_abty, .N, replace = T, prob = c(0.5, 0.25, 0.25))]
  # Late
  d[silo == "late_acute", surg := sample(fct_late_surg, .N, replace = T, prob = c(0.05, 0.95/2, 0.95/4, 0.95/4))]
  d[silo == "late_acute" & surg %in% c("late_na", "late_dair"), abdu := "late_na"]
  d[silo == "late_acute" & surg == "late_one", abdu := sample(fct_late_abdu, .N, replace = T, prob = c(0.05, 0.95/2, 0.95/2, 0, 0))]
  d[silo == "late_acute" & surg == "late_two", abdu := sample(fct_late_abdu, .N, replace = T, prob = c(0.05, 0, 0, 0.95/2, 0.95/2))]
  d[silo == "late_acute", abty := sample(fct_late_abty, .N, replace = T, prob = c(0.5, 0.25, 0.25))]
  # Chronic
  d[silo == "chronic", surg := sample(fct_chronic_surg, .N, replace = T, prob = c(0.05, 0.95/2, 0.95/2))]
  d[silo == "chronic" & surg %in% c("chronic_na"), abdu := "chronic_na"]
  d[silo == "chronic" & surg == "chronic_one", abdu := sample(fct_chronic_abdu, .N, replace = T, prob = c(0.05, 0.95/2, 0.95/2, 0, 0))]
  d[silo == "chronic" & surg == "chronic_two", abdu := sample(fct_chronic_abdu, .N, replace = T, prob = c(0.05, 0, 0, 0.95/2, 0.95/2))]
  d[silo == "chronic", abty := sample(fct_chronic_abty, .N, replace = T, prob = c(0.5, 0.25, 0.25))]
  
  d[silo == "early_acute", surg := factor(surg, levels = "early_na")]
  d[silo == "early_acute", abdu := factor(abdu, levels = "early_na")]
  d[silo == "early_acute", abty := factor(abty, levels = fct_early_abty)]
  
  d[silo == "late_acute", surg := factor(surg, levels = fct_late_surg)]
  d[silo == "late_acute", abdu := factor(abdu, levels = fct_late_abdu)]
  d[silo == "late_acute", abty := factor(abty, levels = fct_late_abty)]
  
  d[silo == "chronic", surg := factor(surg, levels = fct_chronic_surg)]
  d[silo == "chronic", abdu := factor(abdu, levels = fct_chronic_abdu)]
  d[silo == "chronic", abty := factor(abty, levels = fct_chronic_abty)]

  # Outcome models.
  d[, eta := b_silo[silo] + b_joint[joint]]
  
  d[silo == "early_acute", eta := eta + v_early_surg[surg]]
  d[silo == "late_acute", eta := eta + v_late_surg[surg]]
  d[silo == "chronic", eta := eta + v_chronic_surg[surg]]
  
  d[silo == "early_acute", eta := eta + v_early_abdu[abdu]]
  d[silo == "late_acute", eta := eta + v_late_abdu[abdu]]
  d[silo == "chronic", eta := eta + v_chronic_abdu[abdu]]
  
  d[silo == "early_acute", eta := eta + v_early_abty[abty]]
  d[silo == "late_acute", eta := eta + v_late_abty[abty]]
  d[silo == "chronic", eta := eta + v_chronic_abty[abty]]
  
  d[, pr_y := plogis(eta)]
  
  d[, y := rbinom(.N, 1, prob = pr_y)]
  
  d
}


get_data <- function(
      ){
  
  log_info(paste0(match.call()[[1]]))
  

  
  d_all <- rbindlist(parallel::mclapply(1:g_trial_cfg$N_sim, FUN = function(i){
    
    d <- data.table()
    d[, pt := 1:g_trial_cfg$N]
    
    # Stratification
    silo <- factor(sample(
      g_trial_cfg$silo, size = g_trial_cfg$N, 
      replace = T, prob = g_trial_cfg$pr_silo), 
      levels = g_trial_cfg$silo)
    d[, silo := silo]
    d[, g_pos := rbinom(.N, 1, g_trial_cfg$pr_gpos)]
    
    # Non-randomised/constant allocation:
    d[silo == "early_acute", surg := "dair"]
    d[silo == "early_acute", abdu := "long"]
    
    # Randomisation
    d[silo == "late_acute", surg := sample(g_fct_la_surg, .N, replace = T, prob = g_trial_cfg$pr_la_surg)]
    d[silo == "late_acute" & surg != "dair", abdu := sample(g_fct_abdu, .N, replace = T)]
    
    # Non-randomised/constant allocation:
    d[silo == "late_acute" & surg == "dair", abdu := "long"]
    
    # Balanced rand for chronic (surg and abdu domain)
    d[silo == "chronic", surg := sample(g_fct_ch_surg, .N, replace = T)]
    d[silo == "chronic" , abdu := sample(g_fct_abdu, .N, replace = T)]
    
    # Balanced rand for abty across all silo
    d[, abty := sample(g_fct_abty, .N, replace = T)]
 
    # Outcome
    d <- merge(d, g_trial_cfg$X, by = c("silo", "surg", "abdu", "abty"), all.x = T)
    d[, pr_y := g_trial_cfg$pr_y[id_x]]
    d[, y := rbinom(.N, 1, prob = pr_y)]
    
    # Enrolment times
    t_enrol <- c(0, nhpp.event.times(g_trial_cfg$nhpoislambda, 
                                g_trial_cfg$N-1, 
                                g_trial_cfg$nhomogintens))
    
    # Precision unimportant
    d[, t_enrol := as.integer(t_enrol)]
    d[, id_x := NULL]
    
    d
    
  }, mc.cores = g_cfg$data_cores), idcol = "sim")

  setcolorder(d_all, c("sim", "pt", "t_enrol", "g_pos"))
  d_all
}


get_rand_status <- function(){
  
  X_surg <- unique(g_trial_cfg$X[, .(silo, surg)])
  X_surg[silo == "early_acute", note := "Not randomised"]
  X_surg[silo == "late_acute" & surg %in% c("one", "two"), note := "Randomised as revision"]
  X_surg[is.na(note), note := "Randomised"]
  
  X_abdu <- unique(g_trial_cfg$X[, .(silo, surg, abdu)])
  X_abdu[silo == "early_acute", note := "Not randomised"]
  X_abdu[silo == "late_acute" & surg == "dair", note := "Not randomised"]
  X_abdu[is.na(note), note := "Randomised"]
  
  X_abty <- unique(g_trial_cfg$X[, .(silo, abty)])
  X_abty[, note := "Randomised (gram pos)"]
  
  list(X_surg = X_surg, X_abdu = X_abdu, X_abty = X_abty)
}




create_db <- function(){
  log_info(paste0(match.call()[[1]]))
  
  if(file.exists(g_cfg$db_name)){
    log_info("Removing existing db")
    ok <- file.remove(g_cfg$db_name)
    stopifnot("Existing file NOT removed" = ok == T)
  }
  
  db <- dbConnect(RSQLite::SQLite(), g_cfg$db_name)
  
  # Main data table
  # db_cols <- unname(c(
  #   "'id'  INTEGER NOT NULL,",
  #   "'sim'  INTEGER NOT NULL,",
  #   "'pt'  INTEGER NOT NULL,",
  #   "'t_enrol'  INTEGER NOT NULL,",
  #   sapply(g_trial_cfg$y_cols, function(z){
  #     paste0("'",z,"'", " INTEGER NOT NULL,")
  #   }),
  #   "PRIMARY KEY('id')"
  # ))
  # cmd <- paste0("CREATE TABLE IF NOT EXISTS 'trial_data' (" , 
  #        paste0(db_cols, collapse = " "), ");")
  
  cmd <- "CREATE TABLE IF NOT EXISTS 'trial_data' (
    'id' INTEGER NOT NULL,
    'sim' INTEGER NOT NULL,
    'pt' INTEGER NOT NULL,
    't_enrol' INTEGER NOT NULL,
    'g_pos' INTEGER NOT NULL,
    'silo' TEXT NOT NULL,
    'surg' TEXT NOT NULL,
    'abdu' TEXT NOT NULL,
    'abty' TEXT NOT NULL,
    'pr_y' REAL NOT NULL,
    'y' INTEGER NOT NULL,
    PRIMARY KEY('id') );
  "
  dbExecute(db, cmd)
  
  
  
  dbExecute(db, "CREATE TABLE IF NOT EXISTS 'X_surg' (
    'id' INTEGER NOT NULL,
    'silo' TEXT NOT NULL,
    'surg' TEXT NOT NULL,
    'note' TEXT NOT NULL,
    PRIMARY KEY('id') );
  ")
  dbExecute(db, "CREATE TABLE IF NOT EXISTS 'X_abdu' (
    'id' INTEGER NOT NULL,
    'silo' TEXT NOT NULL,
    'surg' TEXT NOT NULL,
    'abdu' TEXT NOT NULL,
    'note' TEXT NOT NULL,
    PRIMARY KEY('id') );
  ")
  dbExecute(db, "CREATE TABLE IF NOT EXISTS 'X_abty' (
    'id' INTEGER NOT NULL,
    'silo' TEXT NOT NULL,
    'abty' TEXT NOT NULL,
    'note' TEXT NOT NULL,
    PRIMARY KEY('id') );
  ")
  
  # Views
  
  dbExecute(db, "CREATE VIEW v_surg AS
    SELECT td.sim, td.silo, td.surg, sum(td.y) as y, count(*) as N, Xs.note
    FROM trial_data as td
  	INNER JOIN X_surg as Xs on Xs.silo = td.silo and Xs.surg = td.surg
	  GROUP BY td.sim, td.silo, td.surg;")
  
  dbExecute(db, "CREATE VIEW v_abdu AS
    SELECT td.sim, td.silo, td.surg, td.abdu, sum(td.y) as y, count(*) as N, Xa.note
    FROM trial_data as td
  	INNER JOIN X_abdu as Xa on Xa.silo = td.silo and Xa.surg = td.surg and Xa.abdu = td.abdu
	  GROUP BY td.sim, td.silo, td.surg, td.abdu;")
  
  dbExecute(db, "CREATE VIEW v_abty AS
    SELECT td.sim, td.silo, td.abty, sum(td.y) as y, count(*) as N, Xa.note
    FROM trial_data as td
  	INNER JOIN X_abty as Xa on Xa.silo = td.silo and Xa.abty = td.abty and td.g_pos == 1
    GROUP BY td.sim, td.silo, td.abty;")
  
  
  
  # END
  
  dbDisconnect(db)
  
}

populate_db <- function(){
  log_info(paste0(match.call()[[1]]))
  d_all <- get_data()
  db <- dbConnect(RSQLite::SQLite(), g_cfg$db_name)
  dbAppendTable(db, "trial_data", d_all)

  l_rand <- get_rand_status()
  dbAppendTable(db, "X_surg", l_rand$X_surg)
  dbAppendTable(db, "X_abdu", l_rand$X_abdu)
  dbAppendTable(db, "X_abty", l_rand$X_abty)
  
  # END
  dbDisconnect(db)
}

run_data <- function(){
  log_info(paste0(match.call()[[1]]))
  create_db()
  populate_db()
}

run_none_data <- function(){
  log_info("run_none_data: Nothing doing here bud.")
}

main_data <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_data()