# load libs, load config, initialise log

source("R/libs.R")

is_html <- knitr::is_html_output()
# is_pdf <- knitr::is_latex_output()
# is_word <- !is_html & !is_pdf


# Config - store of local OS file system
f_cfg <- file.path("./etc", "cfg.yml")
g_cfg <- config::get(file = f_cfg)
stopifnot("Config is null" = !is.null(g_cfg))


# message(Sys.time(), " Config read from ", f_cfg)

# Logs
f_log <- file.path("./logs", "log.txt")
log_appender(appender_file(f_log))
# message(Sys.time(), " Log file initialised ", f_log)
log_info("*** START UP ***")

# pars/effects of interest
g_pars <- c(
  "alpha", "gamma_c", 
  "b_a_l", 
  "b_b1_l", "b_b2_l","b_a_c","b_b1_c", "b_b2_c",
  "b_c"
)
g_effs <- c(
  "b_a_l_2", "b_a_c_2",  # domain a, late (revision) and chronic (two-stage)
  "b_b1_l_2", # domain b, (late/revision one stage pts) wk12p1 (ref is wk6p1)
  "b_b2_l_2", # wk12p2 (ref is day7p2)
  "b_b1_c_2", # wk12p1 (ref is wk6p1)
  "b_b2_c_2", # wk12p2 (ref is day7p2)
  "b_c_2" # rif (ref is no-rif)
)
