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


