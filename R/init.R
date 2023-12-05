# load libs, load config, initialise log

suppressPackageStartupMessages(library("logger"))
suppressPackageStartupMessages(library("DBI"))
suppressPackageStartupMessages(library("RSQLite"))
suppressPackageStartupMessages(library("config"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(
  suppressWarnings(library("gt")))
suppressPackageStartupMessages(library("git2r"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("jsonlite"))
suppressPackageStartupMessages(library("pracma"))
suppressPackageStartupMessages(library("cmdstanr"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("poisson"))

ggplot2::theme_set(theme_bw() + theme(
  legend.position = "bottom"
))

is_html <- knitr::is_html_output()
# is_pdf <- knitr::is_latex_output()
# is_word <- !is_html & !is_pdf


# Config - store of local OS file system
f_cfg <- file.path("./etc", "config.yml")
g_cfg <- config::get(file = f_cfg)
stopifnot("Config is null" = !is.null(g_cfg))


# message(Sys.time(), " Config read from ", f_cfg)

# Logs
f_log <- file.path("./logs", "log.txt")
log_appender(appender_file(f_log))
# message(Sys.time(), " Log file initialised ", f_log)
log_info("*** START UP ***")


