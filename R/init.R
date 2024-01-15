# load libs, load config, initialise log

suppressPackageStartupMessages(library("cmdstanr"))
suppressPackageStartupMessages(library("config"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("DBI"))
suppressPackageStartupMessages(library("DT"))
suppressPackageStartupMessages(library("RSQLite"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(suppressWarnings(library("gt")))
suppressPackageStartupMessages(library("git2r"))
suppressPackageStartupMessages(library("gt"))
suppressPackageStartupMessages(library("jsonlite"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("logger"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("poisson"))
suppressPackageStartupMessages(library("pracma"))
suppressPackageStartupMessages(library("roadmap.data"))


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


