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