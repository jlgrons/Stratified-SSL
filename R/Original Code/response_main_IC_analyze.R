setwd('~/Desktop')
n.t <- 200
my_beta_true_200 <- read.table( paste0(n.t, "-beta_true.csv"))
my_beta_sl_200 <- read.table( paste0(n.t, "-beta_sl.csv"))
my_beta_ssl_200 <- read.table( paste0(n.t, "-beta_ssl.csv"))
my_beta_ssl2_200 <- read.table( paste0(n.t, "-beta_ssl2.csv"))

my_beta_true_200_m <- colMeans(my_beta_true_200)
my_beta_sl_200_m <- colMeans(my_beta_sl_200)
my_beta_ssl_200_m <- colMeans(my_beta_ssl_200)
my_beta_ssl2_200_m <- colMeans(my_beta_ssl2_200)

my_beta_true_200_m
my_beta_sl_200_m
my_beta_ssl_200_m
my_beta_ssl2_200_m

setwd('~/Desktop')
n.t <- 400
my_beta_true_400 <- read.table( paste0(n.t, "-beta_true.csv"))
my_beta_sl_400 <- read.table( paste0(n.t, "-beta_sl.csv"))
my_beta_ssl_400 <- read.table( paste0(n.t, "-beta_ssl.csv"))
my_beta_ssl2_400 <- read.table( paste0(n.t, "-beta_ssl2.csv"))
