setwd('~/Desktop/new_result_S4')
n.t <- 200
my_beta_true_200 <- as.matrix(read.table( paste0(n.t, "-beta_true.csv")))
my_beta_sl_200 <- as.matrix(read.table( paste0(n.t, "-beta_sl.csv")))
my_beta_ssl_200 <- as.matrix(read.table( paste0(n.t, "-beta_ssl.csv")))
my_beta_ssl2_200 <- as.matrix(read.table( paste0(n.t, "-beta_ssl2.csv")))

my_beta_true_200_m <- colMeans(my_beta_true_200)
my_beta_sl_200_m <- colMeans(my_beta_sl_200)
my_beta_ssl_200_m <- colMeans(my_beta_ssl_200)
my_beta_ssl2_200_m <- colMeans(my_beta_ssl2_200)

my_beta_true_200_m
my_beta_sl_200_m
my_beta_ssl_200_m
my_beta_ssl2_200_m

my_beta_true_200_s <- colSds(my_beta_true_200)
my_beta_sl_200_s <- colSds(my_beta_sl_200)
my_beta_ssl_200_s <- colSds(my_beta_ssl_200)
my_beta_ssl2_200_s <- colSds(my_beta_ssl2_200)

my_beta_true_200_s
my_beta_sl_200_s
my_beta_ssl_200_s
my_beta_ssl2_200_s

bias_sl <- my_beta_sl_200_m - my_beta_true_200_m
bias_ssl <- my_beta_ssl_200_m - my_beta_true_200_m
bias_ssl2 <- my_beta_ssl2_200_m - my_beta_true_200_m
bias_200 <- round(cbind(bias_sl, bias_ssl, bias_ssl2), 2)
bias_200

n.t <- 400
my_beta_true_400 <- as.matrix(read.table( paste0(n.t, "-beta_true.csv")))
my_beta_sl_400 <- as.matrix(read.table( paste0(n.t, "-beta_sl.csv")))
my_beta_ssl_400 <- as.matrix(read.table( paste0(n.t, "-beta_ssl.csv")))
my_beta_ssl2_400 <- as.matrix(read.table( paste0(n.t, "-beta_ssl2.csv")))


my_beta_true_400_m <- colMeans(my_beta_true_400)
my_beta_sl_400_m <- colMeans(my_beta_sl_400)
my_beta_ssl_400_m <- colMeans(my_beta_ssl_400)
my_beta_ssl2_400_m <- colMeans(my_beta_ssl2_400)

my_beta_true_400_m
my_beta_sl_400_m
my_beta_ssl_400_m
my_beta_ssl2_400_m

my_beta_true_400_s <- colSds(my_beta_true_400)
my_beta_sl_400_s <- colSds(my_beta_sl_400)
my_beta_ssl_400_s <- colSds(my_beta_ssl_400)
my_beta_ssl2_400_s <- colSds(my_beta_ssl2_400)

my_beta_true_400_s
my_beta_sl_400_s
my_beta_ssl_400_s
my_beta_ssl2_400_s

bias_sl <- my_beta_sl_400_m - my_beta_true_400_m
bias_ssl <- my_beta_ssl_400_m - my_beta_true_400_m
bias_ssl2 <- my_beta_ssl2_400_m - my_beta_true_400_m
bias_400 <- round(cbind(bias_sl, bias_ssl, bias_ssl2), 2)
bias_400


