setwd('~/Desktop/stratSSL_Revision_Results/old_result_S2')
n.t <- 200
my_beta_true_200 <- as.matrix(read.table( paste0(n.t, "-beta_true.csv")))
my_beta_sl_200 <- as.matrix(read.table( paste0(n.t, "-beta_sl.csv")))
my_beta_ssl_200 <- as.matrix(read.table( paste0(n.t, "-beta_ssl.csv")))
my_beta_ssl2_200 <- as.matrix(read.table( paste0(n.t, "-beta_ssl2.csv")))

my_beta_true_200_m <- colMeans(my_beta_true_200)
my_beta_sl_200_m <- colMeans(my_beta_sl_200)
my_beta_ssl_200_m <- colMeans(my_beta_ssl_200)
my_beta_ssl2_200_m <- colMeans(my_beta_ssl2_200)

my_beta_true_200_s <- colSds(my_beta_true_200)
my_beta_sl_200_s <- colSds(my_beta_sl_200)
my_beta_ssl_200_s <- colSds(my_beta_ssl_200)
my_beta_ssl2_200_s <- colSds(my_beta_ssl2_200)

bias_sl_S2 <- my_beta_sl_200_m - my_beta_true_200_m
bias_ssl_S2 <- my_beta_ssl_200_m - my_beta_true_200_m
bias_ssl2_S2 <- my_beta_ssl2_200_m - my_beta_true_200_m
bias_200_S2 <- cbind(bias_sl_S2, bias_ssl_S2, bias_ssl2_S2)
bias_200_S2

n.t <- 400
my_beta_true_400 <- as.matrix(read.table( paste0(n.t, "-beta_true.csv")))
my_beta_sl_400 <- as.matrix(read.table( paste0(n.t, "-beta_sl.csv")))
my_beta_ssl_400 <- as.matrix(read.table( paste0(n.t, "-beta_ssl.csv")))
my_beta_ssl2_400 <- as.matrix(read.table( paste0(n.t, "-beta_ssl2.csv")))


my_beta_true_400_m <- colMeans(my_beta_true_400)
my_beta_sl_400_m <- colMeans(my_beta_sl_400)
my_beta_ssl_400_m <- colMeans(my_beta_ssl_400)
my_beta_ssl2_400_m <- colMeans(my_beta_ssl2_400)


my_beta_true_400_s <- colSds(my_beta_true_400)
my_beta_sl_400_s <- colSds(my_beta_sl_400)
my_beta_ssl_400_s <- colSds(my_beta_ssl_400)
my_beta_ssl2_400_s <- colSds(my_beta_ssl2_400)

bias_sl_S2 <- my_beta_sl_400_m - my_beta_true_400_m
bias_ssl_S2 <- my_beta_ssl_400_m - my_beta_true_400_m
bias_ssl2_S2 <- my_beta_ssl2_400_m - my_beta_true_400_m
bias_400_S2 <- cbind(bias_sl_S2, bias_ssl_S2, bias_ssl2_S2)
bias_400_S2

setwd('~/Desktop/stratSSL_Revision_Results/new_result_S4')
n.t <- 200
my_beta_true_200 <- as.matrix(read.table( paste0(n.t, "-beta_true.csv")))
my_beta_sl_200 <- as.matrix(read.table( paste0(n.t, "-beta_sl.csv")))
my_beta_ssl_200 <- as.matrix(read.table( paste0(n.t, "-beta_ssl.csv")))
my_beta_ssl2_200 <- as.matrix(read.table( paste0(n.t, "-beta_ssl2.csv")))

my_beta_true_200_m <- colMeans(my_beta_true_200)
my_beta_sl_200_m <- colMeans(my_beta_sl_200)
my_beta_ssl_200_m <- colMeans(my_beta_ssl_200)
my_beta_ssl2_200_m <- colMeans(my_beta_ssl2_200)

my_beta_true_200_s <- colSds(my_beta_true_200)
my_beta_sl_200_s <- colSds(my_beta_sl_200)
my_beta_ssl_200_s <- colSds(my_beta_ssl_200)
my_beta_ssl2_200_s <- colSds(my_beta_ssl2_200)

bias_sl_S4 <- my_beta_sl_200_m - my_beta_true_200_m
bias_ssl_S4 <- my_beta_ssl_200_m - my_beta_true_200_m
bias_ssl2_S4 <- my_beta_ssl2_200_m - my_beta_true_200_m
bias_200_S4 <- cbind(bias_sl_S4, bias_ssl_S4, bias_ssl2_S4)
bias_200_S4

n.t <- 400
my_beta_true_400 <- as.matrix(read.table( paste0(n.t, "-beta_true.csv")))
my_beta_sl_400 <- as.matrix(read.table( paste0(n.t, "-beta_sl.csv")))
my_beta_ssl_400 <- as.matrix(read.table( paste0(n.t, "-beta_ssl.csv")))
my_beta_ssl2_400 <- as.matrix(read.table( paste0(n.t, "-beta_ssl2.csv")))


my_beta_true_400_m <- colMeans(my_beta_true_400)
my_beta_sl_400_m <- colMeans(my_beta_sl_400)
my_beta_ssl_400_m <- colMeans(my_beta_ssl_400)
my_beta_ssl2_400_m <- colMeans(my_beta_ssl2_400)


my_beta_true_400_s <- colSds(my_beta_true_400)
my_beta_sl_400_s <- colSds(my_beta_sl_400)
my_beta_ssl_400_s <- colSds(my_beta_ssl_400)
my_beta_ssl2_400_s <- colSds(my_beta_ssl2_400)

bias_sl_S4 <- my_beta_sl_400_m - my_beta_true_400_m
bias_ssl_S4 <- my_beta_ssl_400_m - my_beta_true_400_m
bias_ssl2_S4 <- my_beta_ssl2_400_m - my_beta_true_400_m
bias_400_S4 <- cbind(bias_sl_S4, bias_ssl_S4, bias_ssl2_S4)
bias_400_S4

round(bias_400_S2,4)
round(bias_400_S4,4)



# create the table
options("scipen" = 100, "digits" = 5)
bias_table <- cbind(round(abs(bias_200_S2[, 3]  - bias_200_S2[, 2]), 4),
      round(abs(bias_400_S2[, 3]  - bias_400_S2[, 2]), 4),
      round(abs(bias_200_S4[, 3]  - bias_200_S4[, 2]), 4),
      round(abs(bias_400_S4[, 3]  - bias_400_S4[, 2]), 4)
)
colnames(bias_table) <- rep(c('$n = 100$', '$n = 100$'), 2)
rownames(bias_table) <- paste0('$\\theta_', 0:10, '$')

library(kableExtra)
kable(bias_table, "latex", booktabs = T, escape = FALSE) %>%
  add_header_above(c(" ", "S = 2"=2,"S = 4"=2))



