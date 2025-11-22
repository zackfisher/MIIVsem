library("MIIVsem")

data("reisenzein1986")

model <- '
    #### latent variable definitions (x must come first then y)
Control =~ Z2 + Z3 + Z4
Sympathy =~ Z5 + Z6 + Z7
Anger =~ Z8 + Z9 + Z10
Help =~ Z11 + Z12 + Z13

#### regressions
Sympathy ~ Control
Anger ~ Control
Help ~ Sympathy + Anger
'

# part1 ----
fit <- miive( model, data = reisenzein1986)
# print(fit)
# summary(fit)
# print(fit,meanvar.df=F)
# print(fit,meanvar.df=T)
# summary(fit,meanvar.df=F)
# summary(fit,meanvar.df=T)
# miive( model, data = reisenzein1986,overid = "meanvar")
# miive( model, data = reisenzein1986,overid = "mean")
# miive( model, data = reisenzein1986,overid = "adjusted")
# miive( model, data = reisenzein1986,overid = "classic")

#-------------------------------------------------------#  
context("ex05: overidentification test statistics correct (continuous variables)")
#-------------------------------------------------------# 

sargan <- unlist(lapply(fit$eqn, "[", c("test.stat")))
# names(sargan) <- unlist(lapply(fit$eqn, "[", c("DVlat")))
# sargan   <- sargan[sort(names(sargan))]

expect_equal_to_reference(
  sargan, 
  "rds/ex05_polr_sargan_con.rds"
)

# saveRDS(
#   sargan,
#   file = "rds/ex05_polr_sargan_con.rds")

# part2 ----
####  Make some 10-category variables as 5-category ####
####  Only Z13 is 5-category. Others are 10-category.
sapply(2 : 13, function(z) table(reisenzein1986[,z]))
####  Z5-Z7, Z11-Z13 remain categorical. Others are treated as continuous.
########  Z5
j <- 5
table(reisenzein1986[,j])
reisenzein1986[ which(reisenzein1986[,j] %in% c(1, 2)), j ] <- 1
reisenzein1986[ which(reisenzein1986[,j] %in% c(3, 4, 5)), j ] <- 2
reisenzein1986[ which(reisenzein1986[,j] %in% c(6, 7)), j ] <- 3
reisenzein1986[ which(reisenzein1986[,j] %in% c(8)), j ] <- 4
reisenzein1986[ which(reisenzein1986[,j] %in% c(9)), j ] <- 5
table(reisenzein1986[,j])
########  Z6
j <- 6
table(reisenzein1986[,j])
reisenzein1986[ which(reisenzein1986[,j] %in% c(1, 2)), j ] <- 1
reisenzein1986[ which(reisenzein1986[,j] %in% c(3, 4, 5)), j ] <- 2
reisenzein1986[ which(reisenzein1986[,j] %in% c(6, 7)), j ] <- 3
reisenzein1986[ which(reisenzein1986[,j] %in% c(8)), j ] <- 4
reisenzein1986[ which(reisenzein1986[,j] %in% c(9)), j ] <- 5
table(reisenzein1986[,j])
########  Z7
j <- 7
table(reisenzein1986[,j])
reisenzein1986[ which(reisenzein1986[,j] %in% c(1, 2)), j ] <- 1
reisenzein1986[ which(reisenzein1986[,j] %in% c(3, 4, 5)), j ] <- 2
reisenzein1986[ which(reisenzein1986[,j] %in% c(6, 7)), j ] <- 3
reisenzein1986[ which(reisenzein1986[,j] %in% c(8)), j ] <- 4
reisenzein1986[ which(reisenzein1986[,j] %in% c(9)), j ] <- 5
table(reisenzein1986[,j])
########  Z11
j <- 11
table(reisenzein1986[,j])
reisenzein1986[ which(reisenzein1986[,j] %in% c(1, 2)), j ] <- 1
reisenzein1986[ which(reisenzein1986[,j] %in% c(3, 4, 5)), j ] <- 2
reisenzein1986[ which(reisenzein1986[,j] %in% c(6, 7)), j ] <- 3
reisenzein1986[ which(reisenzein1986[,j] %in% c(8)), j ] <- 4
reisenzein1986[ which(reisenzein1986[,j] %in% c(9)), j ] <- 5
table(reisenzein1986[,j])
########  Z12
j <- 12
table(reisenzein1986[,j])
reisenzein1986[ which(reisenzein1986[,j] %in% c(1, 2)), j ] <- 1
reisenzein1986[ which(reisenzein1986[,j] %in% c(3, 4, 5)), j ] <- 2
reisenzein1986[ which(reisenzein1986[,j] %in% c(6, 7)), j ] <- 3
reisenzein1986[ which(reisenzein1986[,j] %in% c(8)), j ] <- 4
reisenzein1986[ which(reisenzein1986[,j] %in% c(9)), j ] <- 5
table(reisenzein1986[,12])
reisenzein1986 <- data.frame(reisenzein1986[,-1])

fit <- miive( model, data = reisenzein1986, 
                   ordered = paste("Z", c(5 : 7, 11 : 13), sep = ""))

#-------------------------------------------------------#  
context("ex05: overidentification test statistics correct (continuous and categorical variables)")
#-------------------------------------------------------# 

sargan <- unlist(lapply(fit$eqn, "[", c("test.stat")))
# names(sargan) <- unlist(lapply(fit$eqn, "[", c("DVlat")))
# sargan   <- sargan[sort(names(sargan))]

expect_equal_to_reference(
  sargan, 
  "rds/ex05_polr_sargan_cate.rds"
)

# saveRDS(
#   sargan,
#   file = "rds/ex05_polr_sargan_cate.rds")