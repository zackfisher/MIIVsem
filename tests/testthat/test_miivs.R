library("MIIVsem")
context("miivs")

##-----------------------------------------------------------##
##-----------------------------------------------------------##
##--------------------Political Demo Test--------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##
bollen1989a_model <- '
Xi1  =~ x1 + x2 + x3
Eta1 =~ y1 + y2 + y3 + y4
Eta2 =~ y5 + y6 + y7 + y8
Eta1 ~ Xi1
Eta2 ~ Xi1 + Eta1
y1 ~~ y5
y2 ~~ y4 + y6
y3 ~~ y7
y4 ~~ y8
y6 ~~ y8
'

bollen1989a_miivs <-
  list(
    structure(
      list(
        DVobs = "y1",
        IVobs = "x1",
        IV = c("x2", "x3")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "y5",
        IVobs = c("y1", "x1"),
        IV = c("y2", "y3", "y4", "x2", "x3")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "y2",
        IVobs = "y1",
        IV = c("y3", "y7", "y8", "x2", "x3", "x1")
      ),
      .Names = c("DVobs",  "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "y3",
        IVobs = "y1",
        IV = c("y2",  "y4", "y6", "y8", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs",  "IV")
    ),
    structure(
      list(
        DVobs = "y4",
        IVobs = "y1",
        IV = c("y3",  "y6", "y7", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "y6",
        IVobs = "y5",
        IV = c("y3", "y4", "y7", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs",  "IV")
    ),
    structure(
      list(
        DVobs = "y7",
        IVobs = "y5",
        IV = c("y2", "y4", "y6", "y8", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs",  "IV")
    ),
    structure(
      list(
        DVobs = "y8",
        IVobs = "y5",
        IV = c("y2",  "y3", "y7", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs",  "IV")
    ),
    structure(
      list(
        DVobs = "x2",
        IVobs = "x1",
        IV = c("y1", "y5", "y2", "y3", "y4", "y6", "y7", "y8", "x3")
      ),
      .Names = c("DVobs",  "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x3",
        IVobs = "x1",
        IV = c("y1", "y5", "y2", "y3", "y4", "y6", "y7", "y8", "x2")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    )
  )

bollen1989a_miivs_test <- lapply(miivs(bollen1989a_model)$eqns, "[", c("DVobs", "IVobs", "IV"))

test_that("pol demo miivs correct", {
  expect_equal(bollen1989a_miivs, bollen1989a_miivs_test)
})


##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Colom Higher-Order Test------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

colom_model <- '
# General Intelligence
g  =~ 1*WM + PS + CI + SA
# Working Memory
WM =~ 1*x1 + x2 + x3
# Processing Speed
PS =~ 1*x4 + x5 + x6 
# Crystalized Intelligence
CI =~ 1*x7 + x8 + x9
# Spatial Abilities
SA =~ 1*x10 + x11 + x12
'
colom_miivs_test <- lapply(miivs(colom_model)$eqns, "[", c("DVobs", "IVobs", "IV"))

colom_miivs <-
  list(
    structure(
      list(
        DVobs = "x4",
        IVobs = "x1",
        IV = c("x7", "x10", "x8", "x9", "x11", "x12") 
      ),
      .Names = c("DVobs", "IVobs",
                 "IV")
    ),
    structure(
      list(
        DVobs = "x7",
        IVobs = "x1",
        IV = c("x4", "x10", "x5", "x6", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x10",
        IVobs = "x1",
        IV = c("x4", "x7", "x5", "x6", "x8", "x9")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x2",
        IVobs = "x1",
        IV = c("x4", "x7", "x10", "x3", "x5", "x6", "x8", "x9", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x3",
        IVobs = "x1",
        IV = c("x4", "x7", "x10", "x2", "x5", "x6", "x8", "x9", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x5",
        IVobs = "x4",
        IV = c("x1", "x7", "x10", "x2", "x3", "x6", "x8", "x9", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x6",
        IVobs = "x4",
        IV = c("x1", "x7", "x10", "x2", "x3", "x5", "x8", "x9", "x11", "x12")
      ),
      .Names = c("DVobs",
                 "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x8",
        IVobs = "x7",
        IV = c("x1", "x4", "x10", "x2", "x3", "x5", "x6", "x9", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x9",
        IVobs = "x7",
        IV = c("x1", "x4", "x10", "x2", "x3", "x5", "x6", "x8", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x11",
        IVobs = "x10",
        IV = c("x1", "x4", "x7", "x2", "x3", "x5", "x6", "x8", "x9", "x12")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    ),
    structure(
      list(
        DVobs = "x12",
        IVobs = "x10",
        IV = c("x1", "x4", "x7", "x2", "x3", "x5", "x6", "x8", "x9", "x11")
      ),
      .Names = c("DVobs", "IVobs", "IV")
    )
  )

test_that("simul_reg_1 is correct", {
  expect_equal(colom_miivs , colom_miivs_test)
})

##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Simultaneous Eq Reg 1--------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

simul_reg_1 <- '
  educ ~ paeducM + maeducM + paNOeduc + maNOeduc + 
         napaeduc +  namaeduc + papr80M + mapr80M +
         napapr80 + namapr80 + age + lnsibs + female + 
         black + asian + hispanic + othrace
  prestg80 ~ educ
'

simul_reg_1_miivs_test <- lapply(miivs(simul_reg_1)$eqns, "[", c("DVobs", "IVobs", "IV"))

simul_reg_1_miivs <- 
  list(structure(
    list(
      DVobs = "educ",
      IVobs = c(
        "paeducM",
        "maeducM",
        "paNOeduc",
        "maNOeduc",
        "napaeduc",
        "namaeduc",
        "papr80M",
        "mapr80M",
        "napapr80",
        "namapr80",
        "age",
        "lnsibs",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace"
      ),
      IV = c(
        "paeducM",
        "maeducM",
        "paNOeduc",
        "maNOeduc",
        "napaeduc",
        "namaeduc",
        "papr80M",
        "mapr80M",
        "napapr80",
        "namapr80",
        "age",
        "lnsibs",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace"
      )
    ),
    .Names = c("DVobs", "IVobs", "IV")
  ), structure(
    list(
      DVobs = "prestg80",
      IVobs = "educ",
      IV = c(
        "paeducM",
        "maeducM",
        "paNOeduc",
        "maNOeduc",
        "napaeduc",
        "namaeduc",
        "papr80M",
        "mapr80M",
        "napapr80",
        "namapr80",
        "age",
        "lnsibs",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace",
        "educ"
      )
    ),
    .Names = c("DVobs",
               "IVobs", "IV")
  ))

test_that("simul_reg_1 is correct", {
  expect_equal(simul_reg_1_miivs, simul_reg_1_miivs_test)
})

##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Simultaneous Eq Reg 2--------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

simul_reg_2 <- '
  educ ~ paeducM + maeducM + paNOeduc + maNOeduc + napaeduc + 
         namaeduc + papr80M + mapr80M + napapr80 + namapr80 + 
         age + lnsibs + female + black + asian + hispanic + othrace
  prestg80 ~ educ + napapr80 + papr80M + namapr80 + mapr80M + 
             age + female + black +  asian  + hispanic + othrace
'

simul_reg_2_miivs_test <- lapply(miivs(simul_reg_2)$eqns, "[", c("DVobs", "IVobs", "IV"))

simul_reg_2_miivs <- 
list(structure(
  list(
    DVobs = "educ",
    IVobs = c(
      "paeducM",
      "maeducM",
      "paNOeduc",
      "maNOeduc",
      "napaeduc",
      "namaeduc",
      "papr80M",
      "mapr80M",
      "napapr80",
      "namapr80",
      "age",
      "lnsibs",
      "female",
      "black",
      "asian",
      "hispanic",
      "othrace"
    ),
    IV = c(
      "paeducM",
      "maeducM",
      "paNOeduc",
      "maNOeduc",
      "napaeduc",
      "namaeduc",
      "papr80M",
      "mapr80M",
      "napapr80",
      "namapr80",
      "age",
      "lnsibs",
      "female",
      "black",
      "asian",
      "hispanic",
      "othrace"
    )
  ),
  .Names = c("DVobs", "IVobs", "IV")
), structure(
  list(
    DVobs = "prestg80",
    IVobs = c(
      "educ",
      "napapr80",
      "papr80M",
      "namapr80",
      "mapr80M",
      "age",
      "female",
      "black",
      "asian",
      "hispanic",
      "othrace"
    ),
    IV = c(
      "paeducM",
      "maeducM",
      "paNOeduc",
      "maNOeduc",
      "napaeduc",
      "namaeduc",
      "papr80M",
      "mapr80M",
      "napapr80",
      "namapr80",
      "age",
      "lnsibs",
      "female",
      "black",
      "asian",
      "hispanic",
      "othrace",
      "educ"
    )
  ),
  .Names = c("DVobs", "IVobs", "IV")
))

test_that("simul_reg_2 is correct", {
  expect_equal(simul_reg_2_miivs, simul_reg_2_miivs_test)
})

##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Simultaneous Eq Reg 3--------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

simul_reg_3 <- '
educ ~ paeducM + maeducM + paNOeduc + maNOeduc + 
napaeduc +  namaeduc + papr80M + mapr80M +
napapr80 + namapr80 + age + lnsibs + female + 
black + asian + hispanic + othrace
prestg80 ~ educ
prestg80 ~~ educ
'

simul_reg_3_miivs_test <- lapply(miivs(simul_reg_3)$eqns, "[", c("DVobs", "IVobs", "IV"))

simul_reg_3_miivs <-
  list(structure(
    list(
      DVobs = "educ",
      IVobs = c(
        "paeducM",
        "maeducM",
        "paNOeduc",
        "maNOeduc",
        "napaeduc",
        "namaeduc",
        "papr80M",
        "mapr80M",
        "napapr80",
        "namapr80",
        "age",
        "lnsibs",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace"
      ),
      IV = c(
        "paeducM",
        "maeducM",
        "paNOeduc",
        "maNOeduc",
        "napaeduc",
        "namaeduc",
        "papr80M",
        "mapr80M",
        "napapr80",
        "namapr80",
        "age",
        "lnsibs",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace"
      )
    ),
    .Names = c("DVobs", "IVobs", "IV")
  ), structure(
    list(
      DVobs = "prestg80",
      IVobs = "educ",
      IV = c(
        "paeducM",
        "maeducM",
        "paNOeduc",
        "maNOeduc",
        "napaeduc",
        "namaeduc",
        "papr80M",
        "mapr80M",
        "napapr80",
        "namapr80",
        "age",
        "lnsibs",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace"
      )
    ),
    .Names = c("DVobs",
               "IVobs", "IV")
  ))
test_that("simul_reg_3 is correct", {
  expect_equal(simul_reg_3_miivs, simul_reg_3_miivs_test)
})


##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Simultaneous Eq Reg 4--------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

simul_reg_4 <- '
educ ~ paeducM + maeducM + paNOeduc + maNOeduc + napaeduc + 
namaeduc + papr80M + mapr80M + napapr80 + namapr80 + 
age + lnsibs + female + black + asian + hispanic + othrace
prestg80 ~ educ + napapr80 + papr80M + namapr80 + mapr80M + 
age + female + black +  asian  + hispanic + othrace
prestg80 ~~ educ
'

simul_reg_4_miivs_test <- lapply(miivs(simul_reg_4)$eqns, "[", c("DVobs", "IVobs", "IV"))

simul_reg_4_miivs <- 
  list(structure(
    list(
      DVobs = "educ",
      IVobs = c(
        "paeducM",
        "maeducM",
        "paNOeduc",
        "maNOeduc",
        "napaeduc",
        "namaeduc",
        "papr80M",
        "mapr80M",
        "napapr80",
        "namapr80",
        "age",
        "lnsibs",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace"
      ),
      IV = c(
        "paeducM",
        "maeducM",
        "paNOeduc",
        "maNOeduc",
        "napaeduc",
        "namaeduc",
        "papr80M",
        "mapr80M",
        "napapr80",
        "namapr80",
        "age",
        "lnsibs",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace"
      )
    ),
    .Names = c("DVobs", "IVobs", "IV")
  ), structure(
    list(
      DVobs = "prestg80",
      IVobs = c(
        "educ",
        "napapr80",
        "papr80M",
        "namapr80",
        "mapr80M",
        "age",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace"
      ),
      IV = c(
        "paeducM",
        "maeducM",
        "paNOeduc",
        "maNOeduc",
        "napaeduc",
        "namaeduc",
        "papr80M",
        "mapr80M",
        "napapr80",
        "namapr80",
        "age",
        "lnsibs",
        "female",
        "black",
        "asian",
        "hispanic",
        "othrace"
      )
    ),
    .Names = c("DVobs", "IVobs",
               "IV")
  ))
test_that("simul_reg_4 is correct", {
  expect_equal(simul_reg_4_miivs, simul_reg_4_miivs_test)
})

