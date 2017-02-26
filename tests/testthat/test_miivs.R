library("MIIVsem")
context("miivs")

mySort <- function(a){
  a <- lapply(a, function(b){
    lapply(b, sort) 
  })
  
  a[order(unlist(lapply(a,function(b){b$DVobs})))]
}

expect_equal_sort <- function(a,b){
  expect_equal(mySort(a), mySort(b))
}

##-----------------------------------------------------------##
##--------------------Political Demo Test--------------------##
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
        MIIVs = c("x2", "x3")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "y5",
        IVobs = c("x1", "y1"),
        MIIVs = c("y2", "y3", "y4", "x2", "x3")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "y2",
        IVobs = "y1",
        MIIVs = c("y3", "y7", "y8", "x2", "x3", "x1")
      ),
      .Names = c("DVobs",  "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "y3",
        IVobs = "y1",
        MIIVs = c("y2",  "y4", "y6", "y8", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs",  "MIIVs")
    ),
    structure(
      list(
        DVobs = "y4",
        IVobs = "y1",
        MIIVs = c("y3",  "y6", "y7", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "y6",
        IVobs = "y5",
        MIIVs = c("y3", "y4", "y7", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs",  "MIIVs")
    ),
    structure(
      list(
        DVobs = "y7",
        IVobs = "y5",
        MIIVs = c("y2", "y4", "y6", "y8", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs",  "MIIVs")
    ),
    structure(
      list(
        DVobs = "y8",
        IVobs = "y5",
        MIIVs = c("y2",  "y3", "y7", "x2", "x3", "x1")
      ),
      .Names = c("DVobs", "IVobs",  "MIIVs")
    ),
    structure(
      list(
        DVobs = "x2",
        IVobs = "x1",
        MIIVs = c("y1", "y5", "y2", "y3", "y4", "y6", "y7", "y8", "x3")
      ),
      .Names = c("DVobs",  "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x3",
        IVobs = "x1",
        MIIVs = c("y1", "y5", "y2", "y3", "y4", "y6", "y7", "y8", "x2")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    )
  )

bollen1989a_miivs_test <- lapply(miivs(bollen1989a_model)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))

# Sort both lists

test_that("pol demo miivs correct", {
  expect_equal_sort(bollen1989a_miivs, bollen1989a_miivs_test)
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
colom_miivs_test <- lapply(miivs(colom_model)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))

colom_miivs <-
  list(
    structure(
      list(
        DVobs = "x4",
        IVobs = "x1",
        MIIVs = c("x7", "x10", "x8", "x9", "x11", "x12") 
      ),
      .Names = c("DVobs", "IVobs",
                 "MIIVs")
    ),
    structure(
      list(
        DVobs = "x7",
        IVobs = "x1",
        MIIVs = c("x4", "x10", "x5", "x6", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x10",
        IVobs = "x1",
        MIIVs = c("x4", "x7", "x5", "x6", "x8", "x9")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x2",
        IVobs = "x1",
        MIIVs = c("x4", "x7", "x10", "x3", "x5", "x6", "x8", "x9", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x3",
        IVobs = "x1",
        MIIVs = c("x4", "x7", "x10", "x2", "x5", "x6", "x8", "x9", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x5",
        IVobs = "x4",
        MIIVs = c("x1", "x7", "x10", "x2", "x3", "x6", "x8", "x9", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x6",
        IVobs = "x4",
        MIIVs = c("x1", "x7", "x10", "x2", "x3", "x5", "x8", "x9", "x11", "x12")
      ),
      .Names = c("DVobs",
                 "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x8",
        IVobs = "x7",
        MIIVs = c("x1", "x4", "x10", "x2", "x3", "x5", "x6", "x9", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x9",
        IVobs = "x7",
        MIIVs = c("x1", "x4", "x10", "x2", "x3", "x5", "x6", "x8", "x11", "x12")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x11",
        IVobs = "x10",
        MIIVs = c("x1", "x4", "x7", "x2", "x3", "x5", "x6", "x8", "x9", "x12")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "x12",
        IVobs = "x10",
        MIIVs = c("x1", "x4", "x7", "x2", "x3", "x5", "x6", "x8", "x9", "x11")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    )
  )

test_that("simul_reg_1 is correct", {
  expect_equal_sort(colom_miivs , colom_miivs_test)
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

simul_reg_1_miivs_test <- lapply(miivs(simul_reg_1)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))

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
      MIIVs = c(
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
    .Names = c("DVobs", "IVobs", "MIIVs")
  ), structure(
    list(
      DVobs = "prestg80",
      IVobs = "educ",
      MIIVs = c(
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
               "IVobs", "MIIVs")
  ))

test_that("simul_reg_1 is correct", {
  expect_equal_sort(simul_reg_1_miivs, simul_reg_1_miivs_test)
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

simul_reg_2_miivs_test <- lapply(miivs(simul_reg_2)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))

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
      MIIVs = c(
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
    .Names = c("DVobs", "IVobs", "MIIVs")
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
      MIIVs = c(
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
    .Names = c("DVobs", "IVobs", "MIIVs")
  ))

test_that("simul_reg_2 is correct", {
  expect_equal_sort(simul_reg_2_miivs, simul_reg_2_miivs_test)
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

simul_reg_3_miivs_test <- lapply(miivs(simul_reg_3)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))

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
      MIIVs = c(
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
    .Names = c("DVobs", "IVobs", "MIIVs")
  ), structure(
    list(
      DVobs = "prestg80",
      IVobs = "educ",
      MIIVs = c(
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
               "IVobs", "MIIVs")
  ))
test_that("simul_reg_3 is correct", {
  expect_equal_sort(simul_reg_3_miivs, simul_reg_3_miivs_test)
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

simul_reg_4_miivs_test <- lapply(miivs(simul_reg_4)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))

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
      MIIVs = c(
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
    .Names = c("DVobs", "IVobs", "MIIVs")
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
      MIIVs = c(
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
               "MIIVs")
  ))
test_that("simul_reg_4 is correct", {
  expect_equal_sort(simul_reg_4_miivs, simul_reg_4_miivs_test)
})

##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Two Factors Example 1 Eq-----------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

two_cfa_1 <- '
Systolic  =~ Z2 + 1*Z1 + 1*Z3
Diastolic =~ Z5 + 1*Z4 + 1*Z6
Systolic  ~~ Diastolic
'

two_cfa_1_miivs_test <- lapply(miivs(two_cfa_1)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))

two_cfa_1_miivs <-
  list(
    structure(
      list(
        DVobs = "Z1",
        IVobs = "Z2",
        MIIVs = c("Z3","Z4", "Z6", "Z5")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "Z3",
        IVobs = "Z2",
        MIIVs = c("Z1", "Z4", "Z6", "Z5")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "Z4",
        IVobs = "Z5",
        MIIVs = c("Z1", "Z3", "Z6", "Z2")
      ),
      .Names = c("DVobs","IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "Z6",
        IVobs = "Z5",
        MIIVs = c("Z1","Z3", "Z4", "Z2")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    )
  )

test_that("two_cfa_1 is correct", {
  expect_equal_sort(two_cfa_1_miivs, two_cfa_1_miivs_test)
})

##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Kirby Model 1 ---------------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

kirby_model_1 <-'

F1 =~ V2 + V1 + V3 + V4
F2 =~ V5 + V4 + V6 + V7
F3 =~ V8 + V7 + V9 + V6

F2 ~ F1 
F3 ~ F2 
'

kirby_model_1_miivs_test <- lapply(miivs(kirby_model_1)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))


kirby_model_1_miivs <-
  list(
    structure(
      list(
        DVobs = "V5",
        IVobs = "V2",
        MIIVs = c("V1","V3")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V8",
        IVobs = "V5",
        MIIVs = c("V4", "V1", "V3", "V2")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V4",
        IVobs = c("V2","V5"),
        MIIVs = c("V8", "V6", "V7", "V9", "V1", "V3")
      ),
      .Names = c("DVobs","IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V6",
        IVobs = c("V5",
                  "V8"),
        MIIVs = c("V4", "V7", "V9", "V1", "V3", "V2")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V7",
        IVobs = c("V5",
                  "V8"),
        MIIVs = c("V4", "V6", "V9", "V1", "V3", "V2")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V9",
        IVobs = "V8",
        MIIVs = c("V5", "V4", "V6", "V7", "V1", "V3", "V2")
      ),
      .Names = c("DVobs", "IVobs","MIIVs")
    ),
    structure(
      list(
        DVobs = "V1",
        IVobs = "V2",
        MIIVs = c("V5","V8", "V4", "V6", "V7", "V9", "V3")
      ),
      .Names = c("DVobs", "IVobs","MIIVs")
    ),
    structure(
      list(
        DVobs = "V3",
        IVobs = "V2",
        MIIVs = c("V5","V8", "V4", "V6", "V7", "V9", "V1")
      ),
      .Names = c("DVobs", "IVobs","MIIVs")
    )
  )

test_that("kirby_model_1 is correct", {
  expect_equal_sort(kirby_model_1_miivs, kirby_model_1_miivs_test)
})



##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Kirby Model 3 ---------------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

kirby_model_3 <-'

F1 =~ V2 + V1 + V3 + V4
F2 =~ V5 + V4 + V6 + V7
F3 =~ V8 + V7 + V9 + V6

F1 ~ V10 + V11 + V12 + V13
F2 ~ F1 + V10 + V12
F3 ~ F2 + V10 + V12

#V10 ~~ V11 + V12 + V13
#V11 ~~ V12 + V13
#V12 ~~ V13
'

kirby_model_3_miivs_test <- lapply(miivs(kirby_model_3)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))


kirby_model_3_miivs <-
  list(
    structure(
      list(
        DVobs = "V2",
        IVobs = c("V10", "V11", "V12",
                  "V13"),
        MIIVs = c("V10", "V11", "V12", "V13")
      ),
      .Names = c("DVobs",
                 "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V5",
        IVobs = c("V2",
                  "V10", "V12"),
        MIIVs = c("V10", "V11", "V12", "V13", "V1", "V3")
      ),
      .Names = c("DVobs",
                 "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V8",
        IVobs = c("V5",
                  "V10", "V12"),
        MIIVs = c("V10", "V11", "V12", "V13", "V2", "V1",
               "V3", "V4")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V1",
        IVobs = "V2",
        MIIVs = c("V10", "V11", "V12", "V13",
               "V5", "V8", "V3", "V4", "V6", "V7", "V9")
      ),
      .Names = c("DVobs",
                 "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V3",
        IVobs = "V2",
        MIIVs = c("V10",
               "V11", "V12", "V13", "V5", "V8", "V1", "V4", "V6", "V7", "V9")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V4",
        IVobs = c("V2", "V5"),
        MIIVs = c("V10", "V11", "V12", "V13",
               "V8", "V1", "V3", "V6", "V7", "V9")
      ),
      .Names = c("DVobs",
                 "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V6",
        IVobs = c("V5",
                  "V8"),
        MIIVs = c("V10", "V11", "V12", "V13", "V2", "V1", "V3", "V4",
               "V7", "V9")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V7",
        IVobs = c("V5", "V8"),
        MIIVs = c("V10", "V11",
               "V12", "V13", "V2", "V1", "V3", "V4", "V6", "V9")
      ),
      .Names = c("DVobs",
                 "IVobs", "MIIVs")
    ),
    structure(
      list(
        DVobs = "V9",
        IVobs = "V8",
        MIIVs = c("V10",
               "V11", "V12", "V13", "V2", "V5", "V1", "V3", "V4", "V6", "V7")
      ),
      .Names = c("DVobs", "IVobs", "MIIVs")
    )
  )


test_that("kirby_model_3 is correct", {
  expect_equal_sort(kirby_model_3_miivs, kirby_model_3_miivs_test)
})


##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Trust 1 - GMM ---------------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

trust_model_1 <-'
F1 =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8
'

trust_model_1_miivs_test <- lapply(miivs(trust_model_1)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))


trust_model_1_miivs <-
  list(
    structure(
      list(
        DVobs = "V2",
        IVobs = "V1",
        MIIVs = c("V3",
               "V4", "V5", "V6", "V7", "V8")
      ),
      .Names = c("DVobs", "IVobs",
                 "MIIVs")
    ),
    structure(
      list(
        DVobs = "V3",
        IVobs = "V1",
        MIIVs = c("V2",
               "V4", "V5", "V6", "V7", "V8")
      ),
      .Names = c("DVobs", "IVobs",
                 "MIIVs")
    ),
    structure(
      list(
        DVobs = "V4",
        IVobs = "V1",
        MIIVs = c("V2",
               "V3", "V5", "V6", "V7", "V8")
      ),
      .Names = c("DVobs", "IVobs",
                 "MIIVs")
    ),
    structure(
      list(
        DVobs = "V5",
        IVobs = "V1",
        MIIVs = c("V2",
               "V3", "V4", "V6", "V7", "V8")
      ),
      .Names = c("DVobs", "IVobs",
                 "MIIVs")
    ),
    structure(
      list(
        DVobs = "V6",
        IVobs = "V1",
        MIIVs = c("V2",
               "V3", "V4", "V5", "V7", "V8")
      ),
      .Names = c("DVobs", "IVobs",
                 "MIIVs")
    ),
    structure(
      list(
        DVobs = "V7",
        IVobs = "V1",
        MIIVs = c("V2",
               "V3", "V4", "V5", "V6", "V8")
      ),
      .Names = c("DVobs", "IVobs",
                 "MIIVs")
    ),
    structure(
      list(
        DVobs = "V8",
        IVobs = "V1",
        MIIVs = c("V2",
               "V3", "V4", "V5", "V6", "V7")
      ),
      .Names = c("DVobs", "IVobs",
                 "MIIVs")
    )
  )

test_that("trust_model_1 is correct", {
  expect_equal_sort(trust_model_1_miivs, trust_model_1_miivs_test)
})

##-----------------------------------------------------------##
##-----------------------------------------------------------##
##------------------Trust 2 - GMM ---------------------------##
##-----------------------------------------------------------##
##-----------------------------------------------------------##

trust_model_2 <-'
F1 =~ y1 + y2 
F1 ~  x1 + x2 + x3 + x4 + x5 + x6
'

trust_model_2_miivs_test <- lapply(miivs(trust_model_2)$eqns, "[", c("DVobs", "IVobs", "MIIVs"))

trust_model_2_miivs <-
  list(structure(list(DVobs = "y1", 
                      IVobs = c("x1", "x2", "x3", "x4", "x5", "x6"), 
                      MIIVs = c("x1", "x2", "x3", "x4", "x5", "x6")), 
                 .Names = c("DVobs",  "IVobs", "MIIVs")), 
       structure(list(DVobs = "y2", 
                      IVobs = "y1", 
                      MIIVs = c("x1","x2", "x3", "x4", "x5", "x6")), 
                 .Names = c("DVobs", "IVobs", "MIIVs")))

test_that("trust_model_2 is correct", {
  expect_equal_sort(trust_model_2_miivs, trust_model_2_miivs_test)
})

