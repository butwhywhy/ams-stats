context("AMS base functions")

test_that("I know how to test", {
          expect_that(2, equals(1 + 1))
})

# ----------------------------------------
# Auxiliar data for tests
vec_Aux<-c(.5, .5, 0, -1, 0, 0,
           .5, .5, 0, 1, 0, 0,
           1, 0, 0, 0, 0, 0,
           .5, .5, 0, -1, 0, 0,
           .5, .5, 0, 1, 0, 0,
           0, .5, .5, 0, -1, 0,
           0, .5, .5, 0, 1, 0,
           0, 1, 0, 0, 0, 0,
           0, .5, .5, 0, -1, 0,
           0, .5, .5, 0, 1, 0,
           .5, 0, .5, 0, 0, -1,
           .5, 0, .5, 0, 0, 1,
           0, 0, 1, 0, 0, 0,
           .5, 0, .5, 0, 0, -1,
           .5, 0, .5, 0, 0, 1
           )
vec_Aux6x6 <- c(1, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0,
                .5, .5, 0, 1, 0, 0,
                0, .5, .5, 0, 1, 0,
                .5, 0, .5, 0, 0, 1
                )
mat_D_6x6 <-matrix(vec_Aux6x6, 6, 6, byrow=T)
mat_D_default <-matrix(vec_Aux, 15, 6, byrow=T)
# ----------------------------------------

test_that("AMSsetup constructor works", {
          default <- AMSsetup()
          expect_that(default, is_a('AMSsetup'))
          expect_that(design(default), is_identical_to(mat_D_default))

          setup6x6 <- AMSsetup(setup.matrix = mat_D_6x6)
          expect_that(setup6x6, is_a('AMSsetup'))
          expect_that(design(setup6x6), is_identical_to(mat_D_6x6))
})

test_that("AMSmeasures constructors works", {
          reps <- rep(1, 15)
          poss <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
          values <- c(0.86,1.26,1.10,0.86,1.26,1.26,1.26,1.02,1.26,1.26,
                      0.80,1.80,1.50,0.80,1.80)
          measures <- AMSmeasures(reps, poss, values)
          expect_that(measures, is_a('AMSmeasures'))
          
          dataframe <- data.frame(reps, poss, values)
          measures2 <- as.AMSmeasures(dataframe)
          expect_that(measures2, is_a('AMSmeasures'))
          expect_that(measures2, equals(measures))

})

test_that("estimation of AMS tensor and generation of exact measures are consistent", {
          check_consistent <- function(tensor, setup=AMSsetup()) {
              measures <- ExactMeasures(tensor, setup)
              tensor2 <- SuscTensor(measures, setup)
              expect_that(tensor2, equals(tensor))
          }
          
          check_consistent(suscep_matrix(c(1, 1.1, 1.9)))
          check_consistent(suscep_matrix(c(1.9, 1, 1.9)))
          check_consistent(suscep_matrix(c(2, 2, 2)))
          check_consistent(suscep_matrix(c(0, 0, 2)))
          check_consistent(suscep_matrix(c(-1, 1.1, 1.9)))
          check_consistent(suscep_matrix(c(1, 1.1, 1.9), c(pi/3, pi/8, pi/7)))
          check_consistent(suscep_matrix(c(1, 1.1, 1.9), c(pi/3, pi/8, pi/7)))
          check_consistent(suscep_matrix(c(1, 1, 1), c(pi/3, pi/8, pi/7)))
})

