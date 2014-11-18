context("AMS simulations")

set.seed(1000001)

test_that("error_norm_dist_generator works", {
          real_vals = c(0, 0.5, 0.44, -1.3, 7.7)
          dist0 <- error_norm_dist_generator(0)
          expect_that(dist0(real_vals, 1)[1,], equals(real_vals))

          results02 <- dist0(real_vals, 2)
          expect_that(results02[1, ], equals(real_vals))
          expect_that(results02[2, ], equals(real_vals))

          check_sigma <- function(sigma) {
              dist <- error_norm_dist_generator(sigma)
              results <- dist(real_vals, 100000)
              means <- apply(X = results, MARGIN = 2, FUN = mean)
              expect_that(means, equals(real_vals, tolerance=0.01))

              sds <- apply(X = results, MARGIN = 2, FUN = sd)
              expect_that(sds, equals(rep(sigma, length(real_vals)), tolerance=0.01))
          }

          check_sigma(0.1)
          check_sigma(0.5)
          check_sigma(5)
          check_sigma(0)
          check_sigma(2.3)
})

test_that("fake_measurements works", {

          check_fake <- function(tensor, setup) {
              dist0 <- error_norm_dist_generator(0)
              fake0 <- fake_measurements(tensor, 1, dist0, setup)

              expect_that(fake0, equals(AMSmeasures.exact(tensor, setup)))
              expect_that(SuscTensor(fake0, setup), equals(tensor))

              fake0_10 <- fake_measurements(tensor, 10, dist0, setup)
              expect_that(SuscTensor(fake0_10, setup), equals(tensor))

              dist0.2 <- error_norm_dist_generator(0.2)
              fake0.2 <- fake_measurements(tensor, 500, dist0.2, setup)
              expect_that(SuscTensor(fake0.2, setup), equals(tensor, tolerance = 0.05))
          }

          setup <- AMSsetup()

          tensor1 <- suscep_matrix(c(1, 1.3, 1.2), c(pi/11, pi/3, pi/2))
          check_fake(tensor1, setup)
          check_fake(tensor1, setup)
          check_fake(tensor1, setup)
          check_fake(tensor1, setup)
          check_fake(tensor1, setup)

          tensor2 <- vector2symtensor(c(-.5, -1, 19, 3.4, -5, 6))
          check_fake(tensor2, setup)
          check_fake(tensor2, setup)
          check_fake(tensor2, setup)
          check_fake(tensor2, setup)
          check_fake(tensor2, setup)
})
