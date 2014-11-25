context("AMS analysis method: Bootstrap")

load('backwards-data.rda')

alpha <- 0.9

test_that("Bootstrap method is backwards compatible", {
          R <- 100
          setup <- AMSsetup()

          for (i in 1:3) {
              for (j in 1:4) {
                  measure_name <- paste('measures', i, '_', j, sep = '')
                  measures <- get(measure_name)

                  # Bootstrap analysis
                  set.seed(1000001)
                  boot_analysis_name <- paste('constable', i, '_', j, sep = '')
                  boot_analysis <- BootstrapAnalyse(measures, setup, 
                                                 alpha = alpha, R = R)
                  expect_that(boot_analysis, equals(get(boot_analysis_name)))

              }
          }
})

