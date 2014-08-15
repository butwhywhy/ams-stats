context("AMS analysis methods")

load('backwards-data.rda')

alpha <- 0.9

test_that("Bootstrap method is backwards compatible", {
          R <- 100
          setup <- ams.setup()

          for (i in 1:3) {
              for (j in 1:4) {
                  measure_name <- paste('measures', i, '_', j, sep = '')
                  measures <- get(measure_name)

                  # Bootstrap analysis
                  set.seed(1000001)
                  boot_analysis_name <- paste('constable', i, '_', j, sep = '')
                  boot_analysis <- ams.constable(measures, setup, 
                                                 alpha = alpha, R = R)
                  expect_that(boot_analysis, equals(get(boot_analysis_name)))

              }
          }
})

test_that("Hext method is backwards compatible", {
          setup <- ams.setup()

          for (i in 1:3) {
              for (j in 1:4) {
                  measure_name <- paste('measures', i, '_', j, sep = '')
                  measures <- get(measure_name)

                  # Hext analysis
                  set.seed(1000001)
                  hext_analysis_name <- paste('hext', i, '_', j, sep = '')
                  hext_analysis <- ams.hext(measures, setup, 
                                                 alpha = alpha)
                  expect_that(hext_analysis, equals(get(hext_analysis_name)))

                  # Hext without N analysis
                  set.seed(1000001)
                  hextwon_analysis_name <- paste('hextwon', i, '_', j, sep = '')
                  hextwon_analysis <- ams.hext.woN(measures, setup, 
                                                 alpha = alpha)
                  expect_that(hextwon_analysis, equals(get(hextwon_analysis_name)))

              }
          }
})

