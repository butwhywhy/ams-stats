context("AMS analysis method: Hext")

load('backwards-data.rda')

alpha <- 0.9

test_that("Hext method is backwards compatible", {
          setup <- AMSsetup()

          for (i in 1:3) {
              for (j in 1:4) {
                  measure_name <- paste('measures', i, '_', j, sep = '')
                  measures <- get(measure_name)

                  # Hext analysis
                  set.seed(1000001)
                  hext_analysis_name <- paste('hext', i, '_', j, sep = '')
                  hext_analysis <- HextAnalyse(measures, setup, 
                                                 alpha = alpha)
                  expect_that(hext_analysis, equals(get(hext_analysis_name)))

                  # Hext without N analysis
                  set.seed(1000001)
                  hextwon_analysis_name <- paste('hextwon', i, '_', j, sep = '')
                  hextwon_analysis <- HextAnalyse(measures, setup, 
                                               alpha = alpha, withoutN = TRUE)
                  expect_that(hextwon_analysis, equals(get(hextwon_analysis_name)))

              }
          }
})

