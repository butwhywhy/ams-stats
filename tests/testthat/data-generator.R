
stop("This script should not be normally executed")

setup <- ams.setup()

nmeas <- 8
alpha <- 0.9
R <- 100

tensor1 <- suscep_matrix(c(1, 1.3, 1.2), c(pi/11, pi/3, pi/2))
tensor2 <- vector2symtensor(c(-.5, -1, 19, 3.4, -5, 6))
tensor3 <- vector2symtensor(c(1.5, 2.1, 4, 30.1, 0.1, 0))

dist1 <- error_norm_dist_generator(0.05)
dist2 <- error_norm_dist_generator(0.5)
dist3 <- error_norm_dist_generator(1.07)
dist4 <- error_norm_dist_generator(3)

tensors <- list(tensor1, tensor2, tensor3)
dists <- list(dist1, dist2, dist3, dist4)

objects_to_save <- c()

for (i in 1:3) {
    for (j in 1:4) {
        # Construct fake measurements
        set.seed(1000001)
        measure_name <- paste('measures', i, '_', j, sep = '')
        measures <- fake_measurements(tensors[[i]], nmeas, dists[[j]], setup)
        assign(measure_name, measures)

        # Bootstrap analysis
        set.seed(1000001)
        boot_analysis_name <- paste('constable', i, '_', j, sep = '')
        boot_analysis <- ams.constable(measures, setup, 
                                      alpha = alpha, R = R)
        assign(boot_analysis_name, boot_analysis)

        # Hext without N analysis
        set.seed(1000001)
        hextwon_analysis_name <- paste('hextwon', i, '_', j, sep = '')
        hextwon_analysis <- ams.hext.woN(measures, setup, alpha = alpha)
        assign(hextwon_analysis_name, hextwon_analysis)

        # Hext analysis
        set.seed(1000001)
        hext_analysis_name <- paste('hext', i, '_', j, sep = '')
        hext_analysis <- ams.hext(measures, setup, alpha = alpha)
        assign(hext_analysis_name, hext_analysis)

        objects_to_save <- append(objects_to_save, 
                                  c(measure_name, boot_analysis_name, 
                                    hext_analysis_name, hextwon_analysis_name))
    }
}

save(list=objects_to_save, file='tests/testthat/backwards-data.rda')
