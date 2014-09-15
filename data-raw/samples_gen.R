# Script to generate the data saved at data/samples_measures.rda

source('R/utils.R')
source('R/ams.R')

sus <- suscep_matrix(lamb_vector = c(.67, 1, 1.33), rot_vector = c(pi/4, pi/4, pi/4))

set.seed(1000001)
error_dist <- error_norm_dist_generator(0.33)
setup <- ams.setup()

measures <- fake_measurements(sus, 10, error_dist, setup)

sample_measures <- as.data.frame(measures)
save(sample_measures, file = 'data/sample_measures.rda')
