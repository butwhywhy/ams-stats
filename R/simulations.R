#source('utils.R')
#source('ams.R')

#' Generates a synthetic AMS measurements object
#'
#' @param sus_matrix Nueric 3x3 symetric matrix.
#' @param n_measurements Integer. The number of repetitions of the 
#'     simulated experimental protocol.
#' @param error_dist A function accepting two arguments and returning a matrix.
#'     It must accept a numeric vector as the first argument and an integer
#'     as the second one, and return an \code{n} times \code{m} numeric
#'     matrix, being \code{n} the length of its first argument and \code{m}
#'     the second one. This function is interpreted as a random generator
#'     for the simulated AMS measures.
#' @param setup Object of class \code{ams.setup}.
#'
#' @export
#' @examples
#' # Reference ('real') susceptibility matrix
#' sus_matrix <- matrix(c(1.6,.5,.3, .5,1.3,.6, .3,.6,1.9), nrow = 3)
#' sus_matrix
#'
#' # Simulated error distribution with standard deviation 0.5
#' error_dist <- error_norm_dist_generator(0.5)
#' # Experimental setup object
#' setup <- ams.setup()
#'
#' # Fake AMS measurements object with 2 repetitions
#' fake_measurements(sus_matrix, 2, error_dist, setup)
#'
fake_measurements <- function(sus_matrix, n_measurements, error_dist, setup=ams.setup()) {
    # error_dist(real_k, n_measures) takes a vector real_k of real values and returns a matrix 
    # with n_measures rows and length(real_k) columns, where every column containes noisy 
    # measures from a single real measure

    setup_matrix <- design_matrix(setup)
    position_number <- nrow(setup_matrix)

    # Susceptibility column vector with the 6 independent coefficients

    S <- symtensor2vector(matrix(sus_matrix, nrow=3))

    # Vector with the position_number's positions from Jelinek's configuration

    real_k <- setup_matrix %*% S

    # BSus' column from dataframe 

    k_fake <- as.vector(error_dist(real_k, n_measurements))

    # N's column from dataframe

    N_vec <- c(1:n_measurements)

    N_column <- rep(N_vec,position_number)

    # Specimen's column from dataframe

    Specimen <- rep(c(1:position_number),each=n_measurements)

    return(ams.measures(repetitions=N_column, positions=Specimen, values=k_fake))
}

#' Generates a normal noise generator for use in \code{fake_measurements} function.
#'
#' @param sigma Numeric. The standard deviation of the simulated normal noise.
#'
#' @export
#' @examples
#' # Normal error distribution function for synthetic AMS measures with 
#' # standard debiation 0.1
#' error_dist <- error_norm_dist_generator(0.1)
error_norm_dist_generator <- function(sigma) {
    error_norm <- function(real_k, n_measures) {
        position_number <- length(real_k)
        N <- position_number*n_measures
        fake_values <- rnorm(N, real_k, sigma)
        fake_values_mat <- matrix(fake_values,n_measures,position_number,byrow=T)
        return(fake_values_mat)
    }
    return(error_norm)
}

