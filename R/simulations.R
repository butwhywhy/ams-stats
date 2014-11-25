
#' Generates a synthetic AMS measurements object
#'
#' @param suscept.matrix Nueric 3x3 symetric matrix.
#' @param nmeasures Integer. The number of repetitions of the 
#'     simulated experimental protocol.
#' @param error.dist A function accepting two arguments and returning a matrix.
#'     It must accept a numeric vector as the first argument and an integer
#'     as the second one, and return an \code{n} times \code{m} numeric
#'     matrix, being \code{n} the length of its first argument and \code{m}
#'     the second one. This function is interpreted as a random generator
#'     for the simulated AMS measures.
#' @param setup Object of class \code{AMSsetup}.
#'
#' @export
#' @examples
#' # Reference ('real') susceptibility matrix
#' suscept.matrix <- matrix(c(1.6,.5,.3, .5,1.3,.6, .3,.6,1.9), nrow = 3)
#' suscept.matrix
#'
#' # Simulated error distribution with standard deviation 0.5
#' error.dist <- NormalErrorGenerator(0.5)
#' # Experimental setup object
#' setup <- AMSsetup()
#'
#' # Fake AMS measurements object with 2 repetitions
#' FakeMeasures(suscept.matrix, 2, error.dist, setup)
#'
FakeMeasures <- function(suscept.matrix, nmeasures, error.dist, setup=AMSsetup()) {
    # error.dist(real_k, n_measures) takes a vector real_k of real values and returns a matrix 
    # with n_measures rows and length(real_k) columns, where every column containes noisy 
    # measures from a single real measure

    setup_matrix <- design(setup)
    position_number <- nrow(setup_matrix)

    # Susceptibility column vector with the 6 independent coefficients

    S <- symtensor2vector(matrix(suscept.matrix, nrow=3))

    # Vector with the position_number's positions from Jelinek's configuration

    real_k <- setup_matrix %*% S

    # BSus' column from dataframe 

    k_fake <- as.vector(error.dist(real_k, nmeasures))

    # N's column from dataframe

    N_vec <- c(1:nmeasures)

    N_column <- rep(N_vec,position_number)

    # Specimen's column from dataframe

    Specimen <- rep(c(1:position_number),each=nmeasures)

    return(AMSmeasures(repetitions=N_column, positions=Specimen, values=k_fake))
}

#' Generates a normal noise generator for use in \code{FakeMeasures} function.
#'
#' @param sigma Numeric. The standard deviation of the simulated normal noise.
#'
#' @export
#' @examples
#' # Normal error distribution function for synthetic AMS measures with 
#' # standard debiation 0.1
#' error.dist <- NormalErrorGenerator(0.1)
NormalErrorGenerator <- function(sigma) {
    error_norm <- function(real_k, n_measures) {
        position_number <- length(real_k)
        N <- position_number*n_measures
        fake_values <- rnorm(N, real_k, sigma)
        fake_values_mat <- matrix(fake_values,n_measures,position_number,byrow=T)
        return(fake_values_mat)
    }
    return(error_norm)
}

