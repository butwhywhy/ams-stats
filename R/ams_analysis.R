#source('ams.R')
#source('simulations.R')
#source('ellipses.R')

#' Performs simultations of AMS statistics analysis
#'
#' @param methods A character vector. The names of the AMS analysis methods
#'     to be simulated.
#' @param sus_matrix Numeric 3x3 symmetric matrix.
#' @param n_measurements Integer. Number of repetitions of the simulated 
#'     experimental measures.
#' @param error_dist A function for simulating the experimental errors,
#'     as in \code{fake_measurements}.
#' @param setup Object of class \code{ams.setup}.
#' @param ... Extra arguments to be passed to the functions given in the
#'     \code{methods} argument.
#'
#' @export
#' @examples
#' # Initial susceptibility tensor to be used as reference for the generation
#' # of simulated data
#' suscep0 <- matrix(c(1,0,0, 0,2,0, 0,0,3), nrow = 3)
#' 
#' sim_analysis(methods = c('ams.hext.woN', 'ams.constable'), 
#'              sus_matrix = suscep0,
#'              n_measurements = 10, 
#'              error_dist = error_norm_dist_generator(0.3), 
#'              setup = ams.setup())
sim_analysis <- function(methods, sus_matrix,
                         n_measurements, error_dist, setup=ams.setup(), ...){

    extra_args_names <- names(list(...))

    measures <- fake_measurements(sus_matrix, n_measurements, error_dist)
    result <- list()
    for (funname in methods) {
        fun <- get(funname)

        # Select appropiate parameters for current method
        admited_args <- names(formals(fun))
        extra_args <- list(...)
        if (! ('...' %in% admited_args)) {
            extra_args <- extra_args[extra_args_names %in% admited_args]
        }
        all_args <- append(list(setup=setup, measures=measures), extra_args)

        analysis <- do.call(fun, all_args)
        result[[funname]] <- analysis
    }

    return(result)

}

#' Studies the consistency of AMS analysis methods
#' using syntetic data.
#'
#' @param methods A character vector. The names of the AMS analysis methods
#'     to be simulated.
#' @param sus_matrix Numeric 3x3 symmetric matrix.
#' @param n_measurements Integer. Number of repetitions of the simulated 
#'     experimental measures.
#' @param error_dist A function for simulating the experimental errors,
#'     as in \code{fake_measurements}.
#' @param m_iterations Integer. The number of the repetitions of the whole
#'     simulation process, needed to estimate the consistency of the 
#'     AMS analysis methods under study.
#' @param setup Object of class \code{ams.setup}.
#' @param ... Extra arguments to be passed to the functions given in the
#'     \code{methods} argument.
#'
#' @export
#' @examples
#' # Initial susceptibility tensor to be used as reference for the generation
#' # of simulated data
#' suscep0 <- matrix(c(1,0,0, 0,2,0, 0,0,3), nrow = 3)
#' 
#' consist_analysis(methods = c('ams.hext.woN', 'ams.constable'), 
#'                  sus_matrix = suscep0,
#'                  n_measurements = 5, 
#'                  error_dist = error_norm_dist_generator(0.3), 
#'                  m_iterations = 10, 
#'                  setup = ams.setup(),
#'                  R = 100 # additional parameter for ams.constable method`
#'                  )
consist_analysis <- function(methods, 
                             sus_matrix,
                             n_measurements,error_dist,m_iterations, 
                             setup=ams.setup(), ...){

    param <- eigen(sus_matrix)

    dir_ini <- param$vectors
    eigen_ini <- param$values
    #print(dir_ini)
    #print(eigen_ini)

    initstats <- function() {
        list(totalcount=0, consistent_eigenvalues=c(0,0,0), consistent_eigenvectors=c(0,0,0), sum_errors_eigenvalues=c(0,0,0), sum_errors_etas=c(0,0,0), sum_errors_zetas=c(0,0,0), reject1Eq2=0, reject2Eq3=0, reject1Eq2Eq3=0)
    }

    updatestats <- function(stats, test_results) {
        stats$totalcount <- stats$totalcount + 1

        eigenvalues <- eigenvalues.ams.analysis(test_results)
        eigenvectors <- eigenvectors.ams.analysis(test_results)
        anisotropy_test <- anisotropy_test.ams.analysis(test_results)

        stats$sum_errors_eigenvalues <- stats$sum_errors_eigenvalues + (eigenvalues$taus_high - eigenvalues$taus_low)/2
        stats$sum_errors_etas[1] <- stats$sum_errors_etas[1] + eigenvectors$ellip1$eta
        stats$sum_errors_etas[2] <- stats$sum_errors_etas[2] + eigenvectors$ellip2$eta
        stats$sum_errors_etas[3] <- stats$sum_errors_etas[3] + eigenvectors$ellip3$eta

        stats$sum_errors_zetas[1] <- stats$sum_errors_zetas[1] + eigenvectors$ellip1$zeta
        stats$sum_errors_zetas[2] <- stats$sum_errors_zetas[2] + eigenvectors$ellip2$zeta
        stats$sum_errors_zetas[3] <- stats$sum_errors_zetas[3] + eigenvectors$ellip3$zeta

        if (reject1Eq2(anisotropy_test)) {
            stats$reject1Eq2 = stats$reject1Eq2 + 1
        }
        if (reject2Eq3(anisotropy_test)) {
            stats$reject2Eq3 = stats$reject2Eq3 + 1
        }
        if (reject1Eq2Eq3(anisotropy_test)) {
            stats$reject1Eq2Eq3 = stats$reject1Eq2Eq3 + 1
        }

        if (.__isConsist_tau1(eigen_ini[1], eigenvalues)) {
            stats$consistent_eigenvalues[1] <- stats$consistent_eigenvalues[1] + 1
        }
        if (.__isConsist_tau2(eigen_ini[2], eigenvalues)) {
            stats$consistent_eigenvalues[2] <- stats$consistent_eigenvalues[2] + 1
        }
        if (.__isConsist_tau3(eigen_ini[3], eigenvalues)) {
            stats$consistent_eigenvalues[3] <- stats$consistent_eigenvalues[3] + 1
        }
        if (.__isConsist_vec1(car2sph(t(dir_ini[,1])), eigenvectors)) {
            stats$consistent_eigenvectors[1] <- stats$consistent_eigenvectors[1] + 1
        }
        if (.__isConsist_vec2(car2sph(t(dir_ini[,2])), eigenvectors)) {
            stats$consistent_eigenvectors[2] <- stats$consistent_eigenvectors[2] + 1
        }
        if (.__isConsist_vec3(car2sph(t(dir_ini[,3])), eigenvectors)) {
            stats$consistent_eigenvectors[3] <- stats$consistent_eigenvectors[3] + 1
        }
        return(stats)
    }

    closestats <- function(stats) {
        results <- with(stats, 
                        list(detailed_stats=stats,
                             perc_consist_eigenvalues=consistent_eigenvalues/totalcount, 
                             perc_consist_eigenvectors=consistent_eigenvectors/totalcount, 
                             mean_error_eigenvalues=sum_errors_eigenvalues/totalcount, 
                             mean_error_eta_eigenvectors=sum_errors_etas/totalcount, 
                             mean_error_zeta_eigenvectors=sum_errors_zetas/totalcount, 
                             perc_rejected_1Eq2=reject1Eq2/totalcount, 
                             perc_rejected_2Eq3=reject2Eq3/totalcount, 
                             perc_rejected_1Eq2Eq3=reject1Eq2Eq3/totalcount))
        return(results)
    }

    results <- list()

    for (funname in methods) {
        results[[funname]] <- initstats()
    }
 
    for (i in 1:m_iterations) {

        #tesitos <- sim_analysis(methods, lamb_ini, rot_vec, n_measurements, error_dist, setup, ...)
        tesitos <- sim_analysis(methods, sus_matrix, n_measurements, error_dist, setup, ...)

        for (funname in methods) {
            results[[funname]] <- updatestats(results[[funname]], tesitos[[funname]])
        }

    }

    for (funname in methods) {
        results[[funname]] <- closestats(results[[funname]])
    }
    return(results)
}

.__isConsist_tau1 <- function(lamb_ini1,eigenvalues){

   with(eigenvalues, taus_low[1] <= lamb_ini1 && taus_high[1] >= lamb_ini1)

}

.__isConsist_tau2 <- function(lamb_ini2,eigenvalues){

    with(eigenvalues, taus_low[2] <= lamb_ini2 && taus_high[2] >= lamb_ini2)

}

.__isConsist_tau3 <- function(lamb_ini3,eigenvalues){

    with(eigenvalues, taus_low[3] <= lamb_ini3 && taus_high[3] >= lamb_ini3)

}

.__isConsist_vec1 <- function(dir_ini1,eigenvectors) {

    if (! is.matrix(dir_ini1)) {
        dir_ini1 <- t(dir_ini1)
    }
    return (in.ellipse(eigenvectors$ellip1, dir_ini1, TRUE))
}

.__isConsist_vec2 <- function(dir_ini2,eigenvectors) {

    if (! is.matrix(dir_ini2)) {
        dir_ini2 <- t(dir_ini2)
    }
    return (in.ellipse(eigenvectors$ellip2, dir_ini2, TRUE))
}

.__isConsist_vec3 <- function(dir_ini3,eigenvectors) {

    if (! is.matrix(dir_ini3)) {
        dir_ini3 <- t(dir_ini3)
    }
    return (in.ellipse(eigenvectors$ellip3, dir_ini3, TRUE))
}
