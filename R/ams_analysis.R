
#' Performs simultations of AMS statistics analysis
#'
#' @param methods A character vector. The names of the AMS analysis methods
#'     to be simulated.
#' @param suscept.matrix Numeric 3x3 symmetric matrix.
#' @param nmeasures Integer. Number of repetitions of the simulated 
#'     experimental measures.
#' @param error.dist A function for simulating the experimental errors,
#'     as in \code{FakeMeasures}.
#' @param setup Object of class \code{AMSsetup}.
#' @param ... Extra arguments to be passed to the functions given in the
#'     \code{methods} argument.
#'
#' @export
#' @examples
#' # Initial susceptibility tensor to be used as reference for the generation
#' # of simulated data
#' suscep0 <- matrix(c(1,0,0, 0,2,0, 0,0,3), nrow = 3)
#' 
#' AMSsimulations(methods = c('ams.hext', 'ams.constable'), 
#'              suscept.matrix = suscep0,
#'              nmeasures = 10, 
#'              error.dist = NormalErrorGenerator(0.3), 
#'              setup = AMSsetup())
AMSsimulations <- function(methods, suscept.matrix,
                         nmeasures, error.dist, setup=AMSsetup(), ...){

    extra_args_names <- names(list(...))

    measures <- FakeMeasures(suscept.matrix, nmeasures, error.dist)
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
#' @param suscept.matrix Numeric 3x3 symmetric matrix.
#' @param nmeasures Integer. Number of repetitions of the simulated 
#'     experimental measures.
#' @param error.dist A function for simulating the experimental errors,
#'     as in \code{FakeMeasures}.
#' @param nrepetitions Integer. The number of the repetitions of the whole
#'     simulation process, needed to estimate the consistency of the 
#'     AMS analysis methods under study.
#' @param setup Object of class \code{AMSsetup}.
#' @param ... Extra arguments to be passed to the functions given in the
#'     \code{methods} argument.
#'
#' @export
#' @examples
#' # Initial susceptibility tensor to be used as reference for the generation
#' # of simulated data
#' suscep0 <- matrix(c(1,0,0, 0,2,0, 0,0,3), nrow = 3)
#' 
#' AnalyseConsistency(methods = c('ams.hext', 'ams.constable'), 
#'                  suscept.matrix = suscep0,
#'                  nmeasures = 5, 
#'                  error.dist = NormalErrorGenerator(0.3), 
#'                  nrepetitions = 10, 
#'                  setup = AMSsetup(),
#'                  R = 100 # additional parameter for ams.constable method`
#'                  )
AnalyseConsistency <- function(methods, 
                             suscept.matrix,
                             nmeasures,error.dist,nrepetitions, 
                             setup=AMSsetup(), ...){

    param <- eigen(suscept.matrix)

    dir_ini <- param$vectors
    eigen_ini <- param$values
    #print(dir_ini)
    #print(eigen_ini)

    initstats <- function() {
        list(totalcount=0, consistent.eigenvalues=c(0,0,0), consistent.eigenvectors=c(0,0,0), sum.errors.eigenvalues=c(0,0,0), sum.errors.etas=c(0,0,0), sum.errors.zetas=c(0,0,0), reject.oblate=0, reject.prolate=0, reject.isotropy=0)
    }

    updatestats <- function(stats, test.results) {
        stats$totalcount <- stats$totalcount + 1

        eigenvalues <- eigenvalues(test.results)
        eigenvectors <- eigenvectors(test.results)
        anisotropy.test <- anisotropytest(test.results)

        stats$sum.errors.eigenvalues <- stats$sum.errors.eigenvalues + (eigenvalues$upper.limits - eigenvalues$lower.limits)/2
        stats$sum.errors.etas[1] <- stats$sum.errors.etas[1] + eigenvectors$ellip1$eta
        stats$sum.errors.etas[2] <- stats$sum.errors.etas[2] + eigenvectors$ellip2$eta
        stats$sum.errors.etas[3] <- stats$sum.errors.etas[3] + eigenvectors$ellip3$eta

        stats$sum.errors.zetas[1] <- stats$sum.errors.zetas[1] + eigenvectors$ellip1$zeta
        stats$sum.errors.zetas[2] <- stats$sum.errors.zetas[2] + eigenvectors$ellip2$zeta
        stats$sum.errors.zetas[3] <- stats$sum.errors.zetas[3] + eigenvectors$ellip3$zeta

        if (anisotropy.test[kRejectOblate]) {
            stats$reject.oblate = stats$reject.oblate + 1
        }
        if (anisotropy.test[kRejectProlate]) {
            stats$reject.prolate = stats$reject.prolate + 1
        }
        if (anisotropy.test[kRejectIsotropy]) {
            stats$reject.isotropy = stats$reject.isotropy + 1
        }

        if (.__isConsist_tau1(eigen_ini[1], eigenvalues)) {
            stats$consistent.eigenvalues[1] <- stats$consistent.eigenvalues[1] + 1
        }
        if (.__isConsist_tau2(eigen_ini[2], eigenvalues)) {
            stats$consistent.eigenvalues[2] <- stats$consistent.eigenvalues[2] + 1
        }
        if (.__isConsist_tau3(eigen_ini[3], eigenvalues)) {
            stats$consistent.eigenvalues[3] <- stats$consistent.eigenvalues[3] + 1
        }
        if (.__isConsist_vec1(car2sph(t(dir_ini[,1])), eigenvectors)) {
            stats$consistent.eigenvectors[1] <- stats$consistent.eigenvectors[1] + 1
        }
        if (.__isConsist_vec2(car2sph(t(dir_ini[,2])), eigenvectors)) {
            stats$consistent.eigenvectors[2] <- stats$consistent.eigenvectors[2] + 1
        }
        if (.__isConsist_vec3(car2sph(t(dir_ini[,3])), eigenvectors)) {
            stats$consistent.eigenvectors[3] <- stats$consistent.eigenvectors[3] + 1
        }
        return(stats)
    }

    closestats <- function(stats) {
        results <- with(stats, 
                        list(stats.detail=stats,
                             fraction.consist.eigenvalues=consistent.eigenvalues/totalcount, 
                             fraction.consist.eigenvectors=consistent.eigenvectors/totalcount, 
                             mean.error.eigenvalues=sum.errors.eigenvalues/totalcount, 
                             mean.error.eta_eigenvectors=sum.errors.etas/totalcount, 
                             mean.error.zeta_eigenvectors=sum.errors.zetas/totalcount, 
                             fraction.rejected.oblate=reject.oblate/totalcount, 
                             fraction.rejected.prolate=reject.prolate/totalcount, 
                             fraction.rejected.isotropy=reject.isotropy/totalcount))
        return(results)
    }

    results <- list()

    for (funname in methods) {
        results[[funname]] <- initstats()
    }
 
    for (i in 1:nrepetitions) {

        #tesitos <- AMSsimulations(methods, lamb_ini, rot_vec, nmeasures, error.dist, setup, ...)
        tesitos <- AMSsimulations(methods, suscept.matrix, nmeasures, error.dist, setup, ...)

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

   with(eigenvalues, lower.limits[1] <= lamb_ini1 && upper.limits[1] >= lamb_ini1)

}

.__isConsist_tau2 <- function(lamb_ini2,eigenvalues){

    with(eigenvalues, lower.limits[2] <= lamb_ini2 && upper.limits[2] >= lamb_ini2)

}

.__isConsist_tau3 <- function(lamb_ini3,eigenvalues){

    with(eigenvalues, lower.limits[3] <= lamb_ini3 && upper.limits[3] >= lamb_ini3)

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

