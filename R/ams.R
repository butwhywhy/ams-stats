
kVectorDefaultAux <-c(.5, .5, 0, -1, 0, 0,
              .5, .5, 0, 1, 0, 0,
              1, 0, 0, 0, 0, 0,
              .5, .5, 0, -1, 0, 0,
              .5, .5, 0, 1, 0, 0,
              0, .5, .5, 0, -1, 0,
              0, .5, .5, 0, 1, 0,
              0, 1, 0, 0, 0, 0,
              0, .5, .5, 0, -1, 0,
              0, .5, .5, 0, 1, 0,
              .5, 0, .5, 0, 0, -1,
              .5, 0, .5, 0, 0, 1,
              0, 0, 1, 0, 0, 0,
              .5, 0, .5, 0, 0, -1,
              .5, 0, .5, 0, 0, 1
              )

kVector6x6Aux <- c(1, 0, 0, 0, 0, 0,
                   0, 1, 0, 0, 0, 0,
                   0, 0, 1, 0, 0, 0,
                   .5, .5, 0, 1, 0, 0,
                   0, .5, .5, 0, 1, 0,
                   .5, 0, .5, 0, 0, 1
                   )

kSetupMatrix6x6 <-matrix(kVector6x6Aux, 6, 6, byrow=T)

kSetupMatrixDefault <-matrix(kVectorDefaultAux , 15, 6, byrow=T)

#' Constructs a AMS experiment setup
#'
#' @param setup.matrix Numeric matrix.
#' @export
#' @seealso \code{\link{design.AMSsetup}}
#' @examples
#' setup <- AMSsetup()
#' class(setup)
AMSsetup <- function(setup.matrix=kSetupMatrixDefault) {
    if (ncol(setup.matrix) != 6) {
        stop("Illegal argument: setup matrix must have 6 columns")
    }

    result <- list()
    class(result) <- c('AMSsetup', class(result))
    result$.__matrix <- setup.matrix
    result$.__B <- .__positions_B(setup.matrix)
    return(result)
}

.__positions_B <- function(positions_D) {
    positions_V <- solve(t(positions_D) %*% positions_D)
    positions_V %*% t(positions_D)
}

#' Returns the design matrix for setup object
#'
#' @param setup Setup object
#' @export
design <- function(setup) {
    UseMethod('design')
}

#' @rdname design
#' @export
design.AMSsetup <- function(setup) {
    return(setup$.__matrix)
}

#' Constructs an AMS measurements object
#'
#' This method constructs an object with class \code{AMSmeasures}, which uses
#' a \code{data.frame} internally.
#'
#' @param repetitions Integer. It must have the same cardinality as the
#'     \code{values} and represents the number of repetition of the
#'     correponding measurement.
#' @param positions Integer vector or factor. It must have the same
#'     cardinality as the \code{vaues} parameter, and represents
#'     the order of the corresponding measurement in the experimental
#'     setup.
#' @param values Numberic vector.The measured values.
#' @export
#' @seealso \code{\link{as.AMSmeasures}}; \code{\link{reps.AMSmeasures}}, 
#'     \code{\link{poss.AMSmeasures}}, \code{\link{values.AMSmeasures}}
#' @examples
#' reps <- rep(1, 15)
#' poss <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
#' values <- c(0.86,1.26,1.10,0.86,1.26,1.26,1.26,1.02,1.26,1.26,0.80,1.80,
#'             1.50,0.80,1.80)
#' # Construct AMS measures object
#' AMSmeasures(reps, poss, values)
AMSmeasures <- function(repetitions, positions, values) {
    if (! is.integer(repetitions)) {
        if (! all(repetitions == as.integer(repetitions))) {
            stop("Illegal argument: repetitions must be an integer vector")
        }
    }
    if (! is.integer(positions)) {
        if (! all(positions == as.integer(positions))) {
            stop("Illegal argument: positions must be an integer vector")
        }
    }
    result <- data.frame(N=as.factor(as.integer(repetitions)), 
                         Specimen=as.factor(as.integer(positions)), 
                         BSus=values)
    class(result) <- c('AMSmeasures', class(result))
    return(result)
}

#' Constructs an AMS measurements object.
#'
#' This method constructs an object with class \code{AMSmeasures}, which uses
#' a \code{data.frame} internally. It takes as argument a data frame with at
#' least three columns (list, data.frame, matrix,...).
#'
#' @param data Data frame with at least three columns. Each row represents a 
#'     repetition of a measurement in one of the experimental setup positions. 
#'     The first columns must be integer and is interpreted as the number of 
#'     repetition of the measurement. The second columns, also 
#'     integer, is the order of the corresponding measurement in the experimental
#'     setup. And the third column, numeric, the measured values.
#' @export
#' @examples
#' # Use sample dataframe
#' class(sample_measures)
#' 
#' # Convert data frame to \code{AMSmeasures} object
#' ams_measure <- as.AMSmeasures(sample_measures)
#' class(sample_measures)
#' 
as.AMSmeasures <- function(data) {
    return(AMSmeasures(data[[1]], data[[2]], data[[3]]))
}

SingleMeasure <- function(positions, values) {
    reps <- rep(as.integer(1), length(values))
    AMSmeasures(reps, positions, values)
}

#' Returns a vector of integers indicating which measure or repetition
#' the corresponding measured data corresponds to.
#'
#' @param measures AMS measures object
#' @export
reps <- function(measures) {
    UseMethod('reps')
}

#' @rdname reps
#' @export
reps.AMSmeasures <- function(measures) {
    return(as.integer(as.character(measures$N)))
}

#' Returns a vector of integers indicating which position in the experimental
#' setup the corresponding measured data corresponds to.
#'
#' @param measures AMS measures object
#' @export
poss <- function(measures) {
    UseMethod('poss')
}

#' @rdname poss
#' @export
poss.AMSmeasures <- function(measures) {
    return(as.integer(as.character(measures$Specimen)))
}

#' Returns the numeric vector of the measured data.
#'
#' @param measures AMS measures object
#' @export
values <- function(measures) {
    UseMethod('values')
}

#' @rdname poss
#' @export
values.AMSmeasures <- function(measures) {
    if (! inherits(measures, 'AMSmeasures')) {
        stop("Illegal argument: measures must be of class 'AMSmeasures'")
    }
    return(measures$BSus)
}

# MeasuresFromFile takes a data file and return an AMSmeasures object
MeasuresFromFile <- function(datapath) {
    frame <- read.table(datapath, header=T)
    AMSmeasures(repetitions=frame$N, positions=frame$Specimen, values=frame$BSus)
}

#' Computes the mean of AMS measures.
#' 
#' Returns an \code{\link{AMSmeasures}} object with a single repetition and 
#' the mean value for each measured position.
#'
#' @param x An object with class \code{AMSmeasures}.
#' @param ... Unused parameters
#' @export
mean.AMSmeasures <- function(x, ...){
    # mean.AMSmeasures takes an AMSmeasures object and return another AMSmeasures object
    # with a single repetition and the mean value for each measured position
    mean_aux <- function(data, l) {
        mean.data <- mean(data)
        return(mean.data)
    } 
    means <- yysaapply(x$BSus, x$Specimen, mean_aux)

    poss <- as.integer(as.character(means$factors))

    SingleMeasure(positions=poss, values=means$values)
}


#' Computes the susceptibility tensor
#' 
#' SuscTensor takes a AMSmeasures object and its setup and returns the 
#' best estimate for the susceptibility tensor as a 3 times 3 matrix.
#'
#' @param measures An object of class \code{AMSmeasures}.
#' @param setup An object of class \code{AMSsetup}
#' @export
#' @examples
#' # Experimental setup object
#' setup <- AMSsetup()
#' # Sample AMS measures object
#' measures <- as.AMSmeasures(sample_measures)
#' # Least squares estimation of the susceptibility tensor
#' SuscTensor(measures, setup)
SuscTensor <- function(measures, setup){
    
    if (! inherits(measures, 'AMSmeasures')) {
        stop("Illegal argument: measures must be of class 'AMSmeasures'")
    }
    if (! inherits(setup, 'AMSsetup')) {
        stop("Illegal argument: setup must be of class 'AMSsetup'")
    }

    D_setup <- design(setup)
    pos <- poss.AMSmeasures(measures);

    D <- .__repeatD(D_setup, pos)

    B <- .__positions_B(D)
    sus_vec <- B %*% values.AMSmeasures(measures)

    sus_mat <- vector2symtensor(sus_vec)

    return(sus_mat)
}

.__repeatD <- function(D_setup, positions) {
    n <- length(positions)
    D <- matrix(rep(0, n*6), ncol=6)
    for (i in 1:n) {
        D[i,] <- D_setup[positions[i],]
    }
    return(D)
}


ExactMeasures <- function(suscept.matrix, setup, positions = NULL) {
    D_setup <- design(setup)
    if (is.null(positions)) {
        positions <- 1:nrow(D_setup)
    }

    counts <- rep(0, nrow(D_setup))
    reps <- rep(0, length(positions))
    j <- 0
    for (i in positions) {
        j <- j + 1
        counts[i] <- counts[i] + 1
        reps[j] <- counts[i]
    }

    D <- .__repeatD(D_setup, positions)
    vec <- symtensor2vector(suscept.matrix)

    values <- D %*% vec
    return(AMSmeasures(values=values, positions=positions, repetitions=as.integer(reps)))
}

AMSanalysis <- function(eigenvalues, eigenvectors, anisotropy.test) {
    if (!inherits(eigenvalues, 'EigenvaluesCI')) {
        stop('eigenvalues must be of class EigenvaluesCI')
    }
    if (!inherits(eigenvectors, 'EigenvectorsCR')) {
        stop('eigenvectors must be of class EigenvectorsCR')
    }
    if (!inherits(anisotropy.test, 'AnisotropyTest')) {
        stop('anisoropy_test must be of class AnisotropyTest')
    }
    results <- list(eigenvalues=eigenvalues, eigenvectors=eigenvectors, anisotropy_test=anisotropy.test)
    class(results) <- c('AMSanalysis', class(results))
    return(results)
}

EigenvaluesCI <- function(eigenvalues, errors, lower.limits=NULL, upper.limits=NULL) {
    # Eigenvalues confidence intervals
    if (is.null(errors) && (is.null(lower.limits) || is.null(upper.limits))) {
        stop('if lower.limits or upper.limits are not given, errors must be given')
    }
    if (is.null(lower.limits)) {
        lower.limits <- eigenvalues - errors
    }
    if (is.null(upper.limits)) {
        upper.limits <- eigenvalues + errors
    }
    taus <- list(eigenvalues=eigenvalues, lower.limits=lower.limits, upper.limits=upper.limits)
    class(taus) <- c('EigenvaluesCI', class(taus))
    return(taus)
}

#' The eigenvalues confidence intervals from an AMS analisys
#'
#' Get the eigenvalues confidence intervals resulting from the AMS analysis.
#'
#' @param analysis The AMS analysis object
#' @export
eigenvalues <- function(analysis) {
    UseMethod('eigenvalues')
}

#' @rdname eigenvalues
#' @export
eigenvalues.AMSanalysis <- function(analysis) {
    return(analysis$eigenvalues)
}
    
EigenvectorsCR <- function(conf.ellipse1, conf.ellipse2, conf.ellipse3) {
    # Eigenvectors confidence regions
    if (! (inherits(conf.ellipse1, 'SpherEllipse') && inherits(conf.ellipse2, 'SpherEllipse') && inherits(conf.ellipse3, 'SpherEllipse'))) {
        stop('arguments must be of class SpherEllipse')
    }
    vec_ellipse <- list(ellip1=conf.ellipse1, ellip2=conf.ellipse2, ellip3=conf.ellipse3)
    class(vec_ellipse) <- c('EigenvectorsCR', class(vec_ellipse))
    return(vec_ellipse)
}

#' The eigenvectors confidence regions from an AMS analisys
#'
#' Get the eigenvectors confidence regions, typically spherical ellipses,
#' resulting from the AMS analysis.
#'
#' @param analysis The AMS analysis object
#' @export
eigenvectors <- function(analysis) {
    UseMethod('eigenvectors')
}

#' @rdname eigenvectors
#' @export
eigenvectors.AMSanalysis <- function(analysis) {
    if (!inherits(analysis, 'AMSanalysis')) {
        stop('not AMSanalysis object')
    }
    return(analysis$eigenvectors)
}
    
kRejectOblate <- 'reject tau1=tau2'
kRejectProlate <- 'reject tau2=tau3'
kRejectIsotropy <- 'reject tau1=tau2=tau3'

AnisotropyTest <- function(higherLTmedium, mediumLTlower, higherLTlower) {
    somedistinct <- higherLTmedium || mediumLTlower || higherLTlower
    test <- c(higherLTmedium, mediumLTlower, somedistinct)
    class(test) <- c('AnisotropyTest', class(test))
    names(test) <- c(kRejectOblate, kRejectProlate, kRejectIsotropy)
    return(test)
}

#' The results of the anisotropy test from an AMS analisys
#'
#' The anisotropy test specifies if the hipotesis of isotropy,
#' oblate or prolate ellipsoid can be rejected, that is, if
#' the different eigenvalues can be equal.
#'
#' @param analysis The AMS analysis object
#' @export
anisotropytest <- function(analysis) {
    UseMethod('anisotropytest')
}

#' @rdname anisotropytest
#' @export
anisotropytest.AMSanalysis <- function(analysis) {
    if (!inherits(analysis, 'AMSanalysis')) {
        stop('not AMSanalysis object')
    }
    return(analysis$anisotropy_test)
}

#' @export
print.AnisotropyTest <- function(x, ...) {
    print.default(unclass(x), ...)
    if (!x[kRejectIsotropy]) {
        print('Isotropy hypothesis not rejected', ...)
    } else {
        print('Anisotropy', ...)
        if ((x[kRejectOblate] && x[kRejectProlate])
            || (!x[kRejectOblate] && !x[kRejectProlate])) {
            print('Shape not well defined', ...)
        } else if (x[kRejectOblate]) {
            print('Shape: prolate', ...)
        } else if (x[kRejectProlate]) {
            print('Shape: oblate', ...)
        } else {
            stop('assertion error')
        }
    }
}

# TODO export and document eigenParams2ellipsoidParams and ellipsoidParams2eigenParams

# tau1 >= tau2 >= tau3
eigenParams2ellipsoidParams <- function(tau1, tau2, tau3) {
    mean <- (tau1 + tau2 + tau3) / 3
    P <- tau1 / tau3
    U <- (2 * tau2 - tau1 - tau3) / (tau1 - tau3)
    result <- c(mean, P, U)
    names(result) <- c('mean', 'P', 'U')
    return(result)
}

ellipsoidParams2eigenParams <- function(mean, P, U) {
    tau1 <- 6 * P * mean / (U * (P-1) + 3 * (P+1))
    tau3 <- 6 * mean / (U * (P-1) + 3 * (P+1))
    tau2 <- 3 * mean * (1 - 2 * (P+1) / (U * (P-1) + 3 * (P+1)))
    result <- c(tau1, tau2, tau3)
    names(result) <- c('eigevvalue1', 'eigenvalue2', 'eigenvalue3')
    return(result)
}
