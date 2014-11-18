#source('utils.R')

.__vec_Aux<-c(.5, .5, 0, -1, 0, 0,
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

.__vec_Aux6x6 <- c(1, 0, 0, 0, 0, 0,
                   0, 1, 0, 0, 0, 0,
                   0, 0, 1, 0, 0, 0,
                   .5, .5, 0, 1, 0, 0,
                   0, .5, .5, 0, 1, 0,
                   .5, 0, .5, 0, 0, 1
                   )

.__mat_D_6x6 <-matrix(.__vec_Aux6x6, 6, 6, byrow=T)

.__mat_D_default <-matrix(.__vec_Aux, 15, 6, byrow=T)

#' Constructs a AMS experiment setup
#'
#' @param setup.matrix Numeric matrix.
#' @export
#' @examples
#' setup <- AMSsetup()
#' class(setup)
AMSsetup <- function(setup.matrix=.__mat_D_default) {
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

design_matrix <- function(setup) {
    if (! inherits(setup, 'AMSsetup')) {
        stop("Illegal argument: setup must be of class 'AMSsetup'")
    }
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

AMSmeasures.single <- function(positions, values) {
    reps <- rep(as.integer(1), length(values))
    AMSmeasures(reps, positions, values)
}

reps.AMSmeasures <- function(measures) {
    if (! inherits(measures, 'AMSmeasures')) {
        stop("Illegal argument: measures must be of class 'AMSmeasures'")
    }
    return(as.integer(as.character(measures$N)))
}

poss.AMSmeasures <- function(measures) {
    if (! inherits(measures, 'AMSmeasures')) {
        stop("Illegal argument: measures must be of class 'AMSmeasures'")
    }
    return(as.integer(as.character(measures$Specimen)))
}

values.AMSmeasures <- function(measures) {
    if (! inherits(measures, 'AMSmeasures')) {
        stop("Illegal argument: measures must be of class 'AMSmeasures'")
    }
    return(measures$BSus)
}

# AMSmeasures.read takes a data file and return an AMSmeasures object
AMSmeasures.read <- function(datapath) {
    frame <- read.table(datapath, header=T)
    AMSmeasures(repetitions=frame$N, positions=frame$Specimen, values=frame$BSus)
}

# AMSmeasures.mean takes an AMSmeasures object and return another AMSmeasures object
# with a single repetition and the mean value for each measured position
AMSmeasures.mean <- function(measures){
    mean_aux <- function(data, l) {
        mean.data <- mean(data)
        return(mean.data)
    } 
    means <- yysaapply(measures$BSus, measures$Specimen, mean_aux)

    poss <- as.integer(as.character(means$factors))

    AMSmeasures.single(positions=poss, values=means$values)
}


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

    D_setup <- design_matrix(setup)
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


AMSmeasures.exact <- function(sus_tensor, setup, positions = NULL) {
    D_setup <- design_matrix(setup)
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
    vec <- symtensor2vector(sus_tensor)

    values <- D %*% vec
    return(AMSmeasures(values=values, positions=positions, repetitions=as.integer(reps)))
}

AMSanalysis <- function(eigenvalues, eigenvectors, anisotropy_test) {
    if (!inherits(eigenvalues, 'AMSanalysis.eigenvalues')) {
        stop('eigenvalues must be of class AMSanalysis.eigenvalues')
    }
    if (!inherits(eigenvectors, 'AMSanalysis.eigenvectors')) {
        stop('eigenvectors must be of class AMSanalysis.eigenvectors')
    }
    if (!inherits(anisotropy_test, 'AMSanalysis.anisotropy_test')) {
        stop('anisoropy_test must be of class ams.analyisis.anisotropy_test')
    }
    results <- list(eigenvalues=eigenvalues, eigenvectors=eigenvectors, anisotropy_test=anisotropy_test)
    class(results) <- c('AMSanalysis', class(results))
    return(results)
}

AMSanalysis.eigenvalues <- function(taus_means, taus_errors, taus_low=NULL, taus_high=NULL) {
    if (is.null(taus_errors) && (is.null(taus_low) || is.null(taus_high))) {
        stop('if taus_low or taus_high are not given, taus_errors must be given')
    }
    if (is.null(taus_low)) {
        taus_low <- taus_means - taus_errors
    }
    if (is.null(taus_high)) {
        taus_high <- taus_means + taus_errors
    }
    taus <- list(taus_means=taus_means, taus_low=taus_low, taus_high=taus_high)
    class(taus) <- c('AMSanalysis.eigenvalues', class(taus))
    return(taus)
}

eigenvalues.AMSanalysis <- function(analysis) {
    if (!inherits(analysis, 'AMSanalysis')) {
        stop('not AMSanalysis object')
    }
    return(analysis$eigenvalues)
}
    
AMSanalysis.eigenvectors <- function(conf_ellipse1, conf_ellipse2, conf_ellipse3) {
    if (! (inherits(conf_ellipse1, 'spherical_ellipse') && inherits(conf_ellipse2, 'spherical_ellipse') && inherits(conf_ellipse3, 'spherical_ellipse'))) {
        stop('arguments must be of class spherical_ellipse')
    }
    vec_ellipse <- list(ellip1=conf_ellipse1, ellip2=conf_ellipse2, ellip3=conf_ellipse3)
    class(vec_ellipse) <- c('AMSanalysis.eigenvectors', class(vec_ellipse))
    return(vec_ellipse)
}

eigenvectors.AMSanalysis <- function(analysis) {
    if (!inherits(analysis, 'AMSanalysis')) {
        stop('not AMSanalysis object')
    }
    return(analysis$eigenvectors)
}
    
.__reject12 <- 'reject tau1=tau2'
.__reject23 <- 'reject tau2=tau3'
.__reject123 <- 'reject tau1=tau2=tau3'

AMSanalysis.anisotropy_test <- function(higherLTmedium, mediumLTlower, higherLTlower) {
    somedistinct <- higherLTmedium || mediumLTlower || higherLTlower
    test <- c(higherLTmedium, mediumLTlower, somedistinct)
    class(test) <- c('AMSanalysis.anisotropy_test', class(test))
    names(test) <- c(.__reject12, .__reject23, .__reject123)
    return(test)
}

anisotropy_test.AMSanalysis <- function(analysis) {
    if (!inherits(analysis, 'AMSanalysis')) {
        stop('not AMSanalysis object')
    }
    return(analysis$anisotropy_test)
}

reject1Eq2 <- function(anisotropy_test) {
    return(anisotropy_test[.__reject12])
}

reject2Eq3 <- function(anisotropy_test) {
    return(anisotropy_test[.__reject23])
}

reject1Eq2Eq3 <- function(anisotropy_test) {
    return(anisotropy_test[.__reject123])
}

print.AMSanalysis.anisotropy_test <- function(test) {
    print.default(unclass(test))
    if (!test[.__reject123]) {
        print('Isotropy hypothesis not rejected')
    } else {
        print('Anisotropy')
        if ((test[.__reject12] && test[.__reject23])
            || (!test[.__reject12] && !test[.__reject23])) {
            print('Shape not well defined')
        } else if (test[.__reject12]) {
            print('Shape: prolate')
        } else if (test[.__reject23]) {
            print('Shape: oblate')
        } else {
            stop('assertion error')
        }
    }
}

eigenvectors2ellipsoid_parameters <- function(tau1, tau2, tau3) {
    mean <- (tau1 + tau2 + tau3) / 3
    P <- tau1 / tau3
    U <- (2 * tau2 - tau1 - tau3) / (tau1 - tau3)
    result <- c(mean, P, U)
    names(result) <- c('mean', 'P', 'U')
    return(result)
}

ellipsoid_parameters2eigenvectors <- function(mean, P, U) {
    tau1 <- 6 * P * mean / (U * (P-1) + 3 * (P+1))
    tau3 <- 6 * mean / (U * (P-1) + 3 * (P+1))
    tau2 <- 3 * mean * (1 - 2 * (P+1) / (U * (P-1) + 3 * (P+1)))
    result <- c(tau1, tau2, tau3)
    names(result) <- c('eigevvalue1', 'eigenvalue2', 'eigenvalue3')
    return(result)
}
