#source('ams.R')
#source('coordinates.R')
#source('ellipses.R')

#' Analyse AMS data using Linear Perturbation Analysis method proposed by Hext
#'
#' @param measures Object with class \code{AMSmeasures}.
#' @param setup Object with class \code{amssetup}.
#' @param alpha Number between 0 and 1. Confidence level.
#' @param withoutN Boolean indicating if the old version,
#' extracted from some bibliography and not taking into
#' account properly the number of measures,
#' should be used. The default FALSE is recommended.
#' Setting it to TRUE gives unexpected results, with
#' errors much bigger than expected, and not scaling
#' properly with the number of measures. This parameter
#' is only kept with the purpose of being able to
#' reproduce old results.
#'
#' @export
#' @examples
#' # Sample AMS measures object
#' measures <- as.AMSmeasures(sample_measures)
#'
#' # Experimental setup object
#' setup <- AMSsetup()
#' 
#' # Estimate AMS parameters and confidence intervals with 0.99 confidence 
#' # level
#' ams.hext(measures, setup, alpha = 0.99)
#'
ams.hext <- function(measures, setup, alpha=0.95, withoutN=F){

    param <- .__hext_param(measures,setup, withoutN)

    N <- nrow(measures) / max(poss.AMSmeasures(measures))

    tausError <- .__errors_taus(param$eig_param,param$varianza,param$V,N, alpha, withoutN)

    #print(tausError$taus)

    #lower.limits <- tausError$taus-tausError$error
#
    #upper.limits <- tausError$taus+tausError$error
#
    ellipses <- .__errors_vec(param$eig_param, param$varianza, param$xb,nrow(measures),alpha)

    #data <- tausError$taus
    #
    stats <- .__statHext(tausError$eigenvalues, param$eig_param, param$V, param$varianza, param$xb, nrow(measures), alpha, withoutN)

    return(AMSanalysis(tausError, ellipses, stats))
   #return(list(eigenvalues=tausError$eigenvalues,
               #taus_error=tausError$error,
               #lower.limits=tausError$lower.limits, upper.limits=tausError$upper.limits,ellipses=ellipses,stats=stats))
}

# .__hext_param takes a frame with the n measurements and returns a list with a data.frame with the eigenvalues and the eigenvectors, the variance and the trace of susceptibility tensor

.__hext_param <- function(measures,setup, withoutN=F) {


    values <- values.AMSmeasures(measures)

    positions <- poss.AMSmeasures(measures)

    #suscHext, susceptibility matrix
    suscHext <- SuscTensor(measures,setup)

    measur_fit <- ExactMeasures(suscHext,setup,positions)

    # resid, is a vector with the residuals

    resid <- values - measur_fit[,3]

    resid.tot <- sum(resid^2)

    # vari is the variance

    N<- nrow(measures)

    # Degrees of freedom, taken as the total number of measures
    # minus the 6 independent components of the tensor (see Hext1963)
    free_g <- N-6
    if (withoutN) {## old
        vari <- sqrt(resid.tot/free_g)
    }
    else {##yy new
        repetitions <- N/max(positions) 
        vari <- sqrt(resid.tot/(free_g*repetitions))
    }

    #print(vari)
    #print(paste('sigma',vari))
    # eig_param, is a data.frame with the eigenvalues and the eigenvectors

    eig_param <- eigen(suscHext)
    if(! orthogonality_check(eig_param$vectors)) {
        print('not orthogonal matrix')
        print(eig_param$vectors)
    }

    # xb is the susceptibility trace

    xb <- (suscHext[1,1] + suscHext[2,2] + suscHext[3,3])/3

    D_setup <- design(setup)
    V <- solve(t(D_setup) %*% D_setup)

    return (list(eig_param=eig_param,varianza=vari,xb=xb, V=V))
}

# .__errors_taus, takes the eigen-parameters and the variance and returns a list with the eigenvalues and their errors

.__errors_taus <- function(eig_param,vari,positions_V,N, alpha, withoutN=F){
    
    #D_setup <- design(setup)

    #print(D_setup)
    #positions_V <- solve(t(D_setup) %*% D_setup)
    #print('V')
    #print(positions_V)

    alpha_margin <- 1 - alpha

    if (withoutN) {##old
        ts <- qnorm(1-alpha_margin/2)/sqrt(N)
    }
    else {##new
        ts <- qnorm(1-alpha_margin/2)
    }

    vectors <- eig_param$vectors

    a11 <- a_tensor(vectors[,1], vectors[,1])
    error_tau1 <- vari * ts * sqrt( t(a11) %*% positions_V %*% a11)

    a22 <- a_tensor(vectors[,2], vectors[,2])
    error_tau2 <- vari * ts * sqrt( t(a22) %*% positions_V %*% a22)

    a33 <- a_tensor(vectors[,3], vectors[,3])
    error_tau3 <- vari * ts * sqrt( t(a33) %*% positions_V %*% a33)

    error_taus <- c(error_tau1, error_tau2, error_tau3)

    return(EigenvaluesCI(eigenvalues=eig_param$values, errors=error_taus))
    #return(list(taus=eig_param$values, error=error_taus))

   }


# .__errors_vec takes the data.frame with the eigenvalues and the eigenvectors, variance and the trace and returns all ellipses parameters and the eigenvector in polars coordinates

.__errors_vec <- function(param,vari,xb,N,alpha){

    # eliij are the semiangles of ellipses and edir and zdir the directions

    f<- sqrt(2*qf(p=alpha, df1=2,df2=N-6))

    eli12 <-atan(f*vari/(2*(param$values[1]-param$values[2])))

    eli23 <-atan(f*vari/(2*(param$values[2]-param$values[3])))

    eli13 <-atan(f*vari/(2*(param$values[1]-param$values[3])))
    
    coor_polar <- car2sph(t(param$vectors))

    ellip1 <- NorthHemisphEllipse(centerDir=coor_polar[1,], axis1Dir=coor_polar[2,], axis2Dir=coor_polar[3,], semiangle1=eli12, semiangle2=eli13)
    ellip2 <- NorthHemisphEllipse(centerDir=coor_polar[2,], axis1Dir=coor_polar[1,], axis2Dir=coor_polar[3,], semiangle1=eli12, semiangle2=eli23)
    ellip3 <- NorthHemisphEllipse(centerDir=coor_polar[3,], axis1Dir=coor_polar[1,], axis2Dir=coor_polar[2,], semiangle1=eli13, semiangle2=eli23)
    #elip_vec1 <- list(eta=eli12, zeta=eli13, edir=coor_polar[2,],zdir=coor_polar[3,], center=coor_polar[1,])
#
    #elip_vec2 <- list(eta=eli12, zeta=eli23, edir=coor_polar[1,], zdir=coor_polar[3,],center=coor_polar[2,])
#
    #elip_vec3 <- list(eta=eli13, zeta=eli23, edir=coor_polar[1,], zdir=coor_polar[2,],center=coor_polar[3,])
#
    #return (list(elip_vec1=elip_vec1, elip_vec2=elip_vec2, elip_vec3=elip_vec3))
    return(EigenvectorsCR(ellip1, ellip2, ellip3))
}

# .__statHext takes the eigenvalues, the variance and the trace determines if the anisotropy of magnetic susceptibility (prolate, oblate)

    .__statHext <- function(taus, eig_params, V, vari, xb, N, alpha, withoutN=F) {
    tau1 = taus[[1]]
    tau2 = taus[[2]]
    tau3 = taus[[3]]

    vec1 <- eig_params$vectors[,1]
    vec2 <- eig_params$vectors[,2]
    vec3 <- eig_params$vectors[,3]

    # TODO: Differencies between this two versions? Which one is right?
    if (withoutN) {##old
        Fh <- 0.4 * (tau1^2 + tau2^2 + tau3^2 - 3*xb^2)/(vari^2)
        F12 <- 0.5 * ((tau1 - tau2) /vari)^2
        F23 <- 0.5 * ((tau2 - tau3) /vari)^2
        qF2 <- qf(alpha, df1=2, df2=N-6)    #4.2565
        qF5 <- qf(alpha, df1=5, df2=N-6)    #3.4817
        rF12 <- F12 > qF2
        rF23 <- F23 > qF2
        rFh <- Fh > qF5

        test <- AnisotropyTest(rF12, rF23, rFh)
    } else {##new
        a1 <- a_tensor(vec1, vec1)
        a2 <- a_tensor(vec2, vec2)
        a3 <- a_tensor(vec3, vec3)
        vec12 <- a1 - a2
        vec23 <- a2 - a3
        vec13 <- a1 - a3
        F12 <- (tau1 - tau2) / (vari * sqrt(t(vec12) %*% V %*% vec12))
        F23 <- (tau2 - tau3) / (vari * sqrt(t(vec23) %*% V %*% vec23))
        F13 <- (tau1 - tau3) / (vari * sqrt(t(vec13) %*% V %*% vec13))

        qTN <- qt(p=alpha, df=N-6)
        rF12 <- F12 > qTN
        rF23 <- F23 > qTN
        rF13 <- F13 > qTN
        test <- AnisotropyTest(rF12, rF23, rF13)
    }
    return(test)

    #if (F12 > 4.2565) {
        #F12o <- 'different lambdas'
    #} else {
        #F12o <- 'oblate'
    #}
#
    #if (F23 > 4.2565) {
        #F23p <- 'different lambdas'
    #} else {
        #F23p <- 'prolate'
    #}
#
    #if (Fh > 3.4817) {
        #Fs <- 'anisotropy'
    #} else {
        #Fs <- 'sphere'
    #}
    #return (cbind(Fs,F12o,F23p))
}
