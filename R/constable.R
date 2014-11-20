#source('utils.R')
#source('ams.R')
#source('bootstrap.R')
#source('coordinates.R')
#source('ellipses.R')

#library('boot')

#' Analyse AMS data using Bootstrap method proposed by Hext.
#'
#' This function takes a \code{AMSmeasures} object and returns the estimated 
#' eigenvalues and their errors and the estimated eigenvectors with 
#' their errors.
#' 
#' @param measures Object with class \code{AMSmeasures}.
#' @param setup Object with class \code{amssetup}.
#' @param alpha Number between 0 and 1. Confidence level.
#' @param R Integer. Number of resamplings
#' @param normalize Boolean. Should the estimations of the AMS parameters
#'     for the resampled data be normalized before mixing them?
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
#' ams.constable(measures, setup, alpha = 0.99)
#'
#' # Estimate AMS parameters and confidence intervals with 0.95 confidence 
#' # level (default) and only 100 bootstrap resamples
#' ams.constable(measures, setup, R = 100)
#'
ams.constable <- function(measures, setup, alpha=0.95, R=1000, normalize=FALSE){

    # cons_eig_para, is a data.frame with the eigenvalues and eigenvectors from the measurements

    cons_eig_para <- .__constable_eigen_param(measures, setup, R=R, normalize=normalize)

    # n, determines the measurements numbers
    
    #n <- nrow(measures) / max(poss.AMSmeasures(measures))

    # tausWithError, is a list with the three means taus and their errors

    tausError <- .__error_taus(cons_eig_para$t[,c(1,2,3)], alpha) 

    #lower.limits <- tausError$taus - tausError$error
#
    #upper.limits <- tausError$taus + tausError$error
#
    stats <- .__stat_boot(tausError$lower.limits, tausError$upper.limits)

    vec1 <- t(cons_eig_para$t[,c(4,5,6)])
    vec2 <- t(cons_eig_para$t[,c(7,8,9)])
    vec3 <- t(cons_eig_para$t[,c(10,11,12)])

    elip_vec1 <- kent_parameters(vec1, cons_eig_para$t0[c(4,5,6)], alpha=alpha)
    
    elip_vec2 <- kent_parameters(vec2, cons_eig_para$t0[c(7,8,9)], alpha=alpha)

    elip_vec3 <- kent_parameters(vec3, cons_eig_para$t0[c(10,11,12)], alpha=alpha)

    # ellipses is a list with all ellipse's parameters for the three means vectors
   #ellipses <- list(elip_vec1=elip_vec1, elip_vec2=elip_vec2, elip_vec3=elip_vec3) 
   ellipses <- EigenvectorsCR(elip_vec1, elip_vec2, elip_vec3)
    #print(tausError)

   return(AMSanalysis(tausError, ellipses, stats))
    #return (list(eigenvalues=tausError$eigenvalues,
                 #taus_error=tausError$error,
                 #lower.limits=tausError$lower.limits,upper.limits=tausError$upper.limits,ellipses=ellipses,stats=stats))

}
 

#' @importFrom boot boot
.__constable_eigen_param <- function(measures, setup, R, normalize) {

    # susc_nvec takes the susceptibility matrix and returns the susceptibility vector
    susc_nvec <- function(meas,l) {
        poss <- as.integer(meas[,1])
        values <- meas[,2]

        meas_single <- SingleMeasure(positions=poss, values=values)
        tensor <- SuscTensor(meas_single, setup)
        if (normalize) {
            tensor <- tensor/sum(diag(tensor))
        }
        return(as.vector(tensor))
    }

    # suscs return N susceptibility vectors from the N measurements
    positions <- poss.AMSmeasures(measures)
    values <- values.AMSmeasures(measures)
    repetitions <- reps(measures)
    suscs <- yysaapply(cbind(positions, values), as.factor(repetitions), susc_nvec)

    # eigen_param_v calculates the eigenparameters for R mean_susc.
    #if (normalize) {
        #eigen_param_v <- boot(suscs[, -1], lambs_normalized_nsuscs_matrix, R)
    #} else {
        eigen_param_v <- boot::boot(suscs[, -1], lambs_nsuscs_matrix, R, normalize=normalize)
    #}
    return(eigen_param_v)
}


#.__error_taus takes a matrix with 3 columns and R rows with all eigenvalues from boots and returns a list with the three eigenvalues and their errors. Desv_st calculates the sd from the r rows for the three eigenvalues, and error_mean_value is the total error, where qt is the t'student for alpha=95%

.__error_taus <- function(r_eig_values, alpha){

        tau_means <- colMeans(r_eig_values)

        desv_st <- apply(r_eig_values,2,sd)
        alpha_margin <- 1 - alpha

        #print(r_eig_values)
        #error_mean_value <- qt(p=(1-alpha_margin/2), df=(n-1)) * desv_st / sqrt(n)
        error_mean_value <- qnorm(p=1-alpha_margin/2)* desv_st 
        #error_mean_value <- desv_st

        return(EigenvaluesCI(eigenvalues=tau_means, errors=error_mean_value))
        #return(list(taus=tau_means, errors= error_mean_value))

}

#fisher_parameters <- function(nvectors) {
    #
    #sum_vector <- colSums(nvectors)
    #R <- sqrt(sum_vector %*% sum_vector)
    #mean_vector <- sum_vector/R
    #n <- nrow(nvectors)
    #alpha95 <- acos(1 - (n - R) * (20^(1/(n - 1)-1)/R))
    #alpha95_prima <- 140/sqrt((n - 1)/(n - R))
#
    #return(mean_vector,alpha95,R)
#}


.__stat_boot <- function(low_error, high_error){
    low_error_tau1 <- low_error[1]
    low_error_tau2 <- low_error[2]
    low_error_tau3 <- low_error[3]
    high_error_tau1 <- high_error[1]
    high_error_tau2 <- high_error[2]
    high_error_tau3 <- high_error[3]

    # TODO: Too naive. It does not provide a generally righ result for the confidence level. Actually, the confidence level is not even considered.
    # TODO: provide a general statistical method, under some simple asumptions (like normaliy), more similar to what is done fr hext.
    test <- AnisotropyTest(high_error_tau2 < low_error_tau1, high_error_tau3 < low_error_tau2, high_error_tau3 < low_error_tau1)

    return(test)
    
}



# kent_parameters takes a matrix (nvectors) with n columns and 3 rows, where n are the values for one eigenvector. And it returns all parameters from confidence ellipse.

kent_parameters <- function(nvectors, ref_vector, n=1, alpha) {
    nvectors <- .__closest_vect(nvectors, ref_vector)

    mean_vec <- rowMeans(nvectors)
    mean_mat <- nvectors %*% t(nvectors)/ncol(nvectors)
    coor <- car2sph(t(mean_vec))
    lat <- coor[2]
    phi <- coor[1]
    theta <- pi/2 - lat
    H_vec <- c(cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta), -sin(phi), cos(phi), 0, sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
    H <- matrix(3,3, data=H_vec)
    B <- t(H) %*% mean_mat %*% H
   # Bl <- B[-3,-3]
   # l <- eigen(Bl, only.values=T)

    if (B[1,1]==B[2,2]){
        pshi<- pi/2
    } else {
        pshi <- atan(2*B[1,2]/(B[1,1]-B[2,2]))/2
    }
    k_vector <- c(cos(pshi), sin(pshi), 0, -sin(pshi), cos(pshi), 0, 0, 0, 1)
    k <- matrix(k_vector, 3,3)
    gamma_mat <- H %*% k

    r1 <- mean_vec %*% mean_vec
   # r2 <- l$values[1]-l$values[2]
    
    mean_vec <- mean_vec/sqrt(r1)

    xprima <- t(gamma_mat) %*% nvectors
    mu <- mean(xprima[3,])
    sigma1 <- mean(xprima[1,]^2)
    sigma2 <- mean(xprima[2,]^2)
    alpha_margin <- 1 - alpha

    # TODO comprobar si es n del numero de medidas iniciales o r numero de boostraps en la g

    g = qchisq(df=2, p=1-alpha_margin)/(n * mu^2)

    if (sigma1 * g < 1){
        eta = asin(sqrt(sigma1 * g))
    }else{
        eta = pi/2
    }

    if (sigma2 * g < 1){
        zeta = asin(sqrt(sigma2 * g))
    }else{
            zeta = pi/2
    }

    if(! orthogonality_check(gamma_mat)) {
        print('not orthogonal matrix')
        print(gamma_mat)
    }
    edir = car2sph(gamma_mat[1,1], gamma_mat[2,1], gamma_mat[3,1])
    zdir = car2sph(gamma_mat[1,2], gamma_mat[2,2], gamma_mat[3,2])

    #return(list(eta=eta, zeta=zeta, edir=edir, zdir=zdir, center=car2sph(t(mean_vec))))
    return(SpherEllipse(centerDir=car2sph(t(mean_vec)), axis1Dir=edir, axis2Dir=zdir, semiangle1=eta, semiangle2=zeta))
}

.__closest_vect <- function(vects, ref) {
    invert <- t(ref) %*% vects < 0
    auxf <- function(vs) {
        vs[invert] <- - vs[invert]
        return(vs)
    }
    return(t(apply(vects, MARGIN=1, FUN=auxf)))
}

