# mean_nsuscs_matrix takes a matrix (nsuscs_matrix) where each row is an one observation for the boot's function. It returns the mean for each susceptibility tensor's element.  

mean_nsuscs_matrix <- function(nsuscs_matrix, d) {
    return(colMeans(nsuscs_matrix[d,]))
}

# lambs_normalized_nsuscs_matrix takes a matrix like mean_nsuscs_matrix and returns the eigenvalues of the mean susceptibility tensor. 

#lambs_normalized_nsuscs_matrix <- function(nsuscs_matrix, d) {
    #return(.__lambs_nsuscs_matrix_aux(nsuscs_matrix, d, normalize=TRUE))
#}
# lambs_nsuscs_matrix takes a matrix like mean_nsuscs_matrix and returns the eigenvalues of the mean susceptibility tensor. 

#lambs_nsuscs_matrix <- function(nsuscs_matrix, d) {
    #mean_susc <- mean_nsuscs_matrix(nsuscs_matrix, d)
    #eigen_param <- eigen(matrix(mean_susc, nrow=sqrt(length(mean_susc))))
#
    #eigen_param_vect <- .__boot.convert(eigen_param)
#
    #return(eigen_param_vect)
    #return(.__lambs_nsuscs_matrix_aux(nsuscs_matrix, d, normalize=FALSE))
#}
# lambs_nsuscs_matrix takes a matrix like mean_nsuscs_matrix and returns the eigenvalues of the mean susceptibility tensor. 

lambs_nsuscs_matrix <- function(nsuscs_matrix, d, normalize) {
   # print(d)
    mean_susc <- mean_nsuscs_matrix(nsuscs_matrix, d)
    mean_tensor <- matrix(mean_susc, nrow=sqrt(length(mean_susc)))
    if (normalize) {
        mean_tensor <- mean_tensor/sum(diag(mean_tensor))
    }
   # print(mean_tensor)
    eigen_param <- eigen(mean_tensor)

    eigen_param_vect <- .__boot.convert(eigen_param)

    return(eigen_param_vect)
}

# .__boot.convert takes a list from the eigen's function and converts this list in a vector
.__boot.convert <- function(eigen_param) {
    
    eigen_param_vec <- c(eigen_param$values, eigen_param$vectors[,1], eigen_param$vectors[,2], eigen_param$vectors[,3])

    return(eigen_param_vec)
}

# .__boot.invert take a vector with the eigenvalues and eigenvectors and converts it in a list
.__boot.invert <- function(eigen_param_vec) {
    eigen_param_list <- list(values = c(eigen_param_vec[1], eigen_param_vec[2], eigen_param_vec[3]), vectors = matrix(3,3, data = c(eigen_param_vec[4],eigen_param_vec[5],eigen_param_vec[6], eigen_param_vec[7], eigen_param_vec[8],eigen_param_vec[9])))

    return(eigen_param_list)
}
