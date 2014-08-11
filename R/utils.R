
# yysaapply takes a function (FUN), one frame's column (x) and one frame's factor column (by) and returns a dataframe which first column are the factor's levels and the results from the function FUN

yysaapply <- function(x, by, FUN){
    factors <- levels(by)
    values <- c()
    i <- 0
    for (l in factors){
        i <- i + 1
        if (is.vector(x)){
            values <- c (values, FUN(x[by ==l],l))
        }else{
            values <- c(values, FUN(x[by == l, ], l))
        }
    }
    ncol <- length(values) / length(factors)
    if (ncol > 1) {
        mat <- matrix(ncol = ncol, nrow = length(factors), byrow = T, data = values)
        values <- data.frame(mat)
    }
    return(data.frame(factors, values))
}

suscep_matrix <- function(lamb_vector, rot_vector=c(0,0,0)) {

    # Eigenvalues vector

        lamb <- c(lamb_vector[1], 0 ,0,
                  0, lamb_vector[2], 0,
                  0, 0, lamb_vector[3])

        lamb_matrix <- matrix(3,3, byrow = T, data = lamb)

    # Angles for rotation

        theta <- rot_vector[1]
        phi <- rot_vector[2]
        psi <- rot_vector[3]


        R <- c(cos(theta)*cos(phi)*cos(psi)-sin(phi)*sin(psi),
             cos(theta)*sin(phi)*cos(psi)+cos(phi)*sin(psi),
             -sin(theta)*cos(psi), 
          -cos(theta)*cos(phi)*sin(psi)-sin(phi)*cos(psi),
             -cos(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi),
             sin(theta)*sin(psi),
          sin(theta)*cos(phi),
             sin(theta)*sin(phi),
             cos(theta))

        rot_matrix <- matrix(3,3, byrow = T, data = R)

    # Susceptibility matrix

         sus_matrix <- t(rot_matrix) %*% lamb_matrix %*% rot_matrix

    return(sus_matrix)
}

# Transform a 3x3 symmetric tensor to a 6-vector representation.
# Inverse of vector2symtensor
symtensor2vector <- function(tensor) {
    return(c(tensor[1,1], tensor[2,2], tensor[3,3], tensor[1,2], tensor[2,3], tensor[1,3]))
}

# Transform a 6-vector to a 3x3 symmetric matrix representation.
# Inverse of symtensor2vector
vector2symtensor <- function(vec) {
    vec9 <- c(vec[1], vec[4], vec[6], 
              vec[4], vec[2], vec[5],
              vec[6], vec[5], vec[3])
    return(matrix(3, 3, byrow=T, data=vec9))
}

# a_tensor allows products like
# t(dir_i) %*% S %*% dir_j
# where S is a simmetric matrix, to be written as
# a_tensor(dir_i, dir_j) %*% s
# where s = (S11, S22, S33, S12, S23, S13) = symtensor2vector(S)
a_tensor <- function(dir_i, dir_j) {
    return( c(dir_i[1] * dir_j[1], 
              dir_i[2] * dir_j[2], 
              dir_i[3] * dir_j[3],
              dir_i[1] * dir_j[2] + dir_i[2] * dir_j[1],
              dir_i[2] * dir_j[3] + dir_i[3] * dir_j[2],
              dir_i[3] * dir_j[1] + dir_i[1] * dir_j[3]) )
}


# Returns the matrix representing an active rotation of axis Z and the given
# angle
rotate_Z <- function(angle) {
    rot_mat <- cbind(c(cos(angle),sin(angle),0), c(-sin(angle),cos(angle),0), c(0,0,1))
    return(rot_mat)
}

# Returns the matrix representing an active rotation of axis X and the given
# angle
rotate_X <- function(angle) {
    rot_mat <- cbind(c(1,0,0), c(0, cos(angle), sin(angle)), c(0, -sin(angle),cos(angle)))
    return(rot_mat)
}

# Returns the matrix representing an active rotation of axis Y and the given
# angle
rotate_Y <- function(angle) {
    rot_mat <- cbind(c(cos(angle), 0,-sin(angle)), c(0,1,0), c(sin(angle),0,cos(angle)))
    return(rot_mat)
}

orthogonality_check <- function(mat, epsilon=0.001) {
    if (norm(t(mat) %*% mat - diag(3)) > epsilon) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

