data_cons <- function(d1,d2,d3,d4,d5,d6){

    k11 <- read.table(d1, sep=",")
    k22 <- read.table(d2, sep=",")
    k33 <- read.table(d3, sep=",")
    k12 <- read.table(d4, sep=",")
    k13 <- read.table(d5, sep=",")
    k23 <- read.table(d6, sep=",")

#    k12 <- k12s[,2] - .5 * (k11[,2] + k22[,2])
 #   k13 <- k13s[,2] - .5 * (k11[,2] + k33[,2])
  #  k23 <- k23s[,2] - .5 * (k22[,2] + k33[,2])

    xn <- seq(from=0, to=1, len=22)

    xl <-xn[-c(1,22)]

    x <- qnorm(xl)

    xerr <- cbind(k11[,1]-x,k22[,1]-x,k33[,1]-x,k12[,1]-x,k13[,1]-x,k23[,1]-x)
    # print(xerr)

    N <- rep(c(1:nrow(k11)),6)
    Specimen <- rep(c(1:6), each=nrow(k11))
    BSus <- c(k11[,2],k22[,2],k33[,2],k12[,2],k13[,2],k23[,2]) 

    #measures <- data.frame(N=N, Specimen=Specimen, BSus=BSus)
    measures <- ams.measures(repetitions=N, positions=Specimen, values=BSus)

    return(measures)
}

data_vec_cons <- list('boot'=cbind('taus'=c(3.790,3.700,3.460),'dec'=c(278.8,6.3,26.4),'inc'=c(6.7,-20.2,68.6),'eta'=c(16.2,16.2,7.9),'zeta'=c(6.1,6.4,4.4)),'boot_norm'=cbind('taus'=c(0.346,0.338,0.316),'dec'=c(-81.2,-173.7,26.4),'inc'=c(6.7,20.2,68.6),'eta'=c(6.1,6.4,4.4),'zeta'=c(16.2,16.2,7.9),'edir_dec'=c(-145.7,-145.6,-41.5),'edir_inc'=c(-74.6,-67.3,-8.3),'zdir_dec'=c(6.8,-80.4,-128.6),'zdir_inc'=c(-13.8,10.0,19.5)),'hext'=cbind('taus'=c(3.787,3.703,3.457),'dec'=c(-81.1,-173.6,26.5),'inc'=c(6.8,20.3,68.5),'eta'=c(6.1,8.9,6.1),'zeta'=c(17.2,17.2,8.9),'edir_dec'=c(-153.5,-153.5,-81.0),'edir_inc'=c(-68.5,-68.5,6.8),'zdir_dec'=c(6.4,-81.7,-173.6),'zdir_inc'=c(-20.3,6.8,20.3)))
