
# car2sph takes a matrix with the directions in (x,y,z) and returns a matrix which the first column is the longitud, the second columns is the latitud and the third, the radius

car2sph <-
function(x,y,z,deg=FALSE){
	if(is.matrix(x) || is.data.frame(x)){
		if(ncol(x) == 1){x=x[,1]}
        else	if(ncol(x) == 2){y = x[, 2];x = x[, 1]}
        else	if(ncol(x) == 3){z = x[, 3];y = x[, 2];x = x[, 1]}
    }
    if(missing(x) | missing(y) | missing(z)){stop('Missing full cartesian 3D input data.')}
	radius = sqrt(x^2 + y^2 + z^2)
    long = atan2(y, x)
    lat = asin(z/radius)
	if(deg){
		long=long*180/pi
		lat=lat*180/pi
	}
	lat[radius==0]=0
	return=cbind(long=long,lat=lat,radius=radius)
	}

# sph2car takes a matrix with the directions in (long, lat, radius) and returns a matrix with (x,y,z)

sph2car <-
function(long,lat,radius=1,deg=FALSE){
    if(is.matrix(long) || is.data.frame(long)){
        if(ncol(long) == 1){long = long[,1]}
		else	if(ncol(long) == 2){lat = long[, 2];long = long[, 1]}
		else	if(ncol(long) == 3){radius = long[, 3];lat = long[, 2];long = long[, 1]}
		}
	if(missing(long) | missing(lat)){stop('Missing full spherical 3D input data.')}
    if(deg){
    	long = long * pi/180
		lat = lat * pi/180
	}
return=cbind(x=radius*cos(long)*cos(lat),y=radius*sin(long)*cos(lat),z=radius*sin(lat))
}

# spher2lambert  takes two vectors (long,lat) and returns rho (inclination by lambert projection from the south pole) and long.

spher2lambert <- function(long,lat) {
    return(cbind(rho=.__lat2rho(lat), long))
} 

lambert2spher <- function(rho,long) {
    return(cbind(long,lat=.__rho2lat(rho)))
}

# lambert_NorthProlongation transforms directions in the south hemisphere to its prolongation in the north hemisphere

lambert_NorthProlongation <- function(rho,long) {
    cond <- rho > 1
    rho[cond] <- .__lat2rho(-.__rho2lat(rho[cond]))
    long[cond] <- (long[cond] + pi) %% (2*pi)
    return(cbind(rho, long))
}

# .__lambert_divide_hemisphere divides a set of directions in several blocks, each block
# containing consecutive points lying in the same hemisphere. Equator points are 
# asigned to the north hemisphere. The function returns a list of such blocks
.__lambert_divide_hemisphere <- function(rho,long) {
    result <- list()
    i <- 1
    blockstart <- 1
    issouth <- rho[i] > 1
    while (i <= length(rho)) {
        while (i <= length(rho) && issouth == (rho[i] > 1)) {
            i <- i + 1
        }
        result[[length(result) + 1]] <- list(rho=rho[blockstart : (i-1)], long=long[blockstart : (i-1)])
        issouth = !issouth
        blockstart <- i
    }
    return(result)
}

.__lat2rho <- function(lat) {
    return(cos(lat/2) - sin(lat/2))
}

.__rho2lat <- function(rho) {
    colat <- 2 * acos(rho/sqrt(2)) 
    return(-(pi/2) + colat)
}


#' Plot in lambert coordinates
#'
#' Plots 3-D directions given in spherical coordinates, longitud and latitud
#' using a Lambert azimuthal projection (equal area projection). Only the 
#' north hemisphere is represented, direcions in the south hemisphere and 
#' reverted and then represented in the north hemisphere.
#'
#' @importFrom plotrix polar.plot
#'
#' @param long Numeric vector of longitud coordinates in radians
#' @param lat Numeric vector of latitud coordinates in radians
#' @param radial.labels Numeric vector of labels of the longitud coordinates,
#'     in degrees, default is \code{c(90, 60, 30, 0)}
#' @param add Boolean indicating if the plot should be started from zero,
#'     removing objectx previously in the plot, default is False
#' @param rp.type Like in \code{\link{polar.plot}}, default is 's' 
#'     (symbol)
#' @param divideHemispheres Boolean value indicating if the plot should be
#'     divided into blocks of points lying in the same hemisphere, so that
#'     points in different hemisphere are not joined. Default is True
#' @param ... Extra parameters to be passed to \code{\link{polar.plot}}
#'
#' @export
lambert.plot <- function(long,lat,radial.labels=c(90,60,30,0),add=F,rp.type='s',divideHemispheres=T,...) {
    #library(plotrix)
    rad_lim <- .__lat2rho(radial.labels*pi/180)

    lam_vec <- spher2lambert(long,lat)

    if (divideHemispheres) {
        lam_vec_divided <- .__lambert_divide_hemisphere(lam_vec[,1], lam_vec[,2])
        first <- T
        for (block in lam_vec_divided) {
            lam_vec_correc <- lambert_NorthProlongation(block$rho, block$long)
            
            #to avoid the polygon behaviour and simulate lines
            rhos <- lam_vec_correc[,1]
            longs <- lam_vec_correc[,2]
            rhos <- c(rhos, rev(rhos))
            longs <- c(longs, rev(longs))

            if (first) {
                first <- F
                result <- plotrix::polar.plot(rhos,longs*180/pi,radial.labels=radial.labels, radial.lim=rad_lim, add=add, rp.type=rp.type,...)
            } else {
                plotrix::polar.plot(rhos,longs*180/pi,radial.labels=radial.labels, radial.lim=rad_lim, add=T, rp.type=rp.type,...)
            }
        }
        return(result)

    } else {
        lam_vec_correc <- lambert_NorthProlongation(lam_vec[,1],lam_vec[,2])

        return(plotrix::polar.plot(lam_vec_correc[,1],lam_vec_correc[,2]*180/pi,radial.labels=radial.labels, radial.lim=rad_lim, add=add, rp.type=rp.type,...))
    }
}

# plots 3-D directions given in cartesian coordinates using a lambert 
# projection. Only the north hemisphere is represented, direcions in the 
# south hemisphere and reverted and then represented in the north hemisphere.
lambert.plot.xyz <- function(x, y, z, ...) {
    sph <- car2sph(x, y, z)
    lon <- sph[,1]
    lat <- sph[,2]
    return(lambert.plot(lon, lat, ...))
}

