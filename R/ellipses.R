
#' Constructs an spherical ellipse
#'
#' An ellipse in the sphere is described by a list object 'elip' with constituent 
#' elements 'elip$center', 'elip$eta', 'elip$zeta', 'elip$edir', 'elip$zdir'. 
#' Of these, 'elip$center', 'elip$edir', 'elip$zdir', are vectors containing the
#' longitud, latitud and radial coordinates of the center direction and the two
#' semiangles, called eta and zeta. 'elip$eta' and 'elip$zeta' are the value of 
#' the two semiangles. All angular magnitudes in radians.
#'
#' @param center.dir Numeric vector, representing longitud and latitude of the
#'     center of the ellipse
#' @param axis1.dir Numeric vector, representing longitud and latitude of the
#'     first semiaxis of the ellipse
#' @param axis2.dir Numeric vector, representing longitud and latitude of the
#'     second semiaxis of the ellipse
#' @param semiangle1 Numberic value, angular magnitud of the first semiaxis
#' @param semiangle2 Numberic value, angular magnitud of the second semiaxis
#' @export
#' @examples
#' center.dir <- c(0, pi/2)
#' axis1.dir <- c(pi/2, 0)
#' axis2.dir <- c(pi/2, pi/2)
#' ellipse <- SpherEllipse(center.dir, axis1.dir, axis2.dir, 1/2, 1)
#' plot(ellipse)
#'
#' @seealso \code{\link{plot.SpherEllipse}}
SpherEllipse <- function(center.dir, axis1.dir, axis2.dir, semiangle1, semiangle2) {
    if (semiangle1 > semiangle2) {
        eta <- semiangle1
        edir <- axis1.dir
        zeta <- semiangle2
        zdir <- axis2.dir 
    } else {
        eta <- semiangle2
        edir <- axis2.dir
        zeta <- semiangle1
        zdir <- axis1.dir 
    }
    ellip <- list(center=center.dir, eta=eta, zeta=zeta, edir=edir, zdir=zdir)
    class(ellip) <- c('SpherEllipse', class(ellip))
    return(ellip)
}

NorthHemisphEllipse <- function(center.dir, axis1.dir, axis2.dir, semiangle1, semiangle2) {
    return(SpherEllipse(.__toNorth(center.dir),
                            .__toNorth(axis1.dir),
                            .__toNorth(axis2.dir),
                            semiangle1,
                            semiangle2))
}

# The ellipse has a preferred reference system in the sphere, where the center
# lies in the z-axis, eta in the x-axis and zeta in the y-axis. In such a reference
# system, the  only free parameters are 'elip$eta' and 'elip$zeta', with
# 'elip$center' = c(0, pi/2, 1), 'elip$edir' = c(0, 0, 1), 
# 'elip$zdir' = c(pi/2, 0, 1), and the equation of the ellipse simplifies to
# 'cos(lat)^2 * ( (cos(lon)/sin(eta))^2 + (sin(lon)/sin(zeta))^2 ) = 1', with
# 'lon' and 'lat' the longitud and latitud coordinates in the ellipse's reference 
# system. See Hext1963


# Returns the rotation matrix that transforms from the ellipse's preferred reference system
# to the original reference system.
.__rotate_from_ellipse_system <- function(elip) {
    pole_rotate <- with(elip, rotate_Z(center[1]) %*% rotate_Y(pi/2 - center[2]))
    north <- with(elip, c(-cos(center[1])*sin(center[2]), -sin(center[1])*sin(center[2]), cos(center[2])))
    eta_dir <- with(elip, c(cos(edir[1])*cos(edir[2]), sin(edir[1])*cos(edir[2]), sin(edir[2])))
    center <- with(elip, c(cos(center[1])*cos(center[2]), sin(center[1])*cos(center[2]), sin(center[2])))
    
    s_eta <- det( cbind(center,north,eta_dir))
    c_eta <- eta_dir %*% north

    # R handles the infinities properly, so there is no need of checking for zero values
    eta_angle <- atan(s_eta/c_eta)

    total_rotate <- pole_rotate %*% rotate_Z(eta_angle)
    return(total_rotate)
}

# Returns the rotation matrix that transforms to the ellipse's preferred
# reference system
.__rotate_to_ellipse_system <- function(elip) {
    return(t(.__rotate_from_ellipse_system(elip)))
}

# Construct 'npoints' in elipse 'elip'. Useful for plotting.
ellipse_points <- function(elip, npoints=200) {
    rlong <- seq(from=0, to=2*pi, length.out=npoints)
    rlat <- acos( sqrt( 1 / ( (cos(rlong)/sin(elip$eta))^2 + (sin(rlong)/sin(elip$zeta))^2) ) )

    rpoints <- cbind(rlong, rlat, rep(1, times=npoints))
    cart_rpoints <- sph2car(rpoints)
    cart_points <- .__rotate_from_ellipse_system(elip) %*% t(cart_rpoints)
    points <- car2sph(t(cart_points))
    return(points)
}

#' Plot spherical ellipse
#'
#' Plots an spherial ellipse in Lambert azimuthal (equal area)
#' projection. Only the north hemisphere
#' is represented, directions in the south hemisphere are reverted and
#' represented in the north pole. 
#'
#' @param x An SpherEllipse object
#' @param npoints The number of points to be plotted, default is 200
#' @param add Boolean indicating if the plot should be drawn over an existing
#'     plot (\code{True}) or on a new one (\code{False}, default)
#' @param line.col The color of the ellipse, default is 'red'
#' @param ... Other parameters to be passed to \code{\link{lambert.plot}}
#' @export
plot.SpherEllipse <- function(x, npoints=200, add=FALSE, line.col='red', ...) {
    points <- ellipse_points(x, npoints)
    lambert.plot(points[,1], points[,2], rp.type='p', line.col=line.col, add=add, ...)
}

# sph_point as returned and acepted by the conversion functions car2ph and sph2car
# in coordinates.R, with rows lon, lat, radious.
# If 'direction' is TRUE, sph_point is considered as a direction, and the 
# opposite vector is also checked.
in.ellipse <- function(elip, sph_point, direction=FALSE) {
    if (direction) {
        sph_opp <- .__opposite(sph_point)
        #sph_opp[,1] <- (sph_opp[,1] + pi) %% (2*pi)
        #sph_opp[,2] <- -sph_opp[,2]
        return(in.ellipse(elip, sph_point, FALSE) || in.ellipse(elip, sph_opp, FALSE))
    }
    
    car_point <- sph2car(sph_point)
    r_cart_point <- .__rotate_to_ellipse_system(elip) %*% t(car_point)
    r_sph_point <- car2sph(t(r_cart_point))
    
    lon <- r_sph_point[,1]
    lat <- r_sph_point[,2]

    eta <- elip$eta
    zeta <- elip$zeta
    
    distance <- cos(lat)^2 * ( (cos(lon)/sin(eta))^2 + (sin(lon)/sin(zeta))^2 ) 
    if (distance > 1) {
        #print('-----')
        #print(elip)
        #print(sph_point)
        #print(distance)
    }
    return(distance <= 1)
}

.__opposite <- function(sphcoord) {
    sph_opp <- sphcoord
    if (is.vector(sph_opp)) {
        sph_opp[1] <- (sph_opp[1] + pi) %% (2*pi)
        sph_opp[2] <- -sph_opp[2]
    } else {
        sph_opp[,1] <- (sph_opp[,1] + pi) %% (2*pi)
        sph_opp[,2] <- -sph_opp[,2]
    }
    return(sph_opp)
}

.__toNorth <- function(sphcoord) {
    if (is.vector(sphcoord)) {
        if (sphcoord[2] < 0) {
            return(.__opposite(sphcoord))
        }
        return(sphcoord)
    } 
    tochange <- sphcoord[,2] < 0
    sphcoord[tochange,] <- .__opposite(sphcoord[tochange,])
    return(sphcoord)
}
