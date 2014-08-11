#source('coordinates.R')

# An ellipse in the sphere is described by a list object 'elip' with constituent 
# elements 'elip$center', 'elip$eta', 'elip$zeta', 'elip$edir', 'elip$zdir'. 
# Of these, 'elip$center', 'elip$edir', 'elip$zdir', are vectors containing the
# longitud, latitud and radial coordinates of the center direction and the two
# semiangles, called eta and zeta. 'elip$eta' and 'elip$zeta' are the value of 
# the two semiangles. All angular magnitudes in radians.

spherical_ellipse <- function(centerDir, axis1Dir, axis2Dir, semiangle1, semiangle2) {
    center <- .__toNorth(centerDir)
    if (semiangle1 > semiangle2) {
        eta <- semiangle1
        edir <- .__toNorth(axis1Dir)
        zeta <- semiangle2
        zdir <- .__toNorth(axis2Dir) 
    } else {
        eta <- semiangle2
        edir <- .__toNorth(axis2Dir)
        zeta <- semiangle1
        zdir <- .__toNorth(axis1Dir) 
    }
    ellip <- list(center=center, eta=eta, zeta=zeta, edir=edir, zdir=zdir)
    class(ellip) <- c('spherical_ellipse', class(ellip))
    return(ellip)
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

# Plots the ellipse 'elip' in Lambert projection. The 'add' parameter indicates
# if the ellipse should be drawn over an existing plot or on a new one (default).
ellipse_plot <- function(elip, npoints=200, add=FALSE, line.col='red', ...) {
    points <- ellipse_points(elip, npoints)
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
