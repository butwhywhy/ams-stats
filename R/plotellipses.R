# For generating graphs with the ellipses for hext and bootstrap methods applied to some
# fake simulated data. Used for poster
#source('simulations.R')
#source('ellipses.R')
#source('hext.R')
#source('constable.R')

plotellipses <- function(ams_tensor, sigma) {
    eigenparams <- eigen(ams_tensor)
    expected <- car2sph(t(eigenparams$vectors))
    setup <- AMSsetup()

    fake <- FakeMeasures(ams_tensor, n_measurements=15, error_dist=NormalErrorGenerator(sigma), setup)
    result_hext <- ams.hext(fake, setup)
    result_boot <- ams.constable(fake, setup)
    .__plot_ellipse(ellipse_hext=result_hext$eigenvectors$ellip1, ellipse_boot=result_boot$eigenvectors$ellip1, expectedCenter=expected[1,], add=F, symbols_index=1)
    .__plot_ellipse(ellipse_hext=result_hext$eigenvectors$ellip2, ellipse_boot=result_boot$eigenvectors$ellip2, expectedCenter=expected[2,], add=T, symbols_index=2)
    .__plot_ellipse(ellipse_hext=result_hext$eigenvectors$ellip3, ellipse_boot=result_boot$eigenvectors$ellip3, expectedCenter=expected[3,], add=T, symbols_index=3)
}

.__symbols <- data.frame(color=c('blue','green','magenta'), original_symbols=c(0,2,1), hext_symbols=c(15,17,16), bootstrap_symbols=c(22,24,21), hext_lines=c(1,1,1), bootstrap_lines=c(3,3,3), background=c('lightblue', 'lightgreen', 'lightpink'), stringsAsFactors=F)

.__plot_ellipse <- function(ellipse_hext, ellipse_boot, expectedCenter, symbols_index, add=F) {
    
    symbols <- .__symbols[symbols_index,]
    color <- .__symbols[symbols_index,'color']
    print(.__symbols[symbols_index,'hext_lines'])
    ellipse_plot(elip=ellipse_hext, add=add, line.col=color, lty=.__symbols[symbols_index,'hext_lines'])
    ellipse_plot(elip=ellipse_boot, add=T, line.col=color, divideHemispheres=F, lty=.__symbols[symbols_index,'bootstrap_lines'])
    lambert.plot(long=ellipse_hext$center[1], lat=ellipse_hext$center[2], divideHemispheres=F, add=T, point.col=color, point.symbols=.__symbols[symbols_index,'hext_symbols'])
    lambert.plot(long=ellipse_boot$center[1], lat=ellipse_boot$center[2], divideHemispheres=F, add=T, point.col=color, point.symbols=.__symbols[symbols_index,'bootstrap_symbols'], bg=.__symbols[symbols_index,'background'])
    lambert.plot(long=expectedCenter[1], lat=expectedCenter[2], divideHemispheres=F, add=T, point.col=color, point.symbols=.__symbols[symbols_index,'original_symbols'])
}
