library(akima)
plot.nice = function(xy, z, zscale, nlevels = 20, ...){
  df = data.frame("x" = xy[,1], "y" = xy[,2], "z" = z)
  fld = with(df, interp(x = x, y = y, z = z))
  xlim = c(min(xy[,1]), max(xy[,1])+diff(range(xy[,1]))*0.2)
  ylim = range(xy[,2])
  offx = diff(range(xy[,1]))*0.2
  if (missing(zscale))
    zscale = range(fld$z, na.rm = TRUE)
  zlevels = seq(zscale[1], zscale[2], length = nlevels)
  plot(0, type='n', xlim = c(min(xy[,1]), max(xy[,1])+offx),
       ylim = range(xy[,2]),
       xlab = "Longitude", ylab = "Latitude", bty = 'n')
  .filled.contour(x = fld$x, y = fld$y, z = fld$z,
                  levels = zlevels,
                  col = tim.colors(nlevels))
  rect(rep(xlim[2] - offx*0.50 , nlevels), head(seq(ylim[1], ylim[2], length = nlevels+1), nlevels),
       rep(xlim[2] - offx*0.75, nlevels), tail(seq(ylim[1], ylim[2], length = nlevels+1), nlevels),
       col = tim.colors(nlevels))
  text(x = rep(xlim[2] - offx*0.25, 7), y = seq(ylim[1], ylim[2], length = 5),
       labels = round(seq(zscale[1], zscale[2], length = 5), 2))
  title(...)
  map("world", add = TRUE)
}
