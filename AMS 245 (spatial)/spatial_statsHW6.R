library(fields)
library(maptools)
library(magic)

load("~/goes75.RData")
load("~/goes135.RData")

cloglog = function(data){
  return(log(-log(1-data)))
}

clogloginv = function(data){
  return(1-exp(-exp(data)))
}

y1 = data.frame(goes75)
y2 = data.frame(goes135)
Y = cloglog(c(y1$albedo, y2$albedo))

hist(cloglog(y1$albedo), main = 'Transformed G075 Albedo', xlab = 'cloglog(albedo)')
hist(cloglog(y2$albedo), main = 'Transformed G135 Albedo', xlab = 'cloglog(albedo)')

# Observation plot
{par(mfrow = c(1,2))
  quilt.plot(y1$lon, y1$lat, y1$albedo, main = 'GOES 075 Observations')
  map("world", add = TRUE)
  
  quilt.plot(y2$lon, y2$lat, y2$albedo, main = 'GOES 135 Observations')
  map("world", add = TRUE)
}


s1 = as.matrix(y1[,c('lon', 'lat')])
s2 = as.matrix(y2[,c('lon', 'lat')])
S = rbind(s1,s2)
# S = adiag(s1, s2)
X = adiag(cbind(1,s1), cbind(1,s2))

# knots
u = as.matrix(expand.grid(seq(-127, -89, by=2), seq(9.5, 35, by=2)))


## Find which points fall over land
data(wrld_simpl)
pts = SpatialPoints(u, proj4string=CRS(proj4string(wrld_simpl)))
u = pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords
# u = adiag(u,u)

m = nrow(u)

{par(mfrow = c(1,2))
  quilt.plot(y1$lon, y1$lat, cloglog(y1$albedo), main = 'GOES 075 Observations')
  map("world", add = TRUE)
  points(u[,1], u[,2], col = 'cadetblue')
  
  quilt.plot(y2$lon, y2$lat, cloglog(y2$albedo), main = 'GOES 135 Observations')
  map("world", add = TRUE)
  points(u[,1], u[,2], col = 'cadetblue')
}

# significant coefficients + intercept
mod = lm(Y ~ S[,1] + S[,2])
summary(mod)

# Prediction locations
pred_locs = as.matrix(expand.grid(seq(-126, -90, length=100), seq(10.5, 34, length=100)))
## Find which points fall over land
data(wrld_simpl)
pts = SpatialPoints(pred_locs, proj4string=CRS(proj4string(wrld_simpl)))
pred_locs = pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords

X_new = adiag(cbind(1, pred_locs[,1], pred_locs[,2]), cbind(1, pred_locs[,1], pred_locs[,2]))
s_new = adiag(pred_locs, pred_locs) 

# ------------------------- gmrf fit -----------------------------
source("/Users/SarJar/Documents/AMS 245 (spatial)/gmrf.R", chdir=TRUE)
source("/Users/SarJar/Documents/AMS 245 (spatial)/convolution.R", chdir=TRUE)
source("/Users/SarJar/Documents/AMS 245 (spatial)/mcmc.R", chdir=TRUE)


set.seed(1)
system.time(out2 <- gp_conv_gmrf_fit(Y, X, u=u, s=S, 
                         cs_v=.5, B=3000, burn=1000, print_freq=100))

b2 <- t(sapply(out2, function(o) o$beta))
colnames(b2) <- colnames(X)
plot.ts(b2)
colMeans(b2)

z <- t(sapply(out2, function(o) o$z))
colnames(z) <- paste0("w",1:nrow(u))

plot.ts(z[,1:4])
plot.ts(z[,5:8])
plot.ts(z[,9:12])

tau2g <- sapply(out2, function(o) o$tau2)
lam <- sapply(out2, function(o) o$lamz)
v2 <- sapply(out2, function(o) o$v)
plot.ts(v2)

plot.ts(cbind(lam, tau2g, v2))

post2 <- cbind(b2, lam, tau2g, v2)
post_summary(post2)

# -------------------------- predict -------------------------
system.time(pred2 <- gp_conv_gmrf_pred(Y, X, S, X_new, s_new, u, out2))

par(mfrow = c(2,2))
quilt.plot(pred_locs[,1], pred_locs[,2], rowMeans(pred2)[1:3270], main = 'cloglog G075 Predictions')
map('world',add=T)

quilt.plot(pred_locs[,1], pred_locs[,2], rowMeans(pred2)[3271:nrow(s_new)], main = 'cloglog G135 Predictions')
map('world',add=T)

quilt.plot(pred_locs[,1], pred_locs[,2], clogloginv(rowMeans(pred2))[1:3270], main = 'G075 Predictions')
map('world',add=T)

quilt.plot(pred_locs[,1], pred_locs[,2], clogloginv(rowMeans(pred2))[3271:nrow(s_new)], main = 'G135 Predictions')
map('world',add=T)


# -------------------------- predictive process mean --------------------
# don't know about this
# par(mfrow = c(2,2))
# # G075?
# quilt.plot(pred_locs[,1], pred_locs[,2], rowMeans(X_new%*%t(b2))[1:3270], main = 'G075 Predcitive Mean')
# map('world',add=T)
# 
# # G135?
# quilt.plot(pred_locs[,1], pred_locs[,2], rowMeans(X_new%*%t(b2))[3271:nrow(X_new%*%t(b2))], main = 'G135 Predcitive Mean')
# map('world',add=T)

# --------------------- predictive variance -------------------
# pred2_sd <- apply(pred2, 1, sd)
# 
# quilt.plot(pred_locs[,1], pred_locs[,2], pred2_sd[1:n], main = 'G075 Predictive SD')
# quilt.plot(pred_locs[,1], pred_locs[,2], pred2_sd[3271:6540], main = 'G135 Predictive SD')

#---------------------------- residuals -----------------------
# predict where we have observations
 system.time(pred.obs <- gp_conv_gmrf_pred(Y, X, S, X, S, u, out2))
# 
# resid <- apply(pred2, 2, function(x) x-Y)
# 
# par(mfrow = c(1,2))
# quilt.plot(s1[,1], s1[,2], rowMeans(resid)[1:1000], main = 'G075 Residuals')
# map('world', add=T)
# 
# quilt.plot(s2[,1], s2[,2], rowMeans(resid)[1001:2000], main = 'G135 Residuals')
# map('world', add=T)

# ----------------------------- accuracy ------------------------------

pred.int = apply(pred.obs, 1, function(x) quantile(x, c(0.025, 0.975)))

1-(sum(Y > pred.int[2,] | Y < pred.int[1,])/length(Y))

# 0.954 

ind=1:25
matplot(rbind(ind,ind), pred.int[,1:25], type="l", lty=1,col=1, xlab="Observation",ylab="cloglog weighted albedo", main = '95% Intervals for Weighted Albedo')
points(ind, Y[1:25], pch=19)

# ----------------------------- difference d(s)? ------------------------------

diff = cbind(pred_locs, rowMeans(pred2)[1:3270] - rowMeans(pred2)[3271:nrow(s_new)])
Y_diff = diff[sample(nrow(diff), 1000),]

X_diff  = cbind(1, Y_diff[,1:2]) 
S_diff = Y_diff[,1:2]

system.time(outdiff <- gp_conv_gmrf_fit(Y_diff[,3], X_diff, u=u, s=S_diff, 
                                        cs_v=.5, B=3000, burn=1000, print_freq=100))
d <- t(sapply(outdiff, function(o) o$z))
colnames(d) <- paste0("w",1:nrow(u))

d.pred = apply(d, 2, mean)

plot.nice(u, d.pred)

# --------------------------- spatial effects m(s) -----------------------
p = nrow(u)
k = ncol(X)

w.pred = apply(z, 2, mean)
w.pred.s = apply(z, 2, sd)

par(mfrow = c(1,2))
plot.nice(u, w.pred, main = 'True Common Albedo m(s)')
plot.nice(u, w.pred.s, main = 'Uncertainty SD m(s)')

mexico.ms = cbind(u, t(z))
names(mexico.ms[,1]) = 'lon'
names(mexico.ms[,2]) = 'lat'

#---------------- xiao ---------------- 

plot.nice(gridComField[,1:2], rowMeans(gridComField[,3:4002]))

plot.nice(gridDiffField[,1:2], rowMeans(gridDiffField[,3:4002]))

test = gridDiffField[sample(nrow(gridDiffField), 77),]

par(mfrow = c(1,2))
plot.nice(test[,1:2], rowMeans(test[,3:4002]), main = 'Discrepancy Between the Satellites d(s)')
plot.nice(test[,1:2], apply(test[,3:4002], 1, sd), main='Uncertainty SD d(s)')


mexico.ds = data.frame(cbind(test[,1], test[,2], test[,3:3002]))

names(mexico.ds) = c('lon', 'lat')
mexico.ds[1:10,1:10]
rownames(mexico.ds) <- NULL
save(mexico.ds, file = 'sarjar-diff-area-3.rda')

# ----------------------- combining predictions ---------------------
dim(difference.area.2)
difference.area.2[1:10,1:10]
meanfield[1:10,1:10]

one.mean = as.matrix(cbind(meanfield[,1:2], rowMeans(meanfield[,3:ncol(meanfield)])))
two.mean = cbind(mean.area.2[,1:2], rowMeans(mean.area.2[,3:ncol(mean.area.2)]))
three.mean = cbind(mexico.ms[,1:2], rowMeans(mexico.ms[,3:ncol(mexico.ms)]))
whole.mean = rbind(one.mean, two.mean, three.mean)

one.diff = as.matrix(cbind(field[,1:2], rowMeans(meanfield[,3:ncol(meanfield)])))

whole.diff = rbind(cbind(difference.area.2[,1:2], rowMeans(difference.area.2[,3:ncol(difference.area.2)])), 
                   cbind(mexico.ds[,1:2], rowMeans(mexico.ds[,3:ncol(mexico.ds)])))

plot.nice(whole.mean[,1:2], whole.mean[,3], main = 'True Common Albedo m(s)')
