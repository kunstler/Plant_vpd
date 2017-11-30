## Compare FF16 and FF16FvCB
compare_strategy_growth<- function(light = 1, ymax = 19){
require(plant)
derivs <- function(t, y, plant, env) {
 plant$ode_state <- y
 plant$compute_vars_phys(env)
 list(plant$ode_rates)
}

tt <- seq(0, 50, length.out=101)
env <- fixed_environment(light)

p <- scm_base_parameters("FF16")
s1 <- strategy(trait_matrix(1e-3,"narea"), p)
s2 <- strategy(trait_matrix(3e-3,"narea"), p)
s3 <- strategy(trait_matrix(9e-3,"narea"), p)
pl1 <- FF16_PlantPlus(s1)
pl2 <- FF16_PlantPlus(s2)
pl3 <- FF16_PlantPlus(s3)

y01 <- setNames(pl1$ode_state, pl1$ode_names)
yy1 <- deSolve::lsoda(y01, tt, derivs, pl1, env=env)
y02 <- setNames(pl2$ode_state, pl2$ode_names)
yy2 <- deSolve::lsoda(y02, tt, derivs, pl2, env=env)
y03 <- setNames(pl3$ode_state, pl3$ode_names)
yy3 <- deSolve::lsoda(y03, tt, derivs, pl3, env=env)
plot(height ~ time, yy1, type="l", ylim = c(0, ymax), main = paste("light level", light))
lines(yy2[, "time"], yy2[, "height"], lty = 2)
lines(yy3[, "time"], yy3[, "height"], lty = 3)


pF <- scm_base_parameters("FF16FvCB")

sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
plF1 <- FF16FvCB_PlantPlus(sF1)
plF2 <- FF16FvCB_PlantPlus(sF2)
plF3 <- FF16FvCB_PlantPlus(sF3)



yF01 <- setNames(plF1$ode_state, plF1$ode_names)
yyF1 <- deSolve::lsoda(yF01, tt, derivs, plF1, env=env)
yF02 <- setNames(plF2$ode_state, plF2$ode_names)
yyF2 <- deSolve::lsoda(yF02, tt, derivs, plF2, env=env)
yF03 <- setNames(plF3$ode_state, plF3$ode_names)
yyF3 <- deSolve::lsoda(yF03, tt, derivs, plF3, env=env)
lines(height ~ time, yyF1, type="l", lty = 1, col = "red")
lines(height ~ time, yyF2, type="l", lty = 2, col = "red")
lines(height ~ time, yyF3, type="l", lty = 3, col = "red")
legend("bottomright",
       c("Narea = 1e-3", "Narea = 3e-3", "Narea = 9e-3", "FF16", "FF16FvCB"),
       lty=c(1,2,3,1,1), col=c(rep("black", 3), "red", "black"), bty="n")

}


f_g<- function(x, pl) {
  env <- fixed_environment(x)
  pl$compute_vars_phys(env)
  pl$internals$height_dt
}

f_demo_size<- function(x, pl, sizeS) {
  env <- fixed_environment(x)
  pl$compute_vars_phys(env)
  pl$internals$height_dt
}


compare_strategy_light_growth<- function(){
require(plant)
openness <- seq(0, 1, length.out=101)

p <- trait_gradients_base_parameters(site_prod=0.3)

s1 <- strategy(trait_matrix(1e-3,"narea"), p)
s2 <- strategy(trait_matrix(3e-3,"narea"), p)
s3 <- strategy(trait_matrix(9e-3,"narea"), p)
pl1 <- FF16_PlantPlus(s1)
pl2 <- FF16_PlantPlus(s2)
pl3 <- FF16_PlantPlus(s3)


plot(openness, sapply(openness, f_g, pl1), type="l", xlim=c(0, 1),
     las=1, xlab="Canopy openness", ylab="Height growth rate (m / yr)",
     ylim = c(0, 1.4), lwd = 2)
lines(openness, sapply(openness, f_g, pl2), lty = 2, lwd = 2)
lines(openness, sapply(openness, f_g, pl3), lty = 3, lwd = 2)

p <- trait_gradients_base_parameters(site_prod=-0.3)

s1 <- strategy(trait_matrix(1e-3,"narea"), p)
s2 <- strategy(trait_matrix(3e-3,"narea"), p)
s3 <- strategy(trait_matrix(9e-3,"narea"), p)
pl1 <- FF16_PlantPlus(s1)
pl2 <- FF16_PlantPlus(s2)
pl3 <- FF16_PlantPlus(s3)
lines(openness, sapply(openness, f_g, pl1), lty = 1)
lines(openness, sapply(openness, f_g, pl2), lty = 2)
lines(openness, sapply(openness, f_g, pl3), lty = 3)


## FF16FvCB
pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))

sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
plF1 <- FF16FvCB_PlantPlus(sF1)
plF2 <- FF16FvCB_PlantPlus(sF2)
plF3 <- FF16FvCB_PlantPlus(sF3)

lines(openness, sapply(openness, f_g, plF1), lty = 1, col = "red", lwd = 2)
lines(openness, sapply(openness, f_g, plF2), lty = 2, col = "red", lwd = 2)
lines(openness, sapply(openness, f_g, plF3), lty = 3, col = "red", lwd = 2)

pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))

sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
plF1 <- FF16FvCB_PlantPlus(sF1)
plF2 <- FF16FvCB_PlantPlus(sF2)
plF3 <- FF16FvCB_PlantPlus(sF3)

lines(openness, sapply(openness, f_g, plF1), lty = 1, col = "red")
lines(openness, sapply(openness, f_g, plF2), lty = 2, col = "red")
lines(openness, sapply(openness, f_g, plF3), lty = 3, col = "red")

legend("topleft",
       c("Narea = 1e-3", "Narea = 3e-3", "Narea = 9e-3", "FF16", "FF16FvCB",
         "Low water availability", "High water availability"),
       lty=c(1,2,3,1,1, 1,1), lwd = c(1,1,1,1,1,1,2),
       col=c(rep("black", 3), "red", "black", "black", "black"), bty="n")

}


compare_strategy_growth_narea<- function(n_length = 100){
require(plant)

fl1 <- function(narea, lights){
  p <- trait_gradients_base_parameters(site_prod=-0.3)
  s <- strategy(trait_matrix(narea,"narea"), p)
  pl <- FF16_PlantPlus(s)
  pl$height <- 0.5
  sapply(lights, f_g, pl)
}

fh1 <- function(narea, lights){
  p <- trait_gradients_base_parameters(site_prod=0.3)
  s <- strategy(trait_matrix(narea,"narea"), p)
  pl <- FF16_PlantPlus(s)
  pl$height <- 0.5
  sapply(lights, f_g, pl)
}

fl2<- function(narea, lights){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))
  s <- strategy(trait_matrix(narea,"narea"), p)
  pl <- FF16FvCB_PlantPlus(s)
  pl$height <- 0.5
  sapply(lights, f_g, pl)
}

fh2<- function(narea, lights){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))
  s <- strategy(trait_matrix(narea,"narea"), p)
  pl <- FF16FvCB_PlantPlus(s)
  pl$height <- 0.5
  sapply(lights, f_g, pl)
}


fl2b <- function(n, lights){
tryCatch(fl2(n, lights), error=function(err) NA)
}

fh2b <- function(n, lights){
tryCatch(fh2(n, lights), error=function(err) NA)
}

v_narea1 <- 10^seq(log10(3e-4),log10(3e-2), length.out = n_length/2)
v_narea2 <- seq(3e-4,3e-2, length.out = n_length/2)
v_narea <-  sort(c(v_narea1, v_narea2))
v_lg <- t(sapply(v_narea, fl1, lights = c(0.25, 1.0)))
v_hg <- t(sapply(v_narea, fh1, lights = c(0.25, 1.0)))
v_lg_F <- t(sapply(v_narea, fl2b, lights = c(0.25, 1.0)))
v_hg_F <- t(sapply(v_narea, fh2b, lights = c(0.25, 1.0)))

plot(v_narea, v_hg[, 2], type="l", ylim = range(v_lg, v_lg_F, v_hg, v_hg_F),
     las=1, xlab="Nitrogen per area", ylab="Height growth rate", log = "x")
lines(v_narea, v_hg[, 1], lty = 2)
lines(v_narea, v_lg[, 1], lty = 2, col = "grey")
lines(v_narea, v_lg[, 2], lty = 1, col = "grey")
lines(v_narea, v_hg_F[, 2], col = "red")
lines(v_narea, v_hg_F[, 1], col = "red", lty = 2)
lines(v_narea, v_lg_F[, 2], col = "pink")
lines(v_narea, v_lg_F[, 1], col = "pink", lty = 2)

legend("topleft",
       c("High light", "Low light", "FF16 site_prod = 0.3", "FF16 site_prod = -0.3", "FF16FvCB vpd = 0",
         "FF16FvCB vpd = 3"),
       lty=c(1,2,1,1, 1,1),
       col=c("black", "black", "black", "grey", "red", "pink"), bty="n")

}


## growth in function of Narea
compare_strategy_growth_narea_vpd<- function(n_length = 40){
require(plant)

f1 <- function(n, vpd = 0.3){
  p <- trait_gradients_base_parameters(site_prod=vpd)
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16_PlantPlus(s)
  pl$height <- 0.5
  f_g(1.0, pl)
}


f2<- function(n, vpd = 0){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = vpd))
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16FvCB_PlantPlus(s)
  pl$height <- 0.5
  f_g(1.0, pl)
}

f2b <- function(n, vpd){
tryCatch(f2(n, vpd), error=function(err) NA)
}

v_narea <- seq(1e-3,3e-2, length.out = n_length)
v_vpd1 <- seq(-0.3, 0.3, length.out = 10)
v_vpd2 <- seq(0, 3, length.out = 10)

m_growth <- matrix(NA, nrow = 10, ncol = n_length)
m_growth_F<- matrix(NA, nrow = 10, ncol = n_length)

for (i in 1:10) {
m_growth[i, ]<- sapply(v_narea, f1, vpd = v_vpd1[i])
m_growth_F[i, ]<- sapply(v_narea, f2b, vpd = v_vpd2[i])
}
return(list(m_growth, m_growth_F))
}

plot_strategy_growth_narea_vpd<- function(l_data){
m_lcp <-  l_data[[1]]
m_lcp_F <- l_data[[2]]
par(mfrow = c(1, 2))
image(m_lcp[rev(seq_len(nrow(m_lcp))), ], xlab = "water stress", ylab = "Narea")
contour(m_lcp[rev(seq_len(nrow(m_lcp))), ], add =TRUE)
image(m_lcp_F, xlab = "water stress", ylab = "Narea")
contour(m_lcp_F, add =TRUE)
}


# much less accurate than what is normaly done in Plant but just while lcp_whole_plant doesn't work ...
lcp_georges <- function(pl){
  openness <- seq(0, 1, length.out=501)
  v <- sapply(openness, f_g, pl)
  lcp <- mean(openness[c(max(which(v < 1e-10)), min(which(v > 1e-10)))])
  return(lcp)
}


compare_strategy_lcp<- function(){
require(plant)
p <- trait_gradients_base_parameters(site_prod=0.3)
s1 <- strategy(trait_matrix(1e-3,"narea"), p)
s2 <- strategy(trait_matrix(3e-3,"narea"), p)
s3 <- strategy(trait_matrix(9e-3,"narea"), p)
pl1 <- FF16_PlantPlus(s1)
pl2 <- FF16_PlantPlus(s2)
pl3 <- FF16_PlantPlus(s3)

openness <- seq(0, 1, length.out=501)

lcp1 <- lcp_georges(pl1)
lcp2 <- lcp_georges(pl2)
lcp3 <- lcp_georges(pl3)

x <- c(lcp1, openness[openness > lcp1])
plot(x, sapply(x, f_g, pl1), type="l", xlim=c(0, 1), ylim = c(0, 1.5),
     las=1, xlab="Canopy openness", ylab="Height growth rate (m / yr)", lwd = 2)
points(lcp1, 0.0, pch=19)
x <- c(lcp2, openness[openness > lcp2])
lines(x, sapply(x, f_g, pl2), lty = 2, lwd = 2)
points(lcp2, 0.0, pch=19)
x <- c(lcp3, openness[openness > lcp3])
lines(x, sapply(x, f_g, pl3), lty = 3, lwd = 2)
points(lcp3, 0.0, pch=19)


p <- trait_gradients_base_parameters(site_prod=-0.3)
s1 <- strategy(trait_matrix(1e-3,"narea"), p)
s2 <- strategy(trait_matrix(3e-3,"narea"), p)
s3 <- strategy(trait_matrix(9e-3,"narea"), p)
pl1 <- FF16_PlantPlus(s1)
pl2 <- FF16_PlantPlus(s2)
pl3 <- FF16_PlantPlus(s3)

openness <- seq(0, 1, length.out=51)

lcp1 <- lcp_georges(pl1)
lcp2 <- lcp_georges(pl2)
lcp3 <- lcp_georges(pl3)

x <- c(lcp1, openness[openness > lcp1])
lines(x, sapply(x, f_g, pl1), lty = 1, col = "grey")
points(lcp1, 0.0, pch=19, col = "grey")
x <- c(lcp2, openness[openness > lcp2])
lines(x, sapply(x, f_g, pl2), lty = 2, col = "grey")
points(lcp2, 0.0, pch=19, col = "grey")
x <- c(lcp3, openness[openness > lcp3])
lines(x, sapply(x, f_g, pl3), lty = 3, col = "grey")
points(lcp3, 0.0, pch=19, col = "grey")

## FF16FvCB
pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))

sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
plF1 <- FF16FvCB_PlantPlus(sF1)
plF2 <- FF16FvCB_PlantPlus(sF2)
plF3 <- FF16FvCB_PlantPlus(sF3)

lcpF1 <- lcp_georges(plF1)
lcpF2 <- lcp_georges(plF2)
lcpF3 <- lcp_georges(plF3)

x <- c(lcpF1, openness[openness > lcpF1])
lines(x, sapply(x, f_g, plF1), col = 'red', lwd = 2)
points(lcpF1, 0.0, pch=19, col = "red")
x <- c(lcpF2, openness[openness > lcpF2])
lines(x, sapply(x, f_g, plF2), lty = 2, col = 'red', lwd = 2)
points(lcpF2, 0.0, pch=19, col = "red")
x <- c(lcpF3, openness[openness > lcpF3])
lines(x, sapply(x, f_g, plF3), lty = 3, col = 'red', lwd = 2)
points(lcpF3, 0.0, pch=19, col = "red")

pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))

sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
plF1 <- FF16FvCB_PlantPlus(sF1)
plF2 <- FF16FvCB_PlantPlus(sF2)
plF3 <- FF16FvCB_PlantPlus(sF3)

lcpF1 <- lcp_georges(plF1)
lcpF2 <- lcp_georges(plF2)
lcpF3 <- lcp_georges(plF3)

x <- c(lcpF1, openness[openness > lcpF1])
lines(x, sapply(x, f_g, plF1), col = 'pink')
points(lcpF1, 0.0, pch=19, col = "pink")
x <- c(lcpF2, openness[openness > lcpF2])
lines(x, sapply(x, f_g, plF2), lty = 2, col = 'pink')
points(lcpF2, 0.0, pch=19, col = "pink")
x <- c(lcpF3, openness[openness > lcpF3])
lines(x, sapply(x, f_g, plF3), lty = 3, col = 'pink')
points(lcpF3, 0.0, pch=19, col = "pink")

legend("topleft",
       c("Narea = 1e-3", "Narea = 3e-3", "Narea = 9e-3",
         "FF16 site_prod = 0.3", "FF16 site_prod = -0.3", "FF16FvCB vpd = 0",
         "FF16FvCB vpd = 3"),
       lty=c(1,2, 3,1,1, 1,1),
       col=c("black", "black", "black","black", "grey", "red", "pink"), bty="n")

}


## lcp in function of Narea
compare_strategy_lcp_narea<- function(n_length = 100){
require(plant)

fl1 <- function(n){
  p <- trait_gradients_base_parameters(site_prod=-0.3)
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16_PlantPlus(s)
  pl$height <- 0.5
  lcp_georges(pl)
}

fh1 <- function(n){
  p <- trait_gradients_base_parameters(site_prod=0.3)
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16_PlantPlus(s)
  pl$height <- 0.5
  lcp_georges(pl)
}

fl2<- function(n){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16FvCB_PlantPlus(s)
  pl$height <- 0.5
  lcp_georges(pl)
}

fh2<- function(n){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16FvCB_PlantPlus(s)
  pl$height <- 0.5
  lcp_georges(pl)
}
fl2b <- function(n){
tryCatch(fl2(n), error=function(err) NA)
}
fh2b <- function(n){
tryCatch(fh2(n), error=function(err) NA)
}

v_narea1 <- 10^seq(log10(1e-3),log10(3e-2), length.out = n_length/2)
v_narea2 <- seq(1e-3,3e-2, length.out = n_length/2)
v_narea <-  sort(c(v_narea1, v_narea2))
v_lcp_l<- sapply(v_narea, fl1)
v_lcp_h<- sapply(v_narea, fh1)
v_lcp_l_F <- sapply(v_narea, fl2b)
v_lcp_h_F <- sapply(v_narea, fh2b)

plot(v_narea, v_lcp_h, type="l", ylim = c(0, 1),
     las=1, xlab="Nitrogen per area", ylab="WP LCP",
     log = "x")
lines(v_narea, v_lcp_l, col = "grey")
lines(v_narea, v_lcp_l_F, col = "pink")
lines(v_narea, v_lcp_h_F, col = "red")
legend("topleft",
       c("FF16 site_prod = 0.3", "FF16 site_prod = -0.3", "FF16FvCB vpd = 0",
         "FF16FvCB vpd = 3"),
       lty=c(1,1, 1,1),
       col=c("black", "grey", "red", "pink"), bty="n")
}



## lcp in function of Narea
compare_strategy_lcp_narea_vpd<- function(n_length = 40){
require(plant)

f1 <- function(n, vpd = 0.3){
  p <- trait_gradients_base_parameters(site_prod=vpd)
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16_PlantPlus(s)
  pl$height <- 0.5
  lcp_georges(pl)
}


f2<- function(n, vpd = 0){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = vpd))
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16FvCB_PlantPlus(s)
  pl$height <- 0.5
  lcp_georges(pl)
}

f2b <- function(n, vpd){
tryCatch(f2(n, vpd), error=function(err) NA)
}

v_narea <- seq(1e-3,3e-2, length.out = n_length)
v_vpd1 <- seq(-0.3, 0.3, length.out = n_length)
v_vpd2 <- seq(0, 3, length.out = n_length)

m_lcp <- matrix(NA, nrow = n_length, ncol = n_length)
m_lcp_F<- matrix(NA, nrow = n_length, ncol = n_length)

for (i in 1:n_length) {
m_lcp[i, ]<- sapply(v_narea, f1, vpd = v_vpd1[i])
m_lcp_F[i, ]<- sapply(v_narea, f2b, vpd = v_vpd2[i])
}
return(list(m_lcp, m_lcp_F))
}


plot_strategy_lcp_narea_vpd<- function(l_data){
m_lcp <-  l_data[[1]]
m_lcp_F <- l_data[[2]]
par(mfrow = c(1, 2))
image(m_lcp[rev(seq_len(nrow(m_lcp))), ], xlab = "water stress", ylab = "Narea")
contour(m_lcp[rev(seq_len(nrow(m_lcp))), ], add =TRUE)
image(m_lcp_F, xlab = "water stress", ylab = "Narea")
contour(m_lcp_F, add =TRUE)
}




## look at patch dynamcis for two different N value

compare_strategy_patch<- function(v_narea = c(3e-3, 5e-3, 9e-3)){
require(plant)

ctrl <- equilibrium_verbose(fast_control())
ctrl$schedule_eps <- 0.002
ctrl$equilibrium_eps <- 1e-4
ctrl$equilibrium_solver_name <- "iteration"
ctrl$equilibrium_verbose <-  TRUE
ctrl$equilibrium_nsteps <- 80

rel <- function(x, xmin) {
  x[x < xmin] <- xmin
  xmax <- max(x, na.rm=TRUE)
  (x - xmin) / (xmax - xmin)
}

p0 <- trait_gradients_base_parameters(site_prod=0.3)
p0$control <- ctrl
p0$disturbance_mean_interval <- 30.0

## First, with three species:
p <- expand_parameters(trait_matrix(v_narea, "narea"), p0, FALSE)

p_eq <- equilibrium_seed_rain(p)

## Then collect the patch-level dynamics:
data <- run_scm_collect(p_eq)

tt <- data$time
h1 <- data$species[[1]]["height", , ]
h2 <- data$species[[2]]["height", , ]
h3 <- data$species[[3]]["height", , ]

cols <- c("#e41a1c", "#377eb8", "#4daf4a")

## Relativise the log densities onto (-4, max)
d1 <- data$species[[1]]["log_density", , ]
d2 <- data$species[[2]]["log_density", , ]
d3 <- data$species[[3]]["log_density", , ]


rd1 <- rel(d1, -4)
rd2 <- rel(d2, -4)
rd3 <- rel(d3, -4)

n <- length(tt)
x1 <- matrix(rep(tt, ncol(h1)), nrow(h1))
x2 <- matrix(rep(tt, ncol(h2)), nrow(h2))
x3 <- matrix(rep(tt, ncol(h3)), nrow(h3))
col1 <- matrix(make_transparent(cols[[1]], rd1), nrow(d1))
col2 <- matrix(make_transparent(cols[[2]], rd2), nrow(d2))
col3 <- matrix(make_transparent(cols[[3]], rd3), nrow(d3))
par(mfrow = c(2,2))
plot(NA, xlim=range(tt),
     las=1, xlab="Time (years)", ylab="Cohort height (m)",
     ylim = range(h1,h2,h3, na.rm = TRUE), main = "prod = 0.3")
segments(x1[-1, ], h1[-1, ], x1[-n, ], h1[-n, ], col=col1[-n, ], lend="butt")
segments(x2[-1, ], h2[-1, ], x2[-n, ], h2[-n, ], col=col2[-n, ], lend="butt")
segments(x3[-1, ], h3[-1, ], x3[-n, ], h3[-n, ], col=col3[-n, ], lend="butt")


p0 <- trait_gradients_base_parameters(site_prod=-0.3)
p0$control <- ctrl
p0$disturbance_mean_interval <- 30.0

## First, with three species:
p <- expand_parameters(trait_matrix(v_narea, "narea"), p0, FALSE)

p_eq <- equilibrium_seed_rain(p)

## Then collect the patch-level dynamics:
data <- run_scm_collect(p_eq)

tt <- data$time
h1 <- data$species[[1]]["height", , ]
h2 <- data$species[[2]]["height", , ]
h3 <- data$species[[3]]["height", , ]

cols <- c("#e41a1c", "#377eb8", "#4daf4a")

## Relativise the log densities onto (-4, max)
d1 <- data$species[[1]]["log_density", , ]
d2 <- data$species[[2]]["log_density", , ]
d3 <- data$species[[3]]["log_density", , ]

rd1 <- rel(d1, -4)
rd2 <- rel(d2, -4)
rd3 <- rel(d3, -4)

n <- length(tt)
x1 <- matrix(rep(tt, ncol(h1)), nrow(h1))
x2 <- matrix(rep(tt, ncol(h2)), nrow(h2))
x3 <- matrix(rep(tt, ncol(h3)), nrow(h3))
col1 <- matrix(make_transparent(cols[[1]], rd1), nrow(d1))
col2 <- matrix(make_transparent(cols[[2]], rd2), nrow(d2))
col3 <- matrix(make_transparent(cols[[3]], rd3), nrow(d3))
plot(NA, xlim=range(tt),
     las=1, xlab="Time (years)", ylab="Cohort height (m)",
     ylim = range(h1,h2,h3, na.rm = TRUE), main = "prod = -0.3")
segments(x1[-1, ], h1[-1, ], x1[-n, ], h1[-n, ], col=col1[-n, ], lend="butt")
segments(x2[-1, ], h2[-1, ], x2[-n, ], h2[-n, ], col=col2[-n, ], lend="butt")
segments(x3[-1, ], h3[-1, ], x3[-n, ], h3[-n, ], col=col3[-n, ], lend="butt")


## FvCB

pF0 <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))
pF0$control <-  ctrl
pF0$disturbance_mean_interval <- 30.0

## First, with three species:
pF <- expand_parameters(trait_matrix(v_narea, "narea"), pF0, FALSE)

pF_eq <- equilibrium_seed_rain(pF)

## Then collect the patch-level dynamics:
dataF <- run_scm_collect(pF_eq)

tt <- dataF$time
h1 <- dataF$species[[1]]["height", , ]
h2 <- dataF$species[[2]]["height", , ]
h3 <- dataF$species[[3]]["height", , ]

cols <- c("#e41a1c", "#377eb8", "#4daf4a")

## Relativise the log densities onto (-4, max)
d1 <- dataF$species[[1]]["log_density", , ]
d2 <- dataF$species[[2]]["log_density", , ]
d3 <- dataF$species[[3]]["log_density", , ]

rd1 <- rel(d1, -4)
rd2 <- rel(d2, -4)
rd3 <- rel(d3, -4)
n <- length(tt)
x1 <- matrix(rep(tt, ncol(h1)), nrow(h1))
x2 <- matrix(rep(tt, ncol(h2)), nrow(h2))
x3 <- matrix(rep(tt, ncol(h3)), nrow(h3))
col1 <- matrix(make_transparent(cols[[1]], rd1), nrow(d1))
col2 <- matrix(make_transparent(cols[[2]], rd2), nrow(d2))
col3 <- matrix(make_transparent(cols[[3]], rd3), nrow(d3))
plot(NA, xlim=range(tt),
     las=1, xlab="Time (years)", ylab="Cohort height (m)",
     ylim = range(h1,h2,h3, na.rm = TRUE), main = "FF16FvCB vpd = 0")
segments(x1[-1,], h1[-1, ], x1[-n, ], h1[-n, ], col=col1[-n, ], lend="butt")
segments(x2[-1, ], h2[-1, ], x2[-n, ], h2[-n, ], col=col2[-n, ], lend="butt")
segments(x3[-1, ], h3[-1, ], x3[-n, ], h3[-n, ], col=col3[-n, ], lend="butt")


## high
pF0 <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))
pF0$control <-  ctrl
pF0$disturbance_mean_interval <- 30.0

## First, with three species:
pF <- expand_parameters(trait_matrix(v_narea, "narea"), pF0, FALSE)

pF_eq <- equilibrium_seed_rain(pF)

## Then collect the patch-level dynamics:
dataF <- run_scm_collect(pF_eq)

tt <- dataF$time
h1 <- dataF$species[[1]]["height", , ]
h2 <- dataF$species[[2]]["height", , ]
h3 <- dataF$species[[3]]["height", , ]

cols <- c("#e41a1c", "#377eb8", "#4daf4a")

## Relativise the log densities onto (-4, max)
d1 <- dataF$species[[1]]["log_density", , ]
d2 <- dataF$species[[2]]["log_density", , ]
d3 <- dataF$species[[3]]["log_density", , ]

rd1 <- rel(d1, -4)
rd2 <- rel(d2, -4)
rd3 <- rel(d3, -4)
n <- length(tt)
x1 <- matrix(rep(tt, ncol(h1)), nrow(h1))
x2 <- matrix(rep(tt, ncol(h2)), nrow(h2))
x3 <- matrix(rep(tt, ncol(h3)), nrow(h3))
col1 <- matrix(make_transparent(cols[[1]], rd1), nrow(d1))
col2 <- matrix(make_transparent(cols[[2]], rd2), nrow(d2))
col3 <- matrix(make_transparent(cols[[3]], rd3), nrow(d3))
plot(NA, xlim=range(tt),
     las=1, xlab="Time (years)", ylab="Cohort height (m)",
     ylim = range(h1,h2,h3, na.rm = TRUE), main = "FF16FvCB vpd = 3")
segments(x1[-1,], h1[-1, ], x1[-n, ], h1[-n, ], col=col1[-n, ], lend="butt")
segments(x2[-1, ], h2[-1, ], x2[-n, ], h2[-n, ], col=col2[-n, ], lend="butt")
segments(x3[-1, ], h3[-1, ], x3[-n, ], h3[-n, ], col=col3[-n, ], lend="butt")


legend("right", legend = c("Narea = 1e-3", "Narea = 5e-3", "Narea = 9e-3"), lty=1, col = cols, bty="n")
}

#####################################################
#####################################################
### COMPARE FF16 and FF16FvCB along stress graident vpd

demo_ver_site_lma_light <- function(fun_param, .site_prod = 1,  lma = 1,
                                   light = 1.0){
hyper <- fun_param(site_prod = .site_prod)
s <- strategy(trait_matrix(lma, "lma"), hyper)
env<- fixed_environment(light)
p<- FF16_PlantPlus(s)
p$compute_vars_phys(env)
# compute demo
tt <- seq(0, 70, length.out = 101)
res<-  grow_plant_to_time(p, tt, env)
res<-  grow_plant_to_size(pl, c(1, 20), "height", env)
dhdt<- sapply(res$plant, function(x) x$internals$height_dt)
mat_demo <- res$state[, c('height', 'mortality', 'fecundity')]
mat_demo <- cbind(time = tt, mat_demo, dhdt = dhdt)
return(mat_demo)
}

demo_ver_siteS_lmaS_lightS<- function(fun_param,
                                      site_prodS = c(-0.4, 0, 0.4),
                                      lmaS = c(0.5, 1 , 1.5),
                                      lightS = c(0.2, 1.0)){
names(site_prodS) <- paste0('site_prod', site_prodS)
names(lmaS) <- paste0('lma', lmaS)
names(lightS) <- paste0('light', lightS)
res <- array(NA, dim = c(101, 5, length(site_prodS),
                         length(lmaS), length(lightS)),
             dimnames = list(NULL,
                             c('time', 'height', 'mortality',
                               'fecundity', 'dhdt'),
                               names(site_prodS),
                               names(lmaS),
                               names(lightS)))
for(site_prod in names(site_prodS)){
    for(lma in names(lmaS)){
        for(light in names(lightS)){
           res[ , ,
              site_prod,
              lma,
              light] <- demo_ver_site_lma_light(fun_param,
                                            .site_prod = site_prodS[site_prod],
                                            lma = lmaS[lma],
                                            light = lightS[light])
        }
     }
 }
return(res)
}


plot_demo_var <- function(res, site_prod_n, var, lmaS_n, lightS_n){
    require(RColorBrewer)
    cols <- brewer.pal(3,'Blues')
   plot(res[ , 'time', site_prod_n,
               lmaS_n[1],
               lightS_n[2]],
        res[ , var, site_prod_n,
               lmaS_n[2],
               lightS_n[1]], type = "l", lwd = 2,
        ylim = range(res[ -1 , var, , , ]),
        col = cols[2],
        xlab = "time", ylab = var)
   lines(res[ , 'time', site_prod_n,
               lmaS_n[1],
               lightS_n[1]],
         res[ , var, site_prod_n,
               lmaS_n[1],
               lightS_n[1]], col = cols[1], lwd = 2)
   lines(res[ , 'time', site_prod_n,
               lmaS_n[3],
               lightS_n[1]],
         res[ , var, site_prod_n,
               lmaS_n[3],
               lightS_n[1]], col = cols[3], lwd = 2)
# low light
   lines(res[ , 'time', site_prod_n,
               lmaS_n[1],
               lightS_n[2]],
         res[ , var, site_prod_n,
               lmaS_n[1],
               lightS_n[2]], col = cols[1], lty = 2, lwd = 2)
   lines(res[ , 'time', site_prod_n,
               lmaS_n[2],
               lightS_n[2]],
         res[ , var, site_prod_n,
               lmaS_n[2],
               lightS_n[2]], col = cols[2], lty = 2, lwd = 2)
   lines(res[ , 'time', site_prod_n,
               lmaS_n[3],
               lightS_n[2]],
         res[ , var, site_prod_n,
               lmaS_n[3],
               lightS_n[2]], col = cols[3], lty = 2, lwd = 2)
}

plot_demo_site_lma_light <- function(vars, res, site_prodS, lmaS, lightS){
site_prodS_n<- paste0('site_prod', site_prodS)
lmaS_n<- paste0('lma', lmaS)
lightS_n<- paste0('light', lightS)

par(mfrow = c(length(site_prodS_n), length(vars)))
  for (var in vars){
   for (site_prod_n in site_prodS_n){
       plot_demo_var(res, site_prod_n, var, lmaS_n, lightS_n)
   }
  }
}


## res <- demo_ver_siteS_lmaS_lightS(trait_gradients_base_parameters)
## plot_demo_site_lma_light(vars = c("mortality", "height", "fecundity"),
##                          res,
##                          site_prodS = c(-0.4, 0, 0.4),
##                          lmaS = c(0.5, 1 , 1.5),
##                          lightS = c(0.2, 1.0))


lcp_lma_site_prod_height <- function(fun_param, lmaS, .site_prod = 0,
                                    height = 0.3920458, narea=1.87e-3){
  hyper <- fun_param(site_prod = .site_prod)
  lcp_vec <- rep(NA, length(lmaS))
  for (i in seq_len(length(lmaS))){
      s <- strategy(trait_matrix(c(lmaS[i], narea), c("lma", "narea")), hyper)
      p<- FF16_PlantPlus(s)
      p$height <- height
      lcp_vec[i]<- lcp_whole_plant(p)
  }
 return(lcp_vec)
}

dhdt_lma_site_prod_height<- function(fun_param, lmaS, .site_prod = 0,
                                     height = 0.3920458, narea=1.87e-3){
  hyper <- fun_param(site_prod = .site_prod)
  dhdt_vec <- rep(NA, length(lmaS))
  for (i in seq_len(length(lmaS))){
      s <- strategy(trait_matrix(c(lmaS[i], narea), c("lma", "narea")), hyper)
      p<- FF16_PlantPlus(s)
      p$height <- height
      env <- fixed_environment(1)
      p$compute_vars_phys(env)
      dhdt_vec[i] <- p$internals$height_dt
  }
 return(dhdt_vec)
}



lcp_narea_site_prod_height <- function(fun_param, nareaS, .site_prod = 0,
                                    height = 0.3920458, lma= 0.1978791){
  hyper <- fun_param(site_prod = .site_prod)
  lcp_vec <- rep(NA, length(nareaS))
  for (i in seq_len(length(nareaS))){
      s <- strategy(trait_matrix(c(lma, nareaS[i]), c("lma", "narea")), hyper)
      p<- FF16_PlantPlus(s)
      p$height <- height
      lcp_vec[i]<- lcp_whole_plant(p)
  }
 return(lcp_vec)
}

dhdt_narea_site_prod_height<- function(fun_param, nareaS, .site_prod = 0,
                                     height = 0.3920458, lma=0.1978791){
  hyper <- fun_param(site_prod = .site_prod)
  dhdt_vec <- rep(NA, length(nareaS))
  for (i in seq_len(length(nareaS))){
      s <- strategy(trait_matrix(c(lma, nareaS[i]), c("lma", "narea")), hyper)
      p<- FF16_PlantPlus(s)
      p$height <- height
      env <- fixed_environment(1)
      p$compute_vars_phys(env)
      dhdt_vec[i] <- p$internals$height_dt
  }
 return(dhdt_vec)
}


lcp_dhdt_lma_site_prod_height<- function(fun_param, lmaS, .site_prod = 0,
                                    height = 1, narea=1.87e-3){
  dhdt_vec <- dhdt_lma_site_prod_height(fun_param, lmaS, .site_prod,
                                        height, narea)
  lcp_vec <- lcp_lma_site_prod_height(fun_param, lmaS, .site_prod,
                                      height, narea)
  return(data.frame(lcp = lcp_vec,
                    dhdt = dhdt_vec))
}


plot_lcp_version <- function(param_slope_P, param_slope_TP){
 lmaS <- seq(from = 0.05, to = 0.3, length.out = 100)
 par(mfrow = c(1, 3))
 for (site in c(-0.3, 0, 0.3)){
 lcp_base <- lcp_lma_site_prod_height(trait_gradients_base_parameters,
                                     lmaS, .site_prod = site)
 lcp_elev <- lcp_lma_site_prod_height(trait_gradients_elev_parameters,
                                     lmaS, .site_prod = site)
 lcp_slope <- lcp_lma_site_prod_height(trait_gradients_slope_parameters,
                                     lmaS, .site_prod = site)
 plot(lmaS, lcp_base, type = 'l',
      xlab = 'lma',
      ylab = 'Whole-plant light compensation point',
      ylim = c(0, 1))
 lines(lmaS, lcp_elev, col = 'red')
 lines(lmaS, lcp_slope, col = 'green')
 }
}

plot_lma_shade_tradeoff <-  function(lmaS = seq(from = 0.06, to = 0.27, length.out = 100)){
df_dhdt_lcp <- lcp_dhdt_lma_site_prod_height(trait_gradients_base_parameters,
                                             lmaS = lmaS)
par(mfcol = c(2,2))
plot(lmaS, 1 - df_dhdt_lcp$lcp, type = 'l',
     xlab = 'LMA',
     ylab = '1 - light compensation point (%)', cex.lab = 1.6)
plot(lmaS, df_dhdt_lcp$dhdt, type = 'l',
     ylab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
     xlab = 'LMA', cex.lab = 1.6)
plot(df_dhdt_lcp$dhdt, 1 - df_dhdt_lcp$lcp, type = 'l',
     xlab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
     ylab = '1 - light compensation point (%)', cex.lab = 1.6)
}


photo_max_LeafN_Stress<- function(
                                narea=1.87e-3,
                                narea_0=1.87e-3,
                                lma_0=0.1978791,
                                B_kl1=0.4565855,
                                B_kl2=1.71,
                                rho_0=608.0,
                                B_dI1=0.01,
                                B_dI2=0.0,
                                B_ks1=0.2,
                                B_ks2=0.0,
                                B_rs1=4012.0,
                                B_rb1=2.0*4012.0,
                                B_f1 =3.0,
                                B_lf1=5120.738 * 1.87e-3 * 24 * 3600 / 1e+06,
                                B_lf2=0.5,
                                B_lf3=0.04,
                                B_lf4=21000,
                                B_lf5=1,
                                k_I=0.5,
                                latitude=0,
                                site_prod=0) {

    ## rho / sapwood respiration relationship:

    ## Narea, photosynthesis, respiration

    assimilation_rectangular_hyperbolae <- function(I, Amax, theta, QY) {
      x <- QY * I + Amax
      (x - sqrt(x^2 - 4 * theta * QY * I * Amax)) / (2 * theta)
    }

    ## Photosynthesis  [mol CO2 / m2 / yr]
    approximate_annual_assimilation <- function(narea, latitude) {
      E <- seq(0, 1, by=0.02)
      ## Only integrate over half year, as solar path is symmetrical
      D <- seq(0, 365/2, length.out = 10000)
      I <- plant:::PAR_given_solar_angle(plant:::solar_angle(D,
                                                             latitude = abs(latitude)))

      Amax <- (1+site_prod)*B_lf1 * (narea/narea_0) ^  B_lf5
      theta <- B_lf2
      QY <- B_lf3

      AA <- NA * E

      for (i in seq_len(length(E))) {
        AA[i] <- 2 * plant:::trapezium(D, assimilation_rectangular_hyperbolae(
                                    k_I * I * E[i], Amax, theta, QY))
      }
      if(all(diff(AA) < 1E-8)) {
        # line fitting will fail if all have are zero, or potentially same value
        ret <- c(last(AA), 0)
        names(ret) <- c("p1","p2")
      } else {
        fit <- nls(AA ~ p1 * E/(p2 + E), data.frame(E = E, AA = AA),
                   start = list(p1 = 100, p2 = 0.2))
        ret <- coef(fit)
      }
      ret
    }
    # This needed in case narea has length zero, in which case trapezium fails
    a_p1 <- a_p2 <- 0 * narea
    ## TODO: Remove the 0.5 hardcoded default for k_I here, and deal
    ## with this more nicely.
    if (length(narea) > 0 || k_I != 0.5) {
      i <- match(narea, unique(narea))
      y <- vapply(unique(narea), approximate_annual_assimilation,
                  numeric(2), latitude)
      a_p1  <- y["p1", i]
      a_p2  <- y["p2", i]
    }


  y <- vapply(unique(narea), approximate_annual_assimilation,
              numeric(2), latitude)
  return(y["p1", 1])
}



## plot strategy at different vpd


data_contour  <- function(nareaS = seq(from = 0.2, to = 10, length.out = 100)/1000,
                           site_prodS = seq(from = -0.6, to = 0.6, length.out = 100)){
   res_photo <- matrix(NA, nrow = 100, ncol = 100)
   rownames(res_photo) <- paste0('narea', nareaS)
   colnames(res_photo) <- paste0('site_prod', site_prodS)
   for(narea in nareaS){
     for(site_prod in site_prodS){
      res <- photo_max_LeafN_Stress(narea = narea, site_prod = site_prod)
      res_photo[paste0('narea', narea),
                paste0('site_prod', site_prod)] <-  res
     }
   }
   return(res_photo)
}

plot_contour_waterStress_LeafN <- function(res_photo,
                                           nareaS = seq(from = 0.2, to = 10, length.out = 100)/1000,
                                           site_prodS = seq(from = -0.6, to = 0.6, length.out = 100)){
contour(nareaS, site_prodS, res_photo, xlab = 'Leaf N', ylab = 'Water availability')
}


plot_water_LeafN <- function(){
nareaS <-  seq(from = 1.5, to = 8, length.out = 100)/1000
lcp_vec <- lcp_narea_site_prod_height(trait_gradients_base_parameters,
                           nareaS = nareaS)
dhdt_vec_hs<- dhdt_narea_site_prod_height(trait_gradients_base_parameters,
                           nareaS = nareaS, .site_prod = -0.5)
par(mfrow = c(2, 1))
plot(nareaS, 1 - lcp_vec,  type = "l",
     ylab = expression(Height~growth~'in'~full~light~low~water~(m~year^{-1})),
     xlab = 'Narea', cex.lab = 1.2)
plot(nareaS, dhdt_vec_hs,  type = "l",
     xlab = 'Narea',
     ylab = '1 - light compensation point (%)', cex.lab = 1.2)
}



### TODO DO PLOT WITH VPD

## p <- FF16FvCB_Parameters(patch_area=1.0, hyperpar = make_FF16FvCB_hyperpar(vpd = 3))
## sF1b <- strategy(trait_matrix(c(1e-3, 4),c("narea", "vpd")), p)
## plF1b <- FF16FvCB_PlantPlus(sF1b)
## yF01b <- setNames(plF1b$ode_state, plF1b$ode_names)
## yyF1b <- deSolve::lsoda(yF01b, tt, derivs, plF1b, env=env)
## lines(height ~ time, yyF1b, type="l", lty = 1, col = "blue")

