## FF16FvCB only

# Growth over time
FF16FvCB_strategy_growth<- function(light = 1, ymax = 19){
  require(plant)
  derivs <- function(t, y, plant, env) {
   plant$ode_state <- y
   plant$compute_vars_phys(env)
   list(plant$ode_rates)
  }

  tt <- seq(0, 50, length.out=101)
  env <- fixed_environment(light)

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
  plot(height ~ time, yyF1, type="l", ylim = c(0, ymax),
       main = paste("FvCB growth for three light levels", light))
  lines(height ~ time, yyF1, type="l", lty = 1)
  lines(height ~ time, yyF2, type="l", lty = 2)
  lines(height ~ time, yyF3, type="l", lty = 3)
  legend("bottomright",
         c("Narea = 1e-3", "Narea = 3e-3", "Narea = 9e-3"),
         lty=c(1,2,3), col=c(rep("black", 3)), bty="n")
}

# Growth light

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


FF16FvCB_strategy_light_growth<- function(){
  require(plant)
  openness <- seq(0, 1, length.out=101)
  ## FF16FvCB
  pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))

  sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
  sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
  sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
  plF1 <- FF16FvCB_PlantPlus(sF1)
  plF2 <- FF16FvCB_PlantPlus(sF2)
  plF3 <- FF16FvCB_PlantPlus(sF3)

  plot(openness, sapply(openness, f_g, plF1), type="l", xlim=c(0, 1),
       las=1, xlab="Canopy openness", ylab="Height growth rate (m / yr)",
       ylim = c(0, 1.4), lwd = 2)
  lines(openness, sapply(openness, f_g, plF2), lty = 2, lwd = 2)
  lines(openness, sapply(openness, f_g, plF3), lty = 3, lwd = 2)

  pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))

  sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
  sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
  sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
  plF1 <- FF16FvCB_PlantPlus(sF1)
  plF2 <- FF16FvCB_PlantPlus(sF2)
  plF3 <- FF16FvCB_PlantPlus(sF3)

  lines(openness, sapply(openness, f_g, plF1), lty = 1)
  lines(openness, sapply(openness, f_g, plF2), lty = 2)
  lines(openness, sapply(openness, f_g, plF3), lty = 3)

  legend("topleft",
         c("Narea = 1e-3", "Narea = 3e-3", "Narea = 9e-3",
           "Low water availability", "High water availability"),
         lty=c(1,2,3,1,1), lwd = c(1,1,1,1,2),
         col=c(rep("black", 3), "black", "black"), bty="n")
}


# Growth Narea

FF16FvCB_strategy_growth_narea<- function(n_length = 100){
require(plant)


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
v_lg_F <- t(sapply(v_narea, fl2b, lights = c(0.25, 1.0)))
v_hg_F <- t(sapply(v_narea, fh2b, lights = c(0.25, 1.0)))

plot(v_narea, v_hg_F[, 2], type="l",
     ylim = range(v_hg, v_hg_F),
     las=1, xlab="Nitrogen per area",
     ylab="Height growth rate", log = "x")
lines(v_narea, v_hg_F[, 1], col = "black", lty = 2)
lines(v_narea, v_lg_F[, 2], col = "grey")
lines(v_narea, v_lg_F[, 1], col = "grey", lty = 2)

legend("topleft",
       c("High light", "Low light",
         "FF16FvCB vpd = 0",
         "FF16FvCB vpd = 3"),
       lty=c(1,2,1,1),
       col=c("black", "black", "black", "grey"), bty="n")

}


# Growth Narea vs VPD
## generate matrix of growth per narea vs vpd
FF16FvCB_strategy_growth_narea_vpd<- function(n_length = 40){
require(plant)


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
v_vpd2 <- seq(0, 3, length.out = 10)

m_growth_F<- matrix(NA, nrow = 10, ncol = n_length)

for (i in 1:10) {
m_growth_F[i, ]<- sapply(v_narea, f2b, vpd = v_vpd2[i])
}
return(m_growth_F)
}

plot_strategy_growth_narea_vpd_FF16FvCB<- function(m_g){
par(mfrow = c(1, 2))
image(m_g, xlab = "water stress", ylab = "Narea")
contour(m_g, add =TRUE)
}




# Growth lma
FF16FvCB_strategy_growth_lma<- function(n_length = 100){
require(plant)
fl2<- function(lma, lights){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))
  s <- strategy(trait_matrix(lma,"lma"), p)
  pl <- FF16FvCB_PlantPlus(s)
  pl$height <- 0.5
  sapply(lights, f_g, pl)
}

fh2<- function(lma, lights){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))
  s <- strategy(trait_matrix(lma,"lma"), p)
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
v_lma <- seq(0.4,1.1, length.out = n_length)
v_lg_F <- t(sapply(v_lma, fl2b, lights = c(0.25, 1.0)))
v_hg_F <- t(sapply(v_lma, fh2b, lights = c(0.25, 1.0)))

plot(v_lma, v_hg_F[, 2], type="l",
     ylim = range(v_lg_F, v_hg_F),
     las=1, xlab="Leaf mass per area",
     ylab="Height growth rate", log = "x")
lines(v_lma, v_hg_F[, 1], col = "black", lty = 2)
lines(v_lma, v_lg_F[, 2], col = "grey")
lines(v_lma, v_lg_F[, 1], col = "grey", lty = 2)

legend("topleft",
       c("High light", "Low light",
         "FF16FvCB vpd = 0",
         "FF16FvCB vpd = 3"),
       lty=c(1,2,1,1),
       col=c("black", "black", "black", "grey"), bty="n")

}


# LCP

# much less accurate than what is normaly done in Plant but just while lcp_whole_plant doesn't work ...
lcp_georges <- function(pl){
  openness <- seq(0, 1, length.out=501)
  v <- sapply(openness, f_g, pl)
  lcp <- mean(openness[c(max(which(v < 1e-10)), min(which(v > 1e-10)))])
  return(lcp)
}


FF16FvCB_strategy_lcp<- function(){
require(plant)

openness <- seq(0, 1, length.out=501)

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
plot(x, sapply(x, f_g, plF1), type="l", xlim=c(0, 1), ylim = c(0, 1.5),
     las=1, xlab="Canopy openness", ylab="Height growth rate (m / yr)", lwd = 2)
points(lcpF1, 0.0, pch=19, col = "black")
x <- c(lcpF2, openness[openness > lcpF2])
lines(x, sapply(x, f_g, plF2), lty = 2, col = 'black', lwd = 2)
points(lcpF2, 0.0, pch=19, col = "black")
x <- c(lcpF3, openness[openness > lcpF3])
lines(x, sapply(x, f_g, plF3), lty = 3, col = 'black', lwd = 2)
points(lcpF3, 0.0, pch=19, col = "black")

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
lines(x, sapply(x, f_g, plF1), col = 'grey')
points(lcpF1, 0.0, pch=19, col = "grey")
x <- c(lcpF2, openness[openness > lcpF2])
lines(x, sapply(x, f_g, plF2), lty = 2, col = 'grey')
points(lcpF2, 0.0, pch=19, col = "grey")
x <- c(lcpF3, openness[openness > lcpF3])
lines(x, sapply(x, f_g, plF3), lty = 3, col = 'grey')
points(lcpF3, 0.0, pch=19, col = "grey")

legend("topleft",
       c("Narea = 1e-3", "Narea = 3e-3", "Narea = 9e-3",
         "FF16FvCB vpd = 0",
         "FF16FvCB vpd = 3"),
       lty=c(1,2, 3,1,1),
       col=c("black", "black", "black","black", "grey"), bty="n")
}


## lcp in function of Narea
FF16FvCB_strategy_lcp_narea<- function(n_length = 100){
require(plant)

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
v_lcp_l_F <- sapply(v_narea, fl2b)
v_lcp_h_F <- sapply(v_narea, fh2b)

plot(v_narea, v_lcp_h_F, type="l", ylim = c(0, 1),
     las=1, xlab="Nitrogen per area", ylab="WP LCP",
     log = "x")
lines(v_narea, v_lcp_l_F, col = "grey")
legend("topleft",
       c("FF16FvCB vpd = 0",
         "FF16FvCB vpd = 3"),
       lty=c(1,1),
       col=c("black", "grey"), bty="n")
}


### IMAGE PLOT

## lcp in function of Narea
FF16FvCB_strategy_lcp_narea_vpd<- function(n_length = 40){
require(plant)

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
v_vpd2 <- seq(0, 3, length.out = n_length)

m_lcp_F<- matrix(NA, nrow = n_length, ncol = n_length)

for (i in 1:n_length) {
m_lcp_F[i, ]<- sapply(v_narea, f2b, vpd = v_vpd2[i])
}
return( m_lcp_F)
}


plot_strategy_lcp_narea_vpd_FF16FvCB<- function(l_data){
m_lcp_F <- l_data
image(m_lcp_F, xlab = "water stress", ylab = "Narea")
contour(m_lcp_F, add =TRUE)
}


#########################
# TRADEOFF FF16FvCB

FF16FvCB_lcp_lma_vpd_height <- function(lmaS, .vpd = 0,
                                        height = 0.3920458, narea=1.87e-3){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = .vpd))
  lcp_vec <- rep(NA, length(lmaS))
  f_lcp_lma <- function(lma, narea, p, height){
      s <- strategy(trait_matrix(c(lma, narea), c("lma", "narea")), p)
      pF<- FF16FvCB_PlantPlus(s)
      pF$height <- height
      lcp_whole_plant(pF)
  }
 lcp_vec <-  sapply(lmaS, f_lcp_lma, narea, p, height)
 return(lcp_vec)
}

FF16FvCB_dhdt_lma_vpd_height<- function(lmaS, .vpd = 0,
                                     height = 0.3920458, narea=1.87e-3){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = .vpd))
  dhdt_vec <- rep(NA, length(lmaS))
  f_g_lma <- function(lma, narea, p, height){
      s <- strategy(trait_matrix(c(lma, narea), c("lma", "narea")), p)
      pF<- FF16FvCB_PlantPlus(s)
      pF$height <- height
      f_g(1, pF)
  }
 dhdt_vec <-  sapply(lmaS, f_g_lma, narea, p, height)
 return(dhdt_vec)
}

FF16FvCB_lcp_narea_vpd_height <- function(nareaS, .vpd = 0,
                                    height = 0.3920458, lma= 0.1978791){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = .vpd))
  lcp_vec <- rep(NA, length(nareaS))
  f_lcp_narea <- function(narea, lma, p, height){
      s <- strategy(trait_matrix(c(lma, narea), c("lma", "narea")), p)
      pF<- FF16FvCB_PlantPlus(s)
      pF$height <- height
      lcp_whole_plant(pF)
  }
 lcp_vec <-  sapply(nareaS, f_lcp_narea, lma, p, height)
 return(lcp_vec)
}

FF16FvCB_dhdt_narea_vpd_height<- function(nareaS, .vpd = 0,
                                     height = 0.3920458, lma=0.1978791){
  p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = .vpd))
  dhdt_vec <- rep(NA, length(nareaS))
  f_g_narea <- function(narea, lma, p, height){
      s <- strategy(trait_matrix(c(lma, narea), c("lma", "narea")), p)
      pF<- FF16FvCB_PlantPlus(s)
      pF$height <- height
      f_g(1, pF)
  }
 dhdt_vec <-  sapply(nareaS, f_g_narea, lma, p, height)
 return(dhdt_vec)
}


plot_lma_shade_tradeoff_FF16FvCB <- function(lmaS = seq(from = 0.06,
                                                        to = 0.27,
                                                        length.out = 100),
                                             height = 1, narea=1.87e-3){
  dhdt_vec_0 <- FF16FvCB_dhdt_lma_vpd_height(lmaS, .vpd = 0,
                                        height, narea)
  lcp_vec_0 <- FF16FvCB_lcp_lma_vpd_height(lmaS, .vpd = 0,
                                      height, narea)

  dhdt_vec_2 <- FF16FvCB_dhdt_lma_vpd_height(lmaS, .vpd = 2,
                                        height, narea)
  lcp_vec_2 <- FF16FvCB_lcp_lma_vpd_height(lmaS, .vpd = 2,
                                      height, narea)
par(mfcol = c(2,2))
plot(lmaS, 1 - lcp_vec_0, type = 'l',
     xlab = 'LMA',
     ylab = '1 - light compensation point (%)',
     ylim = range(1 - lcp_vec_0, 1 - lcp_vec_2))
lines(lmaS, 1 - lcp_vec_2, lty=2)
plot(lmaS, dhdt_vec_0, type = 'l',
     ylab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
     xlab = 'LMA',
     ylim = range(dhdt_vec_0, dhdt_vec_2))
lines(lmaS, dhdt_vec_2, lty=2)
plot(dhdt_vec_0, 1 - lcp_vec_0, type = 'l',
     xlab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
     ylab = '1 - light compensation point (%)',
     xlim = range(dhdt_vec_0, dhdt_vec_2),
     ylim = range(1 - lcp_vec_0, 1 - lcp_vec_2))
lines(dhdt_vec_2, 1 - lcp_vec_2, lty=2)
}


plot_narea_shade_tradeoff_FF16FvCB <- function(nareaS = seq(3e-4,3e-2,
                                                            length.out = 100),
                                             height = 1, lma=0.1978791){
  dhdt_vec_0 <- FF16FvCB_dhdt_narea_vpd_height(nareaS, .vpd = 0,
                                        height, lma)
  lcp_vec_0 <- FF16FvCB_lcp_narea_vpd_height(nareaS, .vpd = 0,
                                      height, lma)

  dhdt_vec_2 <- FF16FvCB_dhdt_narea_vpd_height(nareaS, .vpd = 2,
                                        height, lma)
  lcp_vec_2 <- FF16FvCB_lcp_narea_vpd_height(nareaS, .vpd = 2,
                                             height, lma)

   dhdt_vec_0[!is.finite( dhdt_vec_0)] <- NA
   dhdt_vec_2[!is.finite( dhdt_vec_2)] <- NA
   lcp_vec_0[!is.finite( lcp_vec_0)] <- NA
   lcp_vec_2[!is.finite( lcp_vec_2)] <- NA

par(mfcol = c(2,2))
plot(nareaS, 1 - lcp_vec_0, type = 'l',
     xlab = 'Narea',
     ylab = '1 - light compensation point (%)',
     ylim = range(1 - lcp_vec_0, 1 - lcp_vec_2, na.rm = TRUE))
lines(nareaS, 1 - lcp_vec_2, lty=2)
plot(nareaS, dhdt_vec_0, type = 'l',
     ylab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
     xlab = 'Narea',
     ylim = range(dhdt_vec_0, dhdt_vec_2, na.rm = TRUE))
lines(nareaS, dhdt_vec_2, lty=2)
plot(dhdt_vec_0, 1 - lcp_vec_0, type = 'l',
     xlab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
     ylab = '1 - light compensation point (%)',
     xlim = range(dhdt_vec_0, dhdt_vec_2, na.rm = TRUE),
     ylim = range(1 - lcp_vec_0, 1 - lcp_vec_2, na.rm = TRUE))
lines(dhdt_vec_2, 1 - lcp_vec_2, lty=2)
}


### PLOT DYNAMICS FOR DIFFERENT LMA AND Narea

## look at patch dynamcis for two different N value


data_patch_FF16FvCB<- function(vpd = 0, cancel_F = FALSE, v_narea = c(1e-3, 3e-3, 5e-3)){
ctrl <- equilibrium_verbose(fast_control())
ctrl$schedule_eps <- 0.002
ctrl$equilibrium_eps <- 1e-4
ctrl$equilibrium_solver_name <- "iteration"
ctrl$equilibrium_verbose <-  TRUE
pF0 <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = vpd))
pF0$control <-  ctrl
pF0$disturbance_mean_interval <- 30.0
if(cancel_F){
pF0$strategy_default$a_f1 <- 0.5
pF0$strategy_default$a_f2 <- 0
}
print(vpd)
## First, with three species:
pF <- expand_parameters(trait_matrix(v_narea, "narea"), pF0, FALSE)
pF_eq <- equilibrium_seed_rain(pF)
## Then collect the patch-level dynamics:
dataF <- run_scm_collect(pF_eq)
return(dataF)
}


plot_patch_data <- function(data, v_narea = c(1e-3, 3e-3, 5e-3),
                              title = "FF16 with prod = 0.3",
                              legend_TF = TRUE){
rel <- function(x, xmin) {
  x[x < xmin] <- xmin
  xmax <- max(x, na.rm=TRUE)
  (x - xmin) / (xmax - xmin)
}
cols <- c("#e41a1c", "#377eb8", "#4daf4a")
tt <- data$time
h1 <- data$species[[1]]["height", , ] # TODO

mh1 <- apply(h1, MARGIN = 1, mean , na.rm = TRUE)
h2 <- data$species[[2]]["height", , ]
mh2 <- apply(h2, MARGIN = 1, mean , na.rm = TRUE)
h3 <- data$species[[3]]["height", , ]
mh3 <- apply(h3, MARGIN = 1, mean , na.rm = TRUE)
## Relativise the log densities onto (-4, max)
d1 <- data$species[[1]]["log_density", , ]
d1[d1 < -5] <- -5
md1 <- apply(exp(d1), MARGIN = 1, sum, na.rm = TRUE)
d2 <- data$species[[2]]["log_density", , ]
d2[d2 < -5] <- -5
md2 <- apply(exp(d2), MARGIN = 1, sum, na.rm = TRUE)
d3 <- data$species[[3]]["log_density", , ]
d3[d3 < -5] <- -5
md3 <- apply(exp(d3), MARGIN = 1, sum, na.rm = TRUE)
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
     ylim = range(h1,h2,h3, na.rm = TRUE), main = title)
segments(x1[-1, ], h1[-1, ], x1[-n, ], h1[-n, ], col=col1[-n, ], lend="butt")
segments(x2[-1, ], h2[-1, ], x2[-n, ], h2[-n, ], col=col2[-n, ], lend="butt")
segments(x3[-1, ], h3[-1, ], x3[-n, ], h3[-n, ], col=col3[-n, ], lend="butt")

# size structure
metapopulation <- function(x){
  plant:::trapezium(t2, x*data2$patch_density)
}

patches <- lapply(seq_along(data$time), scm_patch, data)
lai <- sapply(patches, function(x) x$area_leaf_above(0.0))
lai_1 <- sapply(patches, function(x) x$species[[1]]$area_leaf_above(0.0))
lai_2 <- sapply(patches, function(x) x$species[[2]]$area_leaf_above(0.0))
lai_3 <- sapply(patches, function(x) x$species[[3]]$area_leaf_above(0.0))

ff <- function(x, i){
dd <- x$species[[i]]$log_densities
dd[dd< -5] <- -5
sum(exp(dd), na.rm = TRUE)
}

plot(tt, lai, type="l", las=1, lty=3, lwd = 3,
     xlab="Time (years)", ylab="Leaf area index")
lines(tt, lai_1, col=cols[[1]], lwd = 2)
lines(tt, lai_2, col=cols[[2]], lwd = 2)
lines(tt, lai_3, col=cols[[3]], lwd = 2)

metapopulation <- function(x){
  plant:::trapezium(tt, x*data$patch_density)
}
lai_av <- metapopulation(lai)
lai_1_av <- metapopulation(lai_1)
lai_2_av <- metapopulation(lai_2)
lai_3_av <- metapopulation(lai_3)
axis(4, at=lai_av,   tck=0.1, lty = 2, labels=NA, lwd = 5)
axis(4, at=lai_1_av, tck=0.1, col.ticks=cols[[1]], labels=NA, lwd = 4)
axis(4, at=lai_2_av, tck=0.1, col.ticks=cols[[2]], labels=NA, lwd = 4)
axis(4, at=lai_3_av, tck=0.1, col.ticks=cols[[3]], labels=NA,lwd = 4)
axis(1, at=108, labels="Av")

pd_1 <- sapply(patches, ff, i = 1)
pd_2 <- sapply(patches, ff, i = 2)
pd_3 <- sapply(patches, ff, i = 3)
plot(tt, pd_1, type="l", lty=1, lwd = 2,
     xlab="Time (years)", ylab="Density (1 / m / m2)",
     ylim = range(pd_1, pd_2, pd_3, na.rm = TRUE), col = cols[[1]], log = "y")
lines(tt, pd_2, col=cols[[2]], lwd = 2)
lines(tt, pd_3, col=cols[[3]], lwd = 2)
pd_1_av <- metapopulation(pd_1)
pd_2_av <- metapopulation(pd_2)
pd_3_av <- metapopulation(pd_3)
axis(4, at=pd_1_av, tck=0.1, col.ticks=cols[[1]], labels=NA, lwd = 4)
axis(4, at=pd_2_av, tck=0.1, col.ticks=cols[[2]], labels=NA, lwd = 4)
axis(4, at=pd_3_av, tck=0.1, col.ticks=cols[[3]], labels=NA,lwd = 4)
axis(1, at=108, labels="Av")

if(legend_TF){
 legend("topright", legend = paste0("Narea = ", v_narea), lty=1, col = cols, bty="n")
 }
}

plot_patch_data_grad<- function(data_l, data_h,
                                  title = "FF16",
                                  v_narea = c(1e-3, 3e-3, 5e-3)){

par(mfrow = c(2,3))
plot_patch_data(data_l, title = paste(title, "at low gradient"),
                legend_TF = FALSE)
plot_patch_data(data_h, title = paste(title, "at high gradient"),
                legend_TF = TRUE)
}




## ## Compare FF16 and FF16FvCB

## ## Compare FF16 and FF16FvCB
## compare_strategy_fecundity<- function(light = 1, ymax = 19){
##     require(plant)
## derivs <- function(t, y, plant, env) {
##  plant$ode_state <- y
##  plant$compute_vars_phys(env)
##  list(plant$ode_rates)
## }

## tt <- seq(0, 100, length.out=101)
## env <- fixed_environment(light)

## p <- scm_base_parameters("FF16")

## cols <- c("#e41a1c", "#377eb8", "#4daf4a")

## s1 <- strategy(trait_matrix(1e-3,"narea"), p)
## s2 <- strategy(trait_matrix(2e-3,"narea"), p)
## s3 <- strategy(trait_matrix(5e-3,"narea"), p)
## pl1 <- FF16_PlantPlus(s1)
## pl2 <- FF16_PlantPlus(s2)
## pl3 <- FF16_PlantPlus(s3)

## y01 <- setNames(pl1$ode_state, pl1$ode_names)
## yy1 <- deSolve::lsoda(y01, tt, derivs, pl1, env=env)
## y02 <- setNames(pl2$ode_state, pl2$ode_names)
## yy2 <- deSolve::lsoda(y02, tt, derivs, pl2, env=env)
## y03 <- setNames(pl3$ode_state, pl3$ode_names)
## yy3 <- deSolve::lsoda(y03, tt, derivs, pl3, env=env)
## par(mfrow= c(2,3))
## plot(fecundity ~ time, yy1, type="l",
##      ylim = c(0, max( yy1[, "fecundity"],  yy2[, "fecundity"],  yy3[, "fecundity"])),
##      main = "FF16 fecundity function Hmat", col = cols[[1]], lwd = 2)
## lines(yy2[, "time"], yy2[, "fecundity"], lty = 2, col = cols[[2]], lwd = 2)
## lines(yy3[, "time"], yy3[, "fecundity"], lty = 3, col = cols[[3]], lwd = 2)

## plot(fecundity ~ height, yy1, type="l",
##      ylim = c(0, max( yy1[, "fecundity"],  yy2[, "fecundity"],  yy3[, "fecundity"])),
##      xlim = c(0, 19),
##      main = "FF16 fecundity function Hmat", col = cols[[1]], lwd = 2)
## lines(yy2[, "height"], yy2[, "fecundity"], lty = 2, col = cols[[2]], lwd = 2)
## lines(yy3[, "height"], yy3[, "fecundity"], lty = 3, col = cols[[3]], lwd = 2)

## plot(height ~ time, yy1, type="l",
##      ylim = c(0, max( yy1[, "height"],  yy2[, "height"],  yy3[, "height"])),
##      main = "FF16 fecundity function Hmat", col = cols[[1]], lwd = 2)
## lines(yy2[, "time"], yy2[, "height"], lty = 2, col = cols[[2]], lwd = 2)
## lines(yy3[, "time"], yy3[, "height"], lty = 3, col = cols[[3]], lwd = 2)

## p <- scm_base_parameters("FF16")
## p$strategy_default$a_f1 <- 0.5
## p$strategy_default$a_f2 <- 0

## s1 <- strategy(trait_matrix(1e-3,"narea"), p)
## s2 <- strategy(trait_matrix(3e-3,"narea"), p)
## s3 <- strategy(trait_matrix(9e-3,"narea"), p)
## pl1 <- FF16_PlantPlus(s1)
## pl2 <- FF16_PlantPlus(s2)
## pl3 <- FF16_PlantPlus(s3)

## y01 <- setNames(pl1$ode_state, pl1$ode_names)
## yy1 <- deSolve::lsoda(y01, tt, derivs, pl1, env=env)
## y02 <- setNames(pl2$ode_state, pl2$ode_names)
## yy2 <- deSolve::lsoda(y02, tt, derivs, pl2, env=env)
## y03 <- setNames(pl3$ode_state, pl3$ode_names)
## yy3 <- deSolve::lsoda(y03, tt, derivs, pl3, env=env)
## plot(fecundity ~ time, yy1, type="l",
##      ylim = c(0, max( yy1[, "fecundity"],  yy2[, "fecundity"],  yy3[, "fecundity"])),
##      main = "FF16 fecundity constant", lwd = 2)
## lines(yy2[, "time"], yy2[, "fecundity"], lty = 2, lwd = 2)
## lines(yy3[, "time"], yy3[, "fecundity"], lty = 3, lwd = 2)
## legend("topleft",
##        c("Narea = 1e-3", "Narea = 2e-3", "Narea = 5e-3"),
##        lty=c(1,2,3), lwd = 2, col=cols, bty="n")

## plot(fecundity ~ height, yy1, type="l",
##      ylim = c(0, max( yy1[, "fecundity"],  yy2[, "fecundity"],  yy3[, "fecundity"])),
##      xlim = c(0, 60),
##      main = "FF16 fecundity constant", col = cols[[1]], lwd = 2)
## lines(yy2[, "height"], yy2[, "fecundity"], lty = 2, col = cols[[2]], lwd = 2)
## lines(yy3[, "height"], yy3[, "fecundity"], lty = 3, col = cols[[3]], lwd = 2)

## plot(height ~ time, yy1, type="l",
##      ylim = c(0, max( yy1[, "height"],  yy2[, "height"],  yy3[, "height"])),
##      main = "FF16 fecundity constant", col = cols[[1]], lwd = 2)
## lines(yy2[, "time"], yy2[, "height"], lty = 2, col = cols[[2]], lwd = 2)
## lines(yy3[, "time"], yy3[, "height"], lty = 3, col = cols[[3]], lwd = 2)
## }

## compare_strategy_growth<- function(light = 1, ymax = 19){
## require(plant)
## derivs <- function(t, y, plant, env) {
##  plant$ode_state <- y
##  plant$compute_vars_phys(env)
##  list(plant$ode_rates)
## }

## tt <- seq(0, 50, length.out=101)
## env <- fixed_environment(light)

## p <- scm_base_parameters("FF16")
## s1 <- strategy(trait_matrix(1e-3,"narea"), p)
## s2 <- strategy(trait_matrix(3e-3,"narea"), p)
## s3 <- strategy(trait_matrix(9e-3,"narea"), p)
## pl1 <- FF16_PlantPlus(s1)
## pl2 <- FF16_PlantPlus(s2)
## pl3 <- FF16_PlantPlus(s3)

## y01 <- setNames(pl1$ode_state, pl1$ode_names)
## yy1 <- deSolve::lsoda(y01, tt, derivs, pl1, env=env)
## y02 <- setNames(pl2$ode_state, pl2$ode_names)
## yy2 <- deSolve::lsoda(y02, tt, derivs, pl2, env=env)
## y03 <- setNames(pl3$ode_state, pl3$ode_names)
## yy3 <- deSolve::lsoda(y03, tt, derivs, pl3, env=env)
## plot(height ~ time, yy1, type="l", ylim = c(0, ymax), main = paste("light level", light))
## lines(yy2[, "time"], yy2[, "height"], lty = 2)
## lines(yy3[, "time"], yy3[, "height"], lty = 3)

## pF <- scm_base_parameters("FF16FvCB")

## sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
## sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
## sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
## plF1 <- FF16FvCB_PlantPlus(sF1)
## plF2 <- FF16FvCB_PlantPlus(sF2)
## plF3 <- FF16FvCB_PlantPlus(sF3)

## yF01 <- setNames(plF1$ode_state, plF1$ode_names)
## yyF1 <- deSolve::lsoda(yF01, tt, derivs, plF1, env=env)
## yF02 <- setNames(plF2$ode_state, plF2$ode_names)
## yyF2 <- deSolve::lsoda(yF02, tt, derivs, plF2, env=env)
## yF03 <- setNames(plF3$ode_state, plF3$ode_names)
## yyF3 <- deSolve::lsoda(yF03, tt, derivs, plF3, env=env)
## lines(height ~ time, yyF1, type="l", lty = 1, col = "red")
## lines(height ~ time, yyF2, type="l", lty = 2, col = "red")
## lines(height ~ time, yyF3, type="l", lty = 3, col = "red")
## legend("bottomright",
##        c("Narea = 1e-3", "Narea = 3e-3", "Narea = 9e-3", "FF16", "FF16FvCB"),
##        lty=c(1,2,3,1,1), col=c(rep("black", 3), "red", "black"), bty="n")
## }


## f_g<- function(x, pl) {
##   env <- fixed_environment(x)
##   pl$compute_vars_phys(env)
##   pl$internals$height_dt
## }

## f_demo_size<- function(x, pl, sizeS) {
##   env <- fixed_environment(x)
##   pl$compute_vars_phys(env)
##   pl$internals$height_dt
## }


## compare_strategy_light_growth<- function(){
## require(plant)
## openness <- seq(0, 1, length.out=101)

## p <- trait_gradients_base_parameters(site_prod=0.3)

## s1 <- strategy(trait_matrix(1e-3,"narea"), p)
## s2 <- strategy(trait_matrix(3e-3,"narea"), p)
## s3 <- strategy(trait_matrix(9e-3,"narea"), p)
## pl1 <- FF16_PlantPlus(s1)
## pl2 <- FF16_PlantPlus(s2)
## pl3 <- FF16_PlantPlus(s3)

## plot(openness, sapply(openness, f_g, pl1), type="l", xlim=c(0, 1),
##      las=1, xlab="Canopy openness", ylab="Height growth rate (m / yr)",
##      ylim = c(0, 1.4), lwd = 2)
## lines(openness, sapply(openness, f_g, pl2), lty = 2, lwd = 2)
## lines(openness, sapply(openness, f_g, pl3), lty = 3, lwd = 2)

## p <- trait_gradients_base_parameters(site_prod=-0.3)

## s1 <- strategy(trait_matrix(1e-3,"narea"), p)
## s2 <- strategy(trait_matrix(3e-3,"narea"), p)
## s3 <- strategy(trait_matrix(9e-3,"narea"), p)
## pl1 <- FF16_PlantPlus(s1)
## pl2 <- FF16_PlantPlus(s2)
## pl3 <- FF16_PlantPlus(s3)
## lines(openness, sapply(openness, f_g, pl1), lty = 1)
## lines(openness, sapply(openness, f_g, pl2), lty = 2)
## lines(openness, sapply(openness, f_g, pl3), lty = 3)

## ## FF16FvCB
## pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))

## sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
## sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
## sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
## plF1 <- FF16FvCB_PlantPlus(sF1)
## plF2 <- FF16FvCB_PlantPlus(sF2)
## plF3 <- FF16FvCB_PlantPlus(sF3)

## lines(openness, sapply(openness, f_g, plF1), lty = 1, col = "red", lwd = 2)
## lines(openness, sapply(openness, f_g, plF2), lty = 2, col = "red", lwd = 2)
## lines(openness, sapply(openness, f_g, plF3), lty = 3, col = "red", lwd = 2)

## pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))

## sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
## sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
## sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
## plF1 <- FF16FvCB_PlantPlus(sF1)
## plF2 <- FF16FvCB_PlantPlus(sF2)
## plF3 <- FF16FvCB_PlantPlus(sF3)

## lines(openness, sapply(openness, f_g, plF1), lty = 1, col = "red")
## lines(openness, sapply(openness, f_g, plF2), lty = 2, col = "red")
## lines(openness, sapply(openness, f_g, plF3), lty = 3, col = "red")

## legend("topleft",
##        c("Narea = 1e-3", "Narea = 3e-3", "Narea = 9e-3", "FF16", "FF16FvCB",
##          "Low water availability", "High water availability"),
##        lty=c(1,2,3,1,1, 1,1), lwd = c(1,1,1,1,1,1,2),
##        col=c(rep("black", 3), "red", "black", "black", "black"), bty="n")
## }


## compare_strategy_growth_narea<- function(n_length = 100){
## require(plant)

## fl1 <- function(narea, lights){
##   p <- trait_gradients_base_parameters(site_prod=-0.3)
##   s <- strategy(trait_matrix(narea,"narea"), p)
##   pl <- FF16_PlantPlus(s)
##   pl$height <- 0.5
##   sapply(lights, f_g, pl)
## }

## fh1 <- function(narea, lights){
##   p <- trait_gradients_base_parameters(site_prod=0.3)
##   s <- strategy(trait_matrix(narea,"narea"), p)
##   pl <- FF16_PlantPlus(s)
##   pl$height <- 0.5
##   sapply(lights, f_g, pl)
## }

## fl2<- function(narea, lights){
##   p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))
##   s <- strategy(trait_matrix(narea,"narea"), p)
##   pl <- FF16FvCB_PlantPlus(s)
##   pl$height <- 0.5
##   sapply(lights, f_g, pl)
## }

## fh2<- function(narea, lights){
##   p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))
##   s <- strategy(trait_matrix(narea,"narea"), p)
##   pl <- FF16FvCB_PlantPlus(s)
##   pl$height <- 0.5
##   sapply(lights, f_g, pl)
## }

## fl2b <- function(n, lights){
## tryCatch(fl2(n, lights), error=function(err) NA)
## }

## fh2b <- function(n, lights){
## tryCatch(fh2(n, lights), error=function(err) NA)
## }

## v_narea1 <- 10^seq(log10(3e-4),log10(3e-2), length.out = n_length/2)
## v_narea2 <- seq(3e-4,3e-2, length.out = n_length/2)
## v_narea <-  sort(c(v_narea1, v_narea2))
## v_lg <- t(sapply(v_narea, fl1, lights = c(0.25, 1.0)))
## v_hg <- t(sapply(v_narea, fh1, lights = c(0.25, 1.0)))
## v_lg_F <- t(sapply(v_narea, fl2b, lights = c(0.25, 1.0)))
## v_hg_F <- t(sapply(v_narea, fh2b, lights = c(0.25, 1.0)))

## plot(v_narea, v_hg[, 2], type="l", ylim = range(v_lg, v_lg_F, v_hg, v_hg_F),
##      las=1, xlab="Nitrogen per area", ylab="Height growth rate", log = "x")
## lines(v_narea, v_hg[, 1], lty = 2)
## lines(v_narea, v_lg[, 1], lty = 2, col = "grey")
## lines(v_narea, v_lg[, 2], lty = 1, col = "grey")
## lines(v_narea, v_hg_F[, 2], col = "red")
## lines(v_narea, v_hg_F[, 1], col = "red", lty = 2)
## lines(v_narea, v_lg_F[, 2], col = "pink")
## lines(v_narea, v_lg_F[, 1], col = "pink", lty = 2)

## legend("topleft",
##        c("High light", "Low light", "FF16 site_prod = 0.3", "FF16 site_prod = -0.3", "FF16FvCB vpd = 0",
##          "FF16FvCB vpd = 3"),
##        lty=c(1,2,1,1, 1,1),
##        col=c("black", "black", "black", "grey", "red", "pink"), bty="n")

## }


## ## growth in function of Narea
## compare_strategy_growth_narea_vpd<- function(n_length = 40){
## require(plant)

## f1 <- function(n, vpd = 0.3){
##   p <- trait_gradients_base_parameters(site_prod=vpd)
##   s <- strategy(trait_matrix(n,"narea"), p)
##   pl <- FF16_PlantPlus(s)
##   pl$height <- 0.5
##   f_g(1.0, pl)
## }


## f2<- function(n, vpd = 0){
##   p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = vpd))
##   s <- strategy(trait_matrix(n,"narea"), p)
##   pl <- FF16FvCB_PlantPlus(s)
##   pl$height <- 0.5
##   f_g(1.0, pl)
## }

## f2b <- function(n, vpd){
## tryCatch(f2(n, vpd), error=function(err) NA)
## }

## v_narea <- seq(1e-3,3e-2, length.out = n_length)
## v_vpd1 <- seq(-0.3, 0.3, length.out = 10)
## v_vpd2 <- seq(0, 3, length.out = 10)

## m_growth <- matrix(NA, nrow = 10, ncol = n_length)
## m_growth_F<- matrix(NA, nrow = 10, ncol = n_length)

## for (i in 1:10) {
## m_growth[i, ]<- sapply(v_narea, f1, vpd = v_vpd1[i])
## m_growth_F[i, ]<- sapply(v_narea, f2b, vpd = v_vpd2[i])
## }
## return(list(m_growth, m_growth_F))
## }

## plot_strategy_growth_narea_vpd<- function(l_data){
## m_lcp <-  l_data[[1]]
## m_lcp_F <- l_data[[2]]
## par(mfrow = c(1, 2))
## image(m_lcp[rev(seq_len(nrow(m_lcp))), ], xlab = "water stress", ylab = "Narea")
## contour(m_lcp[rev(seq_len(nrow(m_lcp))), ], add =TRUE)
## image(m_lcp_F, xlab = "water stress", ylab = "Narea")
## contour(m_lcp_F, add =TRUE)
## }


## # much less accurate than what is normaly done in Plant but just while lcp_whole_plant doesn't work ...
## lcp_georges <- function(pl){
##   openness <- seq(0, 1, length.out=501)
##   v <- sapply(openness, f_g, pl)
##   lcp <- mean(openness[c(max(which(v < 1e-10)), min(which(v > 1e-10)))])
##   return(lcp)
## }


## compare_strategy_lcp<- function(){
## require(plant)
## p <- trait_gradients_base_parameters(site_prod=0.3)
## s1 <- strategy(trait_matrix(1e-3,"narea"), p)
## s2 <- strategy(trait_matrix(3e-3,"narea"), p)
## s3 <- strategy(trait_matrix(9e-3,"narea"), p)
## pl1 <- FF16_PlantPlus(s1)
## pl2 <- FF16_PlantPlus(s2)
## pl3 <- FF16_PlantPlus(s3)

## openness <- seq(0, 1, length.out=501)

## lcp1 <- lcp_georges(pl1)
## lcp2 <- lcp_georges(pl2)
## lcp3 <- lcp_georges(pl3)

## x <- c(lcp1, openness[openness > lcp1])
## plot(x, sapply(x, f_g, pl1), type="l", xlim=c(0, 1), ylim = c(0, 1.5),
##      las=1, xlab="Canopy openness", ylab="Height growth rate (m / yr)", lwd = 2)
## points(lcp1, 0.0, pch=19)
## x <- c(lcp2, openness[openness > lcp2])
## lines(x, sapply(x, f_g, pl2), lty = 2, lwd = 2)
## points(lcp2, 0.0, pch=19)
## x <- c(lcp3, openness[openness > lcp3])
## lines(x, sapply(x, f_g, pl3), lty = 3, lwd = 2)
## points(lcp3, 0.0, pch=19)


## p <- trait_gradients_base_parameters(site_prod=-0.3)
## s1 <- strategy(trait_matrix(1e-3,"narea"), p)
## s2 <- strategy(trait_matrix(3e-3,"narea"), p)
## s3 <- strategy(trait_matrix(9e-3,"narea"), p)
## pl1 <- FF16_PlantPlus(s1)
## pl2 <- FF16_PlantPlus(s2)
## pl3 <- FF16_PlantPlus(s3)

## openness <- seq(0, 1, length.out=51)

## lcp1 <- lcp_georges(pl1)
## lcp2 <- lcp_georges(pl2)
## lcp3 <- lcp_georges(pl3)

## x <- c(lcp1, openness[openness > lcp1])
## lines(x, sapply(x, f_g, pl1), lty = 1, col = "grey")
## points(lcp1, 0.0, pch=19, col = "grey")
## x <- c(lcp2, openness[openness > lcp2])
## lines(x, sapply(x, f_g, pl2), lty = 2, col = "grey")
## points(lcp2, 0.0, pch=19, col = "grey")
## x <- c(lcp3, openness[openness > lcp3])
## lines(x, sapply(x, f_g, pl3), lty = 3, col = "grey")
## points(lcp3, 0.0, pch=19, col = "grey")

## ## FF16FvCB
## pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))

## sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
## sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
## sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
## plF1 <- FF16FvCB_PlantPlus(sF1)
## plF2 <- FF16FvCB_PlantPlus(sF2)
## plF3 <- FF16FvCB_PlantPlus(sF3)

## lcpF1 <- lcp_georges(plF1)
## lcpF2 <- lcp_georges(plF2)
## lcpF3 <- lcp_georges(plF3)

## x <- c(lcpF1, openness[openness > lcpF1])
## lines(x, sapply(x, f_g, plF1), col = 'red', lwd = 2)
## points(lcpF1, 0.0, pch=19, col = "red")
## x <- c(lcpF2, openness[openness > lcpF2])
## lines(x, sapply(x, f_g, plF2), lty = 2, col = 'red', lwd = 2)
## points(lcpF2, 0.0, pch=19, col = "red")
## x <- c(lcpF3, openness[openness > lcpF3])
## lines(x, sapply(x, f_g, plF3), lty = 3, col = 'red', lwd = 2)
## points(lcpF3, 0.0, pch=19, col = "red")

## pF <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))

## sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
## sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
## sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
## plF1 <- FF16FvCB_PlantPlus(sF1)
## plF2 <- FF16FvCB_PlantPlus(sF2)
## plF3 <- FF16FvCB_PlantPlus(sF3)

## lcpF1 <- lcp_georges(plF1)
## lcpF2 <- lcp_georges(plF2)
## lcpF3 <- lcp_georges(plF3)

## x <- c(lcpF1, openness[openness > lcpF1])
## lines(x, sapply(x, f_g, plF1), col = 'pink')
## points(lcpF1, 0.0, pch=19, col = "pink")
## x <- c(lcpF2, openness[openness > lcpF2])
## lines(x, sapply(x, f_g, plF2), lty = 2, col = 'pink')
## points(lcpF2, 0.0, pch=19, col = "pink")
## x <- c(lcpF3, openness[openness > lcpF3])
## lines(x, sapply(x, f_g, plF3), lty = 3, col = 'pink')
## points(lcpF3, 0.0, pch=19, col = "pink")

## legend("topleft",
##        c("Narea = 1e-3", "Narea = 3e-3", "Narea = 9e-3",
##          "FF16 site_prod = 0.3", "FF16 site_prod = -0.3", "FF16FvCB vpd = 0",
##          "FF16FvCB vpd = 3"),
##        lty=c(1,2, 3,1,1, 1,1),
##        col=c("black", "black", "black","black", "grey", "red", "pink"), bty="n")

## }


## ## lcp in function of Narea
## compare_strategy_lcp_narea<- function(n_length = 100){
## require(plant)

## fl1 <- function(n){
##   p <- trait_gradients_base_parameters(site_prod=-0.3)
##   s <- strategy(trait_matrix(n,"narea"), p)
##   pl <- FF16_PlantPlus(s)
##   pl$height <- 0.5
##   lcp_georges(pl)
## }

## fh1 <- function(n){
##   p <- trait_gradients_base_parameters(site_prod=0.3)
##   s <- strategy(trait_matrix(n,"narea"), p)
##   pl <- FF16_PlantPlus(s)
##   pl$height <- 0.5
##   lcp_georges(pl)
## }

## fl2<- function(n){
##   p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 3))
##   s <- strategy(trait_matrix(n,"narea"), p)
##   pl <- FF16FvCB_PlantPlus(s)
##   pl$height <- 0.5
##   lcp_georges(pl)
## }

## fh2<- function(n){
##   p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = 0))
##   s <- strategy(trait_matrix(n,"narea"), p)
##   pl <- FF16FvCB_PlantPlus(s)
##   pl$height <- 0.5
##   lcp_georges(pl)
## }
## fl2b <- function(n){
## tryCatch(fl2(n), error=function(err) NA)
## }
## fh2b <- function(n){
## tryCatch(fh2(n), error=function(err) NA)
## }

## v_narea1 <- 10^seq(log10(1e-3),log10(3e-2), length.out = n_length/2)
## v_narea2 <- seq(1e-3,3e-2, length.out = n_length/2)
## v_narea <-  sort(c(v_narea1, v_narea2))
## v_lcp_l<- sapply(v_narea, fl1)
## v_lcp_h<- sapply(v_narea, fh1)
## v_lcp_l_F <- sapply(v_narea, fl2b)
## v_lcp_h_F <- sapply(v_narea, fh2b)

## plot(v_narea, v_lcp_h, type="l", ylim = c(0, 1),
##      las=1, xlab="Nitrogen per area", ylab="WP LCP",
##      log = "x")
## lines(v_narea, v_lcp_l, col = "grey")
## lines(v_narea, v_lcp_l_F, col = "pink")
## lines(v_narea, v_lcp_h_F, col = "red")
## legend("topleft",
##        c("FF16 site_prod = 0.3", "FF16 site_prod = -0.3", "FF16FvCB vpd = 0",
##          "FF16FvCB vpd = 3"),
##        lty=c(1,1, 1,1),
##        col=c("black", "grey", "red", "pink"), bty="n")
## }



## ## lcp in function of Narea
## compare_strategy_lcp_narea_vpd<- function(n_length = 40){
## require(plant)

## f1 <- function(n, vpd = 0.3){
##   p <- trait_gradients_base_parameters(site_prod=vpd)
##   s <- strategy(trait_matrix(n,"narea"), p)
##   pl <- FF16_PlantPlus(s)
##   pl$height <- 0.5
##   lcp_georges(pl)
## }


## f2<- function(n, vpd = 0){
##   p <- FF16FvCB_Parameters(hyperpar = make_FF16FvCB_hyperpar(vpd = vpd))
##   s <- strategy(trait_matrix(n,"narea"), p)
##   pl <- FF16FvCB_PlantPlus(s)
##   pl$height <- 0.5
##   lcp_georges(pl)
## }

## f2b <- function(n, vpd){
## tryCatch(f2(n, vpd), error=function(err) NA)
## }

## v_narea <- seq(1e-3,3e-2, length.out = n_length)
## v_vpd1 <- seq(-0.3, 0.3, length.out = n_length)
## v_vpd2 <- seq(0, 3, length.out = n_length)

## m_lcp <- matrix(NA, nrow = n_length, ncol = n_length)
## m_lcp_F<- matrix(NA, nrow = n_length, ncol = n_length)

## for (i in 1:n_length) {
## m_lcp[i, ]<- sapply(v_narea, f1, vpd = v_vpd1[i])
## m_lcp_F[i, ]<- sapply(v_narea, f2b, vpd = v_vpd2[i])
## }
## return(list(m_lcp, m_lcp_F))
## }


## plot_strategy_lcp_narea_vpd<- function(l_data){
## m_lcp <-  l_data[[1]]
## m_lcp_F <- l_data[[2]]
## par(mfrow = c(1, 2))
## image(m_lcp[rev(seq_len(nrow(m_lcp))), ], xlab = "water stress", ylab = "Narea")
## contour(m_lcp[rev(seq_len(nrow(m_lcp))), ], add =TRUE)
## image(m_lcp_F, xlab = "water stress", ylab = "Narea")
## contour(m_lcp_F, add =TRUE)
## }





## #####################################################
## #####################################################
## ### COMPARE FF16 and FF16FvCB along stress graident vpd

## demo_ver_site_lma_light <- function(fun_param, .site_prod = 1,  lma = 1,
##                                    light = 1.0){
## hyper <- fun_param(site_prod = .site_prod)
## s <- strategy(trait_matrix(lma, "lma"), hyper)
## env<- fixed_environment(light)
## p<- FF16_PlantPlus(s)
## p$compute_vars_phys(env)
## # compute demo
## tt <- seq(0, 70, length.out = 101)
## res<-  grow_plant_to_time(p, tt, env)
## res<-  grow_plant_to_size(pl, c(1, 20), "height", env)
## dhdt<- sapply(res$plant, function(x) x$internals$height_dt)
## mat_demo <- res$state[, c('height', 'mortality', 'fecundity')]
## mat_demo <- cbind(time = tt, mat_demo, dhdt = dhdt)
## return(mat_demo)
## }

## demo_ver_siteS_lmaS_lightS<- function(fun_param,
##                                       site_prodS = c(-0.4, 0, 0.4),
##                                       lmaS = c(0.5, 1 , 1.5),
##                                       lightS = c(0.2, 1.0)){
## names(site_prodS) <- paste0('site_prod', site_prodS)
## names(lmaS) <- paste0('lma', lmaS)
## names(lightS) <- paste0('light', lightS)
## res <- array(NA, dim = c(101, 5, length(site_prodS),
##                          length(lmaS), length(lightS)),
##              dimnames = list(NULL,
##                              c('time', 'height', 'mortality',
##                                'fecundity', 'dhdt'),
##                                names(site_prodS),
##                                names(lmaS),
##                                names(lightS)))
## for(site_prod in names(site_prodS)){
##     for(lma in names(lmaS)){
##         for(light in names(lightS)){
##            res[ , ,
##               site_prod,
##               lma,
##               light] <- demo_ver_site_lma_light(fun_param,
##                                             .site_prod = site_prodS[site_prod],
##                                             lma = lmaS[lma],
##                                             light = lightS[light])
##         }
##      }
##  }
## return(res)
## }


## plot_demo_var <- function(res, site_prod_n, var, lmaS_n, lightS_n){
##     require(RColorBrewer)
##     cols <- brewer.pal(3,'Blues')
##    plot(res[ , 'time', site_prod_n,
##                lmaS_n[1],
##                lightS_n[2]],
##         res[ , var, site_prod_n,
##                lmaS_n[2],
##                lightS_n[1]], type = "l", lwd = 2,
##         ylim = range(res[ -1 , var, , , ]),
##         col = cols[2],
##         xlab = "time", ylab = var)
##    lines(res[ , 'time', site_prod_n,
##                lmaS_n[1],
##                lightS_n[1]],
##          res[ , var, site_prod_n,
##                lmaS_n[1],
##                lightS_n[1]], col = cols[1], lwd = 2)
##    lines(res[ , 'time', site_prod_n,
##                lmaS_n[3],
##                lightS_n[1]],
##          res[ , var, site_prod_n,
##                lmaS_n[3],
##                lightS_n[1]], col = cols[3], lwd = 2)
## # low light
##    lines(res[ , 'time', site_prod_n,
##                lmaS_n[1],
##                lightS_n[2]],
##          res[ , var, site_prod_n,
##                lmaS_n[1],
##                lightS_n[2]], col = cols[1], lty = 2, lwd = 2)
##    lines(res[ , 'time', site_prod_n,
##                lmaS_n[2],
##                lightS_n[2]],
##          res[ , var, site_prod_n,
##                lmaS_n[2],
##                lightS_n[2]], col = cols[2], lty = 2, lwd = 2)
##    lines(res[ , 'time', site_prod_n,
##                lmaS_n[3],
##                lightS_n[2]],
##          res[ , var, site_prod_n,
##                lmaS_n[3],
##                lightS_n[2]], col = cols[3], lty = 2, lwd = 2)
## }

## plot_demo_site_lma_light <- function(vars, res, site_prodS, lmaS, lightS){
## site_prodS_n<- paste0('site_prod', site_prodS)
## lmaS_n<- paste0('lma', lmaS)
## lightS_n<- paste0('light', lightS)

## par(mfrow = c(length(site_prodS_n), length(vars)))
##   for (var in vars){
##    for (site_prod_n in site_prodS_n){
##        plot_demo_var(res, site_prod_n, var, lmaS_n, lightS_n)
##    }
##   }
## }


## ## res <- demo_ver_siteS_lmaS_lightS(trait_gradients_base_parameters)
## ## plot_demo_site_lma_light(vars = c("mortality", "height", "fecundity"),
## ##                          res,
## ##                          site_prodS = c(-0.4, 0, 0.4),
## ##                          lmaS = c(0.5, 1 , 1.5),
## ##                          lightS = c(0.2, 1.0))


## lcp_lma_site_prod_height <- function(fun_param, lmaS, .site_prod = 0,
##                                     height = 0.3920458, narea=1.87e-3){
##   hyper <- fun_param(site_prod = .site_prod)
##   lcp_vec <- rep(NA, length(lmaS))
##   for (i in seq_len(length(lmaS))){
##       s <- strategy(trait_matrix(c(lmaS[i], narea), c("lma", "narea")), hyper)
##       p<- FF16_PlantPlus(s)
##       p$height <- height
##       lcp_vec[i]<- lcp_whole_plant(p)
##   }
##  return(lcp_vec)
## }

## dhdt_lma_site_prod_height<- function(fun_param, lmaS, .site_prod = 0,
##                                      height = 0.3920458, narea=1.87e-3){
##   hyper <- fun_param(site_prod = .site_prod)
##   dhdt_vec <- rep(NA, length(lmaS))
##   for (i in seq_len(length(lmaS))){
##       s <- strategy(trait_matrix(c(lmaS[i], narea), c("lma", "narea")), hyper)
##       p<- FF16_PlantPlus(s)
##       p$height <- height
##       env <- fixed_environment(1)
##       p$compute_vars_phys(env)
##       dhdt_vec[i] <- p$internals$height_dt
##   }
##  return(dhdt_vec)
## }



## lcp_narea_site_prod_height <- function(fun_param, nareaS, .site_prod = 0,
##                                     height = 0.3920458, lma= 0.1978791){
##   hyper <- fun_param(site_prod = .site_prod)
##   lcp_vec <- rep(NA, length(nareaS))
##   for (i in seq_len(length(nareaS))){
##       s <- strategy(trait_matrix(c(lma, nareaS[i]), c("lma", "narea")), hyper)
##       p<- FF16_PlantPlus(s)
##       p$height <- height
##       lcp_vec[i]<- lcp_whole_plant(p)
##   }
##  return(lcp_vec)
## }

## dhdt_narea_site_prod_height<- function(fun_param, nareaS, .site_prod = 0,
##                                      height = 0.3920458, lma=0.1978791){
##   hyper <- fun_param(site_prod = .site_prod)
##   dhdt_vec <- rep(NA, length(nareaS))
##   for (i in seq_len(length(nareaS))){
##       s <- strategy(trait_matrix(c(lma, nareaS[i]), c("lma", "narea")), hyper)
##       p<- FF16_PlantPlus(s)
##       p$height <- height
##       env <- fixed_environment(1)
##       p$compute_vars_phys(env)
##       dhdt_vec[i] <- p$internals$height_dt
##   }
##  return(dhdt_vec)
## }


## lcp_dhdt_lma_site_prod_height<- function(fun_param, lmaS, .site_prod = 0,
##                                     height = 1, narea=1.87e-3){
##   dhdt_vec <- dhdt_lma_site_prod_height(fun_param, lmaS, .site_prod,
##                                         height, narea)
##   lcp_vec <- lcp_lma_site_prod_height(fun_param, lmaS, .site_prod,
##                                       height, narea)
##   return(data.frame(lcp = lcp_vec,
##                     dhdt = dhdt_vec))
## }


## plot_lcp_version <- function(param_slope_P, param_slope_TP){
##  lmaS <- seq(from = 0.05, to = 0.3, length.out = 100)
##  par(mfrow = c(1, 3))
##  for (site in c(-0.3, 0, 0.3)){
##  lcp_base <- lcp_lma_site_prod_height(trait_gradients_base_parameters,
##                                      lmaS, .site_prod = site)
##  lcp_elev <- lcp_lma_site_prod_height(trait_gradients_elev_parameters,
##                                      lmaS, .site_prod = site)
##  lcp_slope <- lcp_lma_site_prod_height(trait_gradients_slope_parameters,
##                                      lmaS, .site_prod = site)
##  plot(lmaS, lcp_base, type = 'l',
##       xlab = 'lma',
##       ylab = 'Whole-plant light compensation point',
##       ylim = c(0, 1))
##  lines(lmaS, lcp_elev, col = 'red')
##  lines(lmaS, lcp_slope, col = 'green')
##  }
## }

## plot_lma_shade_tradeoff <-  function(lmaS = seq(from = 0.06, to = 0.27, length.out = 100)){
## df_dhdt_lcp <- lcp_dhdt_lma_site_prod_height(trait_gradients_base_parameters,
##                                              lmaS = lmaS)
## par(mfcol = c(2,2))
## plot(lmaS, 1 - df_dhdt_lcp$lcp, type = 'l',
##      xlab = 'LMA',
##      ylab = '1 - light compensation point (%)', cex.lab = 1.1)
## plot(lmaS, df_dhdt_lcp$dhdt, type = 'l',
##      ylab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
##      xlab = 'LMA', cex.lab = 1.1)
## plot(df_dhdt_lcp$dhdt, 1 - df_dhdt_lcp$lcp, type = 'l',
##      xlab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
##      ylab = '1 - light compensation point (%)', cex.lab = 1.1)
## }


## photo_max_LeafN_Stress<- function(
##                                 narea=1.87e-3,
##                                 narea_0=1.87e-3,
##                                 lma_0=0.1978791,
##                                 B_kl1=0.4565855,
##                                 B_kl2=1.71,
##                                 rho_0=608.0,
##                                 B_dI1=0.01,
##                                 B_dI2=0.0,
##                                 B_ks1=0.2,
##                                 B_ks2=0.0,
##                                 B_rs1=4012.0,
##                                 B_rb1=2.0*4012.0,
##                                 B_f1 =3.0,
##                                 B_lf1=5120.738 * 1.87e-3 * 24 * 3600 / 1e+06,
##                                 B_lf2=0.5,
##                                 B_lf3=0.04,
##                                 B_lf4=21000,
##                                 B_lf5=1,
##                                 k_I=0.5,
##                                 latitude=0,
##                                 site_prod=0) {

##     ## rho / sapwood respiration relationship:

##     ## Narea, photosynthesis, respiration

##     assimilation_rectangular_hyperbolae <- function(I, Amax, theta, QY) {
##       x <- QY * I + Amax
##       (x - sqrt(x^2 - 4 * theta * QY * I * Amax)) / (2 * theta)
##     }

##     ## Photosynthesis  [mol CO2 / m2 / yr]
##     approximate_annual_assimilation <- function(narea, latitude) {
##       E <- seq(0, 1, by=0.02)
##       ## Only integrate over half year, as solar path is symmetrical
##       D <- seq(0, 365/2, length.out = 10000)
##       I <- plant:::PAR_given_solar_angle(plant:::solar_angle(D,
##                                                              latitude = abs(latitude)))

##       Amax <- (1+site_prod)*B_lf1 * (narea/narea_0) ^  B_lf5
##       theta <- B_lf2
##       QY <- B_lf3

##       AA <- NA * E

##       for (i in seq_len(length(E))) {
##         AA[i] <- 2 * plant:::trapezium(D, assimilation_rectangular_hyperbolae(
##                                     k_I * I * E[i], Amax, theta, QY))
##       }
##       if(all(diff(AA) < 1E-8)) {
##         # line fitting will fail if all have are zero, or potentially same value
##         ret <- c(last(AA), 0)
##         names(ret) <- c("p1","p2")
##       } else {
##         fit <- nls(AA ~ p1 * E/(p2 + E), data.frame(E = E, AA = AA),
##                    start = list(p1 = 100, p2 = 0.2))
##         ret <- coef(fit)
##       }
##       ret
##     }
##     # This needed in case narea has length zero, in which case trapezium fails
##     a_p1 <- a_p2 <- 0 * narea
##     ## TODO: Remove the 0.5 hardcoded default for k_I here, and deal
##     ## with this more nicely.
##     if (length(narea) > 0 || k_I != 0.5) {
##       i <- match(narea, unique(narea))
##       y <- vapply(unique(narea), approximate_annual_assimilation,
##                   numeric(2), latitude)
##       a_p1  <- y["p1", i]
##       a_p2  <- y["p2", i]
##     }


##   y <- vapply(unique(narea), approximate_annual_assimilation,
##               numeric(2), latitude)
##   return(y["p1", 1])
## }



## ## plot strategy at different vpd


## data_contour  <- function(nareaS = seq(from = 0.2, to = 10, length.out = 100)/1000,
##                            site_prodS = seq(from = -0.6, to = 0.6, length.out = 100)){
##    res_photo <- matrix(NA, nrow = 100, ncol = 100)
##    rownames(res_photo) <- paste0('narea', nareaS)
##    colnames(res_photo) <- paste0('site_prod', site_prodS)
##    for(narea in nareaS){
##      for(site_prod in site_prodS){
##       res <- photo_max_LeafN_Stress(narea = narea, site_prod = site_prod)
##       res_photo[paste0('narea', narea),
##                 paste0('site_prod', site_prod)] <-  res
##      }
##    }
##    return(res_photo)
## }

## plot_contour_waterStress_LeafN <- function(res_photo,
##                                            nareaS = seq(from = 0.2, to = 10, length.out = 100)/1000,
##                                            site_prodS = seq(from = -0.6, to = 0.6, length.out = 100)){
## contour(nareaS, site_prodS, res_photo, xlab = 'Leaf N', ylab = 'Water availability')
## }


## plot_water_LeafN <- function(){
## nareaS <-  seq(from = 1.5, to = 8, length.out = 100)/1000
## lcp_vec <- lcp_narea_site_prod_height(trait_gradients_base_parameters,
##                            nareaS = nareaS)
## dhdt_vec_hs<- dhdt_narea_site_prod_height(trait_gradients_base_parameters,
##                            nareaS = nareaS, .site_prod = -0.5)
## par(mfrow = c(2, 1))
## plot(nareaS, 1 - lcp_vec,  type = "l",
##      ylab = expression(Height~growth~'in'~full~light~low~water~(m~year^{-1})),
##      xlab = 'Narea', cex.lab = 1.2)
## plot(nareaS, dhdt_vec_hs,  type = "l",
##      xlab = 'Narea',
##      ylab = '1 - light compensation point (%)', cex.lab = 1.2)
## }



