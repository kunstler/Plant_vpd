## Explore equilibrium seed rain and impact on growth and fecundity function between original version and new version based on vignettes code from MEE
library(plant)
library(parallel)
source('R/assembly.R')


## explore alternative model of plant growth along gradients


test_elev<-  demo_ver_site_lma_light(fun_param =
                                 trait_gradients_slope_parameters,
                                 .site_prod = -0.3,
                                 lma = 0.08,
                                 light = 0.5)

test_base<-  demo_ver_site_lma_light(fun_param =
                                 trait_gradients_base_parameters,
                                 .site_prod = -0.3,
                                 lma = 0.08,
                                 light = 0.5)

quartz()
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1))
plot(test_elev[ , 'time'], test_elev[ , 'dhdt'], type = "l",
     ylim = range(test_elev[ , 'dhdt'], test_base[ , 'dhdt']))
title(xlab = 'time', ylab = 'dhdt', line=2, cex.lab=1.2)
lines(test_base[ , 'time'], test_base[ , 'dhdt'], col = 'red')

plot(test_elev[ , 'time'], test_elev[ , 'height'], type = "l",
     ylim = range(test_elev[ , 'height'], test_base[ , 'height']))
title(xlab = 'time', ylab = 'height', line=2, cex.lab=1.2)
lines(test_base[ , 'time'], test_base[ , 'height'], col = 'red')

plot(test_elev[ , 'time'], test_elev[ , 'mortality'], type = "l",
     ylim = range(test_elev[ , 'mortality'], test_base[ , 'mortality']))
title(xlab = 'time', ylab = 'mortality', line=2, cex.lab=1.2)
lines(test_base[ , 'time'], test_base[ , 'mortality'], col = 'red')

plot(test_elev[ , 'time'], test_elev[ , 'fecundity'], type = "l",
     ylim = range(test_elev[ , 'fecundity'], test_base[ , 'fecundity']))
title(xlab = 'time', ylab = 'fecundity', line=2, cex.lab=1.2)
lines(test_base[ , 'time'], test_base[ , 'fecundity'], col = 'red')


res_elev <- demo_ver_siteS_lmaS_lightS(fun_param =
                                     trait_gradients_elev_parameters,
                                     site_prodS = c( -0.3, 0.0, 0.3),
                                     lmaS = c(0.0825, 0.18 , 0.2625),
                                     lightS = c(0.5, 1.0))

res_slope <- demo_ver_siteS_lmaS_lightS(fun_param =
                                     trait_gradients_slope_parameters,
                                     site_prodS = c( -0.3, 0.0, 0.3),
                                     lmaS = c(0.0825, 0.18 , 0.2625),
                                     lightS = c(0.5, 1.0))

res_base<- demo_ver_siteS_lmaS_lightS(fun_param =
                                   trait_gradients_base_parameters,
                                   site_prodS = c(-0.3, 0.0, 0.3),
                                   lmaS = c(0.0825, 0.18 , 0.2625),
                                   lightS = c(0.5, 1.0))


##
quartz()
par(mfrow = c(3, 3), mar = c(1, 2, 1, 1))
plot_demo_site_lma_light(var = 'dhdt', res_elev,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'mortality', res_elev,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'fecundity', res_elev,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))


quartz()
par(mfrow = c(3, 3), mar = c(1, 2, 1, 1))
plot_demo_site_lma_light(var = 'dhdt', res_slope,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'mortality', res_slope,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'fecundity', res_slope,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))


quartz()
par(mfrow = c(3, 3), mar = c(1, 2, 1, 1))
plot_demo_site_lma_light(var = 'dhdt', res_base,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'mortality', res_base,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'fecundity', res_base,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))

quartz()
par(mfrow = c(2, 3), mar = c(1, 2, 1, 1))
plot_demo_site_lma_light(var = 'dhdt', res_elev,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'dhdt', res_base,
                         site_prodS = c(-0.3, 0.0, 0.3),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))


### NEED TO LOOK AT LIGHT COMPENSATION POINT
plot_lcp_version()
#

### IMPACT ON POP DYNAMICS

# run_scm_collect see https://github.com/traitecoevo/plant_paper/blob/master/R/figures_patch.R

### EXPLORE IMPACT on POPULATION EQUILIBRIUM SEED PRODUCTION

run <- function(seed_rain_in, p) {
  p$seed_rain <- seed_rain_in
  run_scm(p)$seed_rains
}

cobweb <- function(m, ...) {
  lines(rep(m[,1], each=2), c(t(m)), ...)
}


run_new_schedule <- function(seed_rain_in, p) {
  p$seed_rain <- seed_rain_in
  res <- build_schedule(p)
  attr(res, "seed_rain_out")
}


## Set model parameters original version
p0 <- scm_base_parameters("FF16")
p <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FALSE)

## Some output seed rains, given an input seed rain:
run(1.0, p)
run(10.0, p)
run(50.0, p)

## # 1: Approach to equilibrium:
p_eq <- equilibrium_seed_rain(p)

## Equilibrium seed rain is right around where the last two values
## from run were producing
p_eq$seed_rain

## From a distance, these both hone in nicely on the equilibrium, and
## rapidly, too.
##+ equilibrium_approach
approach <- attr(p_eq, "progress")
r <- range(approach)
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach, pch=19, cex=.5, type="o")


## Set model parameters with site_prod set to 1

p0_s<- trait_gradients_base_parameters(site_prod=0.1)

p0_s$disturbance_mean_interval <- 30

p_s<- expand_parameters(trait_matrix(0.0825, "lma"), p0_s, FALSE)



## Some output seed rains, given an input seed rain:
run(1.0, p_s)
run(241.0, p_s)
run(369.0, p_s)
run(490.0, p_s)
run(481.0, p_s)
run(478.0, p_s)


## # 1: Approach to equilibrium:
p_s_eq <- equilibrium_seed_rain(p_s)

## Equilibrium seed rain is right around where the last two values
## from run were producing
p_s_eq$seed_rain

## From a distance, these both hone in nicely on the equilibrium, and
## rapidly, too.
##+ equilibrium_approach
approach <- attr(p_s_eq, "progress")
r <- range(approach)
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach, pch=19, cex=.5, type="o")



## Set model parameters with site_prod set to 1.2

p0_s<- trait_gradients_base_parameters(site_prod=1.2)

p0_s$disturbance_mean_interval <- 30

p_s<- expand_parameters(trait_matrix(0.0825, "lma"), p0_s, FALSE)

## Some output seed rains, given an input seed rain:
run(1.0, p_s)
run(352.0, p_s)
run(571.0, p_s)
run(700.0, p_s)
run(658.0, p_s)


## # 1: Approach to equilibrium:
p_s_eq <- equilibrium_seed_rain(p_s)

## Equilibrium seed rain is right around where the last two values
## from run were producing
p_s_eq$seed_rain

## From a distance, these both hone in nicely on the equilibrium, and
## rapidly, too.
##+ equilibrium_approach
approach <- attr(p_s_eq, "progress")
r <- range(approach)
plot(approach, type="n", las=1, xlim=r, ylim=r)
abline(0, 1, lty=2, col="grey")
cobweb(approach, pch=19, cex=.5, type="o")

