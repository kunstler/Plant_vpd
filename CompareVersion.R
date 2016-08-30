## Explore equilibrium seed rain and impact on growth and fecundity function between original version and new version based on vignettes code from MEE
library(plant)
library(parallel)
source('R/assembly.R')


# PREDICT SEED PRODUCTION IN FUNCTION OF SIZE IN A GIVEN ENVIRONMENT
p0<- trait_gradients_base_parameters(site_prod=1)
s <- strategy(trait_matrix(0.1, "lma"), p0)

pl <- FF16_PlantPlus(s)
env <- fixed_environment(1.0)
pl$compute_vars_phys(env)

## plot growth and fecundity for new version with neutralised fecundity curve
tt <- seq(0, 70, length.out = 101)
res <-  grow_plant_to_time(pl, tt, env)
dhdt <- sapply(res$plant, function(x) x$internals$height_dt)
res_FF16<-  grow_plant_to_time(FF16_PlantPlus(FF16_Strategy()), tt, env)
dhdt_FF16<- sapply(res_FF16$plant, function(x) x$internals$height_dt)

par(mfrow = c(2, 2))
plot(height ~ tt, res$state, type = "l", las = 1,
     xlab = "Time (years)", ylab = "Height (m)",
     ylim = range(res$state[, 'height']))
plot(tt, dhdt, type = "l", las = 1,
     xlab = "Time (years)", ylab = "Height growth rate(m / yr)",
     ylim = range(dhdt, dhdt_FF16))
plot(res$state[, 'height'], 1 - 0.5/(1+exp(0*(1-res$state[, 'height']/s$hmat))),
     xlab = "Height (m)", ylab = "dMa_dB Growth fraction",
     type = "l", ylim = c(0, 1))
plot(fecundity ~ tt[-1], res$state[-1, ], type = "l", las = 1,
     xlab = "Time (years)", ylab = "Fecundity (seed number)", log = "y",
     ylim = range(res$state[-1, 'fecundity'],
                  res_FF16$state[-1, 'fecundity']))

x11()
par(mfrow = c(2, 2))
plot(height ~ tt, res_FF16$state, type = "l", las = 1,
     xlab = "Time (years)", ylab = "Height (m)",
     ylim = range(res$state[, 'height']))
plot(tt, dhdt_FF16, type = "l", las = 1,
     xlab = "Time (years)", ylab = "Height growth rate(m / yr)",
     ylim = range(dhdt, dhdt_FF16))
plot(res_FF16$state[, 'height'], 1 - 1/(1+exp(50*(1-res_FF16$state[, 'height']/s$hmat))),
     xlab = "Height (m)", ylab = "dMa_dB Growth fraction",
     type = "l", ylim = c(0, 1))
plot(fecundity ~ tt[-1], res_FF16$state[-1, ], type = "l", las = 1,
     xlab = "Time (years)", ylab = "Fecundity (seed number)", log = "y",
     ylim = range(res$state[-1, 'fecundity'],
                  res_FF16$state[-1, 'fecundity']))


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

p0_s<- trait_gradients_base_parameters(site_prod=1)

p0_s$disturbance_mean_interval <- 30

p_s<- expand_parameters(trait_matrix(0.0825, "lma"), p0_s, FALSE)

## Some output seed rains, given an input seed rain:
run(1.0, p_s)
run(241.0, p_s)
run(369.0, p_s)
run(400.0, p_s)
run(384.0, p_s)


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

