## Explore equilibrium seed rain and impact on growth and fecundity function between original version and new version based on vignettes code from MEE
library(plant)
library(parallel)
source('R/assembly.R')


## explore alternative model of plant growth along gradients

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
dhdt<- sapply(res$plant, function(x) x$internals$height_dt)
mat_demo <- res$state[, c('height', 'mortality', 'fecundity')]
mat_demo <- cbind(time = tt, mat_demo, dhdt = dhdt)
return(mat_demo)
}

test <-  demo_ver_site_lma_light(fun_param =
                                 trait_gradients_slope_elev_parameters,
                                 .site_prod = 0.5,
                                 lma = 0.2,
                                 light = 0.5)

test_b<-  demo_ver_site_lma_light(fun_param =
                                 trait_gradients_base_parameters,
                                 .site_prod = 0.5,
                                 lma = 0.2,
                                 light = 0.5)
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
plot(test[ , 'time'], test[ , 'dhdt'], type = "l",
     ylim = range(test[ , 'dhdt'], test_b[ , 'dhdt']))
lines(test_b[ , 'time'], test_b[ , 'dhdt'], col = 'red')

plot(test[ , 'time'], test[ , 'height'], type = "l",
     ylim = range(test[ , 'height'], test_b[ , 'height']))
lines(test_b[ , 'time'], test_b[ , 'height'], col = 'red')

plot(test[ , 'time'], test[ , 'mortality'], type = "l",
     ylim = range(test[ , 'mortality'], test_b[ , 'mortality']))
lines(test_b[ , 'time'], test_b[ , 'mortality'], col = 'red')

plot(test[ , 'time'], test[ , 'fecundity'], type = "l",
     ylim = range(test[ , 'fecundity'], test_b[ , 'fecundity']))
lines(test_b[ , 'time'], test_b[ , 'fecundity'], col = 'red')

demo_ver_siteS_lmaS_lightS<- function(fun_param,
                                      site_prodS = c(0.8, 1, 1.2),
                                      lmaS = c(0.5, 1 , 1.5),
                                      lightS = c(0.2, 1.0)){
res <- array(NA, dim = c(101, 5, length(site_prodS),
                         length(lmaS), length(lightS)),
             dimnames = list(NULL,
                             c('time', 'height', 'mortality',
                               'fecundity', 'dhdt'),
                             site_prodS, lmaS, lightS))
for(site_prod in site_prodS){
    for(lma in lmaS){
        for(light in lightS){
           res[ , ,
              as.character(site_prod),
              as.character(lma),
              as.character(light)] <- demo_ver_site_lma_light(fun_param,
                                                         .site_prod = site_prod,
                                                         lma,
                                                         light)
        }
     }
 }
return(res)
}



plot_demo_site_lma_light <- function(var, res, site_prodS, lmaS, lightS){
   for (site_prod in site_prodS){
   plot(res[ , 'time', as.character(site_prod),
               as.character(lmaS[1]),
               as.character(lightS[2])],
        res[ , var, as.character(site_prod),
               as.character(lmaS[2]),
               as.character(lightS[1])], type = "l", las = 1,
        xlab = NA, ylab = NA,
        ylim = range(res[ -1 , var, , , ]), col = rev(topo.colors(6))[4])
   lines(res[ , 'time', as.character(site_prod),
               as.character(lmaS[1]),
               as.character(lightS[1])],
         res[ , var, as.character(site_prod),
               as.character(lmaS[1]),
               as.character(lightS[1])], col = rev(topo.colors(6))[1])
   lines(res[ , 'time', as.character(site_prod),
               as.character(lmaS[3]),
               as.character(lightS[1])],
         res[ , var, as.character(site_prod),
               as.character(lmaS[3]),
               as.character(lightS[1])], col = rev(topo.colors(6))[6])

   lines(res[ , 'time', as.character(site_prod),
               as.character(lmaS[1]),
               as.character(lightS[2])],
         res[ , var, as.character(site_prod),
               as.character(lmaS[1]),
               as.character(lightS[2])], col = rev(topo.colors(6))[1], lty = 2)
   lines(res[ , 'time', as.character(site_prod),
               as.character(lmaS[2]),
               as.character(lightS[2])],
         res[ , var, as.character(site_prod),
               as.character(lmaS[2]),
               as.character(lightS[2])], col = rev(topo.colors(6))[4], lty = 2)
   lines(res[ , 'time', as.character(site_prod),
               as.character(lmaS[3]),
               as.character(lightS[2])],
         res[ , var, as.character(site_prod),
               as.character(lmaS[3]),
               as.character(lightS[2])], col = rev(topo.colors(6))[6], lty = 2)
   }
}

res_SE <- demo_ver_siteS_lmaS_lightS(fun_param =
                                     trait_gradients_slope_elev_parameters,
                                     site_prodS = c( -0.5, 0.0, 0.5),
                                     lmaS = c(0.0825, 0.18 , 0.2625),
                                     lightS = c(0.5, 1.0))

res_b<- demo_ver_siteS_lmaS_lightS(fun_param =
                                   trait_gradients_base_parameters,
                                   site_prodS = c(-0.5, 0.0, 0.5),
                                   lmaS = c(0.0825, 0.18 , 0.2625),
                                   lightS = c(0.5, 1.0))


##
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
plot_demo_site_lma_light(var = 'dhdt', res_SE,
                         site_prodS = c(-0.5, 0.0, 0.5),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'mortality', res_SE,
                         site_prodS = c(-0.5, 0.0, 0.5),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'fecundity', res_SE,
                         site_prodS = c(-0.5, 0.0, 0.5),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))




x11()
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
plot_demo_site_lma_light(var = 'dhdt', res_b,
                         site_prodS = c(-0.5, 0.0, 0.5),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'mortality', res_b,
                         site_prodS = c(-0.5, 0.0, 0.5),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'fecundity', res_b,
                         site_prodS = c(-0.5, 0.0, 0.5),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))


par(mfrow = c(2, 3), mar = c(1, 1, 1, 1))
plot_demo_site_lma_light(var = 'dhdt', res_SE,
                         site_prodS = c(-0.5, 0.0, 0.5),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))
plot_demo_site_lma_light(var = 'dhdt', res_SE_b,
                         site_prodS = c(-0.5, 0.0, 0.5),
                         lmaS = c(0.0825, 0.18 , 0.2625),
                         lightS = c(0.5, 1.0))






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

