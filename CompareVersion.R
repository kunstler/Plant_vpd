## Explore equilibrium seed rain and impact on growth and fecundity function between original version and new version based on vignettes code from MEE
library(plant)
library(parallel)
library(plant.assembly)
source('R/assembly.R')


### TEST DANIEL

### TODO
# - OTHER TRAITS for survival with respiration cost

## FOR THE PRESENTATION

## - SLA EFFECT ON SUCCESSIONAL COEXISTENCE
## - Growth in high light vs light compensation point

## - PLOT OF PREDICTION


# IS the modle of Wright et al. 2003 already implemented in plant ?? increasing narea will compensate for decrease in lower prod in dryer site but will have a respiration cost

## See Prentice to see if we can derive estimate of relative cost of leaf N and water

# Do plots of tradeoffs for slides.


library(plant)
library(plant.assembly)
source("R/assembly.R")

p <- trait_gradients_base_parameters(site_prod=0.0)
p$disturbance_mean_interval <- 40

# Let's estimate the bounds for our trait before setting up the assembler
bounds_narea <- viable_fitness(bounds(narea=c(1E-5, 1E2)), p, x = 0.01)
bounds_lma <- viable_fitness(bounds(lma=c(0.001, 3)), p, x = 0.01)

# now pass in the bounds -- previously we just passed in infinite bounds
 sys0 <- community(p, rbind(bounds_narea, bounds_lma),
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
obj_m0 <- assembler(sys0, list(birth_move_tol=0.5,
                               compute_viable_fitness = FALSE))

system.time(
ret <- assembler_run(obj_m0, 20)
)



p <- trait_gradients_base_parameters(site_prod=0.0)
p$disturbance_mean_interval <- 40

# Let's estimate the bounds for our trait before setting up the assembler
bounds_narea <- viable_fitness(bounds(narea=c(1E-5, 1E2)), p, x = 0.01)
bounds_lma <- viable_fitness(bounds(lma=c(0.001, 3)), p, x = 0.01)

# now pass in the bounds -- previously we just passed in infinite bounds
 sys0 <- community(p, bounds_narea,
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
obj_m0 <- assembler(sys0, list(birth_move_tol=0.5, compute_viable_fitness = FALSE))

system.time(
ret1 <- assembler_run(obj_m0, 20)
)

# now pass in the bounds -- previously we just passed in infinite bounds
 sys0 <- community(p, bounds_lma,
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
obj_m0 <- assembler(sys0, list(birth_move_tol=0.5, compute_viable_fitness = FALSE))

system.time(
ret2 <- assembler_run(obj_m0, 20)
)

## Explore photosynthesis model and Leaf N


## That seems similar to Wright et al. 2003 Am. Nat.


nareaS <-  seq(from = 0.2, to = 10, length.out = 3)/1000
prodS <- seq(from = -0.3, to = 0.3, length = 20)
list_t <-  vector('list')
for (i in 1:20){

list_t[[i]]<- dhdt_narea_site_prod_height(trait_gradients_base_parameters,
                           nareaS = nareaS, .site_prod = prodS[i])
}

par(mfrow = c(1, 2))
plot(nareaS, lcp_vec, type = 'l')
plot(nareaS, dhdt_vec, type = 'l')


df_dhdt_lcp_sl<- lcp_dhdt_lma_site_prod_height(trait_gradients_base_parameters,
                                             lmaS = lmaS, .site_prod = -0.3)



lmaS <-  seq(from = 0.06, to = 0.27, length.out = 100)
df_dhdt_lcp <- lcp_dhdt_lma_site_prod_height(trait_gradients_base_parameters,
                                             lmaS = lmaS)
df_dhdt_lcp_hn<- lcp_dhdt_lma_site_prod_height(trait_gradients_base_parameters,
                                             lmaS = lmaS, narea = 0.01)
df_dhdt_lcp_sl<- lcp_dhdt_lma_site_prod_height(trait_gradients_base_parameters,
                                             lmaS = lmaS, .site_prod = -0.3)
df_dhdt_lcp_hn_sl<- lcp_dhdt_lma_site_prod_height(
                                             trait_gradients_base_parameters,
                                             lmaS = lmaS, .site_prod = -0.3,
                                             narea = 0.01)
plot(df_dhdt_lcp$dhdt, 1 - df_dhdt_lcp$lcp, type = 'l',
     xlab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
     ylab = '1 - light compensation point (%)', xlim = c(0.1, 1.8),
     ylim = c(0.55, 0.9))
lines(df_dhdt_lcp_hn$dhdt, 1 - df_dhdt_lcp_hn$lcp, lty = 2)
lines(df_dhdt_lcp_hn_sl$dhdt, 1 - df_dhdt_lcp_hn_sl$lcp, lty = 2, col = 'red')
lines(df_dhdt_lcp_sl$dhdt, 1 - df_dhdt_lcp_sl$lcp, col = 'red')


plot_lma_shade_tradeoff <-  function(lmaS = seq(from = 0.06, to = 0.27, length.out = 100)){
df_dhdt_lcp <- lcp_dhdt_lma_site_prod_height(trait_gradients_base_parameters,
                                             lmaS = lmaS)
plot(df_dhdt_lcp$dhdt, 1 - df_dhdt_lcp$lcp, type = 'l',
     xlab = expression(Height~growth~'in'~full~light~(m~year^{-1})),
     ylab = '1 - light compensation point (%)')
}

## TODO PLOT LCP in function of leaf N

###########
## Competitive ability in plant are LCP and Photo max in full light or seedling growth in full light and WHAT IS STRESS TOLERANCE ?? GROWTH IN LOW PROD SITE

## NOT SURE WHAT IS THE EFFECT FOR SLOPE AND ELEV VARIATION TO PLOT AND FOR BASE MODEL



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

