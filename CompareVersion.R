## Explore equilibrium seed rain and impact on growth and fecundity function between original version and new version based on vignettes code from MEE
library(plant)
library(parallel)
library(plant.assembly)
source('R/assembly.R')


## Explore photosynthesis model and Leaf N




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

nareaS <-  seq(from = 0.2, to = 10, length.out = 100)/1000
# in Wright et al. 2004 N area vary roughly from 0.1 to 10 g per m2
site_prodS <- seq(from = -0.6, to = 0.6, length.out = 100)
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

contour(nareaS, site_prodS, res_photo)
## That seems similar to Wright et al. 2003 Am. Nat.

# IS the modle of Wright et al. 2003 already implemented in plant ?? increasing narea will compensate for decrease in lower prod in dryer site but will have a respiration cost

## See Prentice to see if we can derive estimate of relative cost of leaf N and water

lcp_vec <- lcp_narea_site_prod_height(trait_gradients_base_parameters,
                           nareaS = nareaS)
dhdt_vec <- dhdt_narea_site_prod_height(trait_gradients_base_parameters,
                           nareaS = nareaS)
par(mfrow = c(1, 2))
plot(nareaS, lcp_vec, type = 'l')
plot(nareaS, dhdt_vec, type = 'l')


lmaS <-  seq(from = 0.04, to = 0.3, length.out = 100)
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
     xlab = 'Height growth in full light',
     ylab = '1 -lcp', xlim = c(0.1, 1.8),
     ylim = c(0.55, 0.9))
lines(df_dhdt_lcp_hn$dhdt, 1 - df_dhdt_lcp_hn$lcp, lty = 2)
lines(df_dhdt_lcp_hn_sl$dhdt, 1 - df_dhdt_lcp_hn_sl$lcp, lty = 2, col = 'red')
lines(df_dhdt_lcp_sl$dhdt, 1 - df_dhdt_lcp_sl$lcp, col = 'red')

#ask Daniel Why ?

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

