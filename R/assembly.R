
plot_trait_gradient <- function(list_assembly_lma, vec_site_prod, var = "lma",
                                varlab = expression(paste("Leaf-mass per area (kg ", m^-2,")"))) {
  list_data <- list_assembly_lma
  grad <- vec_site_prod
  y <- lapply(list_data, function(x) x$community$traits[ ,var])

  plot(NA, type="n", log="y", las=1, xlim=c(-0.3, 0.3),
       ylim= range(unlist(y)),
    ylab= varlab,
    xlab="Site productivity")


  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], y[[i]], type='p', col="black", pch=16)
  }
}

plot_trait_gradient_narea_lma<- function(list_assembly_lma,
                                         vec_site_prod,
                                         var = "narea",
               varlab = expression(paste("Nitrogen per area (kg ", m^-2,")"))) {
  list_data <- list_assembly_lma
  grad <- vec_site_prod
  y <- lapply(list_data, function(x) x$community$traits[ ,var])
  z <- lapply(list_data, function(x) x$community$traits[ ,"lma"])


  plot(NA, type="n", log="y", las=1, xlim=c(-0.33, 0.33),
       ylim= range(unlist(y)),
    ylab= varlab,
    xlab="Site productivity")

  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], y[[i]], type='p', col="black",
           pch=1, cex = 0.5+10*z[[i]])
  }
}


plot_cor_narea_lma<- function(list_assembly_lma,
                                         vec_site_prod,
                                         var = "narea") {
  list_data <- list_assembly_lma
  y <- lapply(list_data, function(x) x$community$traits[ ,var])
  z <- lapply(list_data, function(x) x$community$traits[ ,"lma"])

  plot(unlist(y), unlist(z), type="n", log="y", las=1,
    xlab = expression(paste("Nitrogen per area (kg ", m^-2,")")),
    ylab = expression(paste("Leaf-mass per area (kg ", m^-2,")")))
colfunc <- colorRampPalette(c("red", "green"))
cols <- colfunc(length(list_data))
  for(i in seq_along(list_data)) {
    points(y[[i]] , z[[i]], type='p', col=cols[i],
           pch=16)
  }
}


plot_trait_gradient_narea_lma2<- function(list_assembly_lma,
                                         vec_site_prod,
                                         var = "narea",
               varlab = expression(paste("Nitrogen per area (kg ", m^-2,")"))) {
  list_data <- list_assembly_lma
  grad <- vec_site_prod
  y <- lapply(list_data, function(x) x$community$traits[ ,var])
  z <- lapply(list_data, function(x) x$community$traits[ ,"lma"])

  par(mfrow = c(2,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)
  plot(NA, type="n", log="y", las=1, xlim=c(-0.33, 0.33),
       ylim= range(unlist(y)),  axes = FALSE, ylab = varlab)
  axis(side = 1, labels = FALSE, at = grad)
  axis(side = 2)
  box()

  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], y[[i]], type='p', col="black",
           pch=16)
  }
  plot(NA, type="n", log="y", las=1, xlim=c(-0.33, 0.33), ylim= range(unlist(z)),ylab= expression(paste("Leaf-mass per area (kg ", m^-2,")")),
    xlab=NA)

  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], z[[i]], type='p', col="black",
           pch=16)
  }

title(xlab = "Site productivity",
      outer = TRUE, line = 3)
}


plot_fitness_assembly <- function(list_assembly_lma, vec_site_prod) {
  list_data <- list_assembly_lma
  grad <- vec_site_prod


  cols <- rev(brewer.pal(5, "Blues"))

  plot(NA, log="x", xlim=c(0.01, 1), ylim=c(-2,0.5),
    xlab= expression(paste("Leaf-mass per area (kg ", m^-2,")")),
    ylab="Fitness")
  abline(h=0, col="grey")

  for(i in seq_along(list_data)) {
    data <- list_data[[i]]
    community <- last(data$history)
    ff <- community_fitness_approximate(community)
    lma <- sort(c(seq_log_range(data$community$bounds, 400), data$community$traits))
    w <- ff(lma)
    lma_res <- data$community$traits[,"lma"]
    points(lma,w, type='l', col=cols[i])
    points(lma_res,ff(lma_res), type='p', col=cols[i], pch=19)
  }
  legend("topleft", bty="n", legend=c("low", "med", "high"), col=cols, pch=19)
}


run_assembly_list<- function(vec_site_prod,FUN = run_assembly,
                             disturbance_mean_interval=10, name_data_param = NA) {
require(parallel)
mclapply(vec_site_prod,
         FUN,
         disturbance_mean_interval = disturbance_mean_interval,
         name_data_param = name_data_param,
         mc.cores = 7)
}

run_assembly_vpd_list<- function(vec_vpd,
                             disturbance_mean_interval=10,
                             FUN = run_assembly_FvCB_narea_lma) {
require(parallel)
START <- Sys.time()
res <- mclapply(vec_vpd,
                FUN,
                disturbance_mean_interval = disturbance_mean_interval,
                mc.cores = 2)
END <- Sys.time()
print(END-START)
return(res)
}


run_assembly <- function(site_prod=0.0, disturbance_mean_interval=10, ...) {

  p <- trait_gradients_base_parameters(site_prod=site_prod)

  p$disturbance_mean_interval <- disturbance_mean_interval
  sys0 <- community(p, bounds_infinite("lma"),
                     fitness_approximate_control=list(type="gp"))

  obj_m0 <- assembler(sys0, list(birth_move_tol=0.5))
  assembler_run(obj_m0, 20)
}

run_assembly_narea<- function(site_prod=0.0, disturbance_mean_interval=10, ...) {

  p <- trait_gradients_base_parameters(site_prod=site_prod)
  # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0

  p$disturbance_mean_interval <- disturbance_mean_interval

  bounds_narea <- viable_fitness(bounds(narea=c(1E-5, 1E2)), p, x = 0.01)

# now pass in the bounds -- previously we just passed in infinite bounds
  sys0 <- community(p, bounds_narea,
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
obj_m0 <- assembler(sys0, list(birth_move_tol=0.5, compute_viable_fitness = FALSE))

assembler_run(obj_m0, 20)

}

run_assembly_narea_lma<- function(site_prod=0.0, disturbance_mean_interval=10, ...) {

  p <- trait_gradients_base_parameters(site_prod=site_prod)
  # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0

  p$disturbance_mean_interval <- disturbance_mean_interval

  bounds_narea <- viable_fitness(bounds(narea=c(1E-5, 1E2)), p, x = 0.01)
  bounds_lma <- viable_fitness(bounds(lma=c(0.001, 3)), p, x = 0.01)

# now pass in the bounds -- previously we just passed in infinite bounds
  sys0 <- community(p, rbind(bounds_narea, bounds_lma),
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
obj_m0 <- assembler(sys0, list(birth_move_tol=0.5, compute_viable_fitness = FALSE))

assembler_run(obj_m0, 20)

}


run_assembly_FvCB_narea_lma<- function(vpd=0.0, disturbance_mean_interval=10, ...) {
  p <- trait_gradients_FvCB_parameters(vpd=vpd)

  p$disturbance_mean_interval <- disturbance_mean_interval
  # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0

  sys0 <- community(p, rbind(bounds(narea=c(1E-4, 1E-2)), bounds(lma=c(0.001, 3))),
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
obj_m0 <- assembler(sys0, list(birth_move_tol=0.5, compute_viable_fitness = FALSE))
print('initialized')
assembler_run(obj_m0, 20)
}


run_assembly_FvCB_narea_lma_Tleaf<- function(vpd=0.0, name_param_Tleaf, disturbance_mean_interval=10, ...) {
  coefs <- readRDS(name_param_Tleaf)
  Tleaf <- coefs['intercept'] + coefs['intercept']*log(-(vpd - coefs['max_trans']))

  p <- trait_gradients_FvCB_parameters(vpd=vpd,Tleaf = Tleaf)

  p$disturbance_mean_interval <- disturbance_mean_interval
  # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0

  sys0 <- community(p, rbind(bounds(narea=c(1E-4, 1E-2)), bounds(lma=c(0.001, 3))),
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
obj_m0 <- assembler(sys0, list(birth_move_tol=0.5, compute_viable_fitness = FALSE))
print('initialized')
assembler_run(obj_m0, 20)
}


run_assembly_elev_slope <- function(site_prod=0.0, name_data_param, disturbance_mean_interval=10) {
  p <- trait_gradients_elev_slope_parameters(name_data_param = name_data_param, site_prod = site_prod)
  # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0

  p$disturbance_mean_interval <- disturbance_mean_interval
  sys0 <- community(p, bounds_infinite("lma"),
                     fitness_approximate_control=list(type="gp"))

  obj_m0 <- assembler(sys0, list(birth_move_tol=0.5))
  assembler_run(obj_m0, 20)
}

##' Hopefully sensible set of parameters for use with the EBT.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use.
##' @title Sensible, fast (ish) EBT parameters
##' @authorDaniel Falster
##' @export
trait_gradients_base_parameters <- function(...) {
  #plant_log_console()
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.001
  ctrl$equilibrium_eps <- 1e-5

  ctrl$equilibrium_nsteps  <- 80
  ctrl$equilibrium_solver_name <- "iteration"
  ctrl$equilibrium_verbose <-  TRUE

  FF16_trait_gradient_hyperpar <- make_FF16_trait_gradient_hyperpar(...)
  p <- FF16_Parameters(patch_area=1.0, control=ctrl,
                   hyperpar=FF16_trait_gradient_hyperpar)
  p
}


##' Hopefully sensible set of parameters for use with the EBT.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use.
##' @title Sensible, fast (ish) EBT parameters
##' @authorDaniel Falster
##' @export
trait_gradients_FvCB_parameters <- function(...) {
  #plant_log_console()
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.001
  ctrl$equilibrium_eps <- 1e-5

  ctrl$equilibrium_nsteps  <- 80
  ctrl$equilibrium_solver_name <- "iteration"
  ## ctrl$equilibrium_solver_name <- "hybrid"
                                 #"hybrid" # in default this is "iteration"
  ctrl$equilibrium_verbose <-  TRUE

  FF16_FvCB_trait_gradient_hyperpar <- make_FF16FvCB_hyperpar(...)
  p <- FF16FvCB_Parameters(patch_area=1.0, control=ctrl,
                   hyperpar=FF16_FvCB_trait_gradient_hyperpar)
  p
}



##' hyperparameterisation that include a variation of LMA LT elevation and or slope of the tradeoff changing with MAP or MAT over MAP
##' @title Sensible, fast (ish) EBT parameters
##' @author Georges Kunstler
##' @export
trait_gradients_elev_slope_parameters<- function(name_data_param, ...) {
  #plant_log_console()
  param <- read.csv(name_data_param)

  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.001
  ctrl$equilibrium_eps <- 1e-5

  ctrl$equilibrium_nsteps  <- 80
  ctrl$equilibrium_solver_name <- "iteration" #"hybrid" # in default this is "iteration"
  ctrl$equilibrium_verbose <-  TRUE

B_kl1_1_log10b <- log10(0.4565855)
B_kl1_2_log10b <- 0
B_kl2_1b <- 1.71
B_kl2_2b <- 0

if (!is.na(param[1, 'LMAelev'])){
B_kl1_1_log10b <- param[1, 'LMAelev']
B_kl1_2_log10b <- param[2, 'LMAelev']*(2/0.6)
}
if (!is.na(param[1, 'LMAslope'])){
B_kl2_1b <- -param[1, 'LMAslope']
B_kl2_2b <- -param[2, 'LMAslope']*(2/0.6)
}
#this slope is for a gradient of precipitation from -1 to 1 but we can not vary that much the productivity (otherwise growth become zero) so I rescale it to have the same slope on a range only between -0.3 and 0.3 (which is a range where growth is non zero)

  FF16_trait_gradient_hyperpar <-
      make_FF16_trait_gradient_slope_elev_hyperpar(...,
                                B_kl1_1_log10 = B_kl1_1_log10b,
                                B_kl1_2_log10 = B_kl1_2_log10b,
                                B_kl2_1 = B_kl2_1b,
                                B_kl2_2 = B_kl2_2b)
  p <- FF16_Parameters(patch_area=1.0, control=ctrl,
                   hyperpar=FF16_trait_gradient_hyperpar)
  p
}




##' Hyperparameters for FF16 physiological model
##' @title Hyperparameters for FF16 physiological model
##' @param lma_0 Central (mean) value for leaf mass per area [kg /m2]
##' @param B_kl1 Rate of leaf turnover at lma_0 [/yr]
##' @param B_kl2 Scaling slope for phi in leaf turnover [dimensionless]
##' @param rho_0 Central (mean) value for wood density [kg /m3]
##' @param B_dI1 Rate of instantaneous mortality at rho_0 [/yr]
##' @param B_dI2 Scaling slope for wood density in intrinsic mortality [dimensionless]
##' @param B_ks1 Rate of sapwood turnover at rho_0 [/yr]
##' @param B_ks2 Scaling slope for rho in sapwood turnover [dimensionless]
##' @param B_rs1 CO_2 respiration per unit sapwood volume [mol / yr / m3 ]
##' @param B_rb1 CO_2 respiration per unit sapwood volume [mol / yr / m3 ]
##' @param B_f1 Cost of seed accessories per unit seed mass [dimensionless]
##' @param narea nitrogen per leaf area [kg / m2]
##' @param narea_0 central (mean) value for nitrogen per leaf area [kg / m2]
##' @param B_lf1 Potential CO_2 photosynthesis at average leaf nitrogen [mol / d / m2]
##' @param B_lf2 Curvature of leaf photosynthetic light response curve [dimensionless]
##' @param B_lf3 Quantum yield of leaf photosynthetic light response curve [dimensionless]
##' @param B_lf4 CO_2 respiration per unit leaf nitrogen [mol / yr / kg]
##' @param B_lf5 Scaling exponent for leaf nitrogen in maximum leaf photosynthesis [dimensionless]
##' @param k_I light extinction coefficient [dimensionless]
##' @param latitude degrees from equator (0-90), used in solar model [deg]
##' @param site_prod productivity/aridity of the site from -1 to 1.
##' @export
##' @rdname FF16_hyperpar
make_FF16_trait_gradient_hyperpar <- function(
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
                                narea=1.87e-3,
                                narea_0=1.87e-3,
                                B_lf1=5120.738 * 1.87e-3 * 24 * 3600 / 1e+06,
                                B_lf2=0.5,
                                B_lf3=0.04,
                                B_lf4=21000,
                                B_lf5=1,
                                k_I=0.5,
                                latitude=0,
                                site_prod=0) {
  assert_scalar <- function(x, name=deparse(substitute(x))) {
    if (length(x) != 1L) {
      stop(sprintf("%s must be a scalar", name), call. = FALSE)
    }
  }
  assert_scalar(lma_0)
  assert_scalar(B_kl1)
  assert_scalar(B_kl2)
  assert_scalar(rho_0)
  assert_scalar(B_dI1)
  assert_scalar(B_dI2)
  assert_scalar(B_ks1)
  assert_scalar(B_ks2)
  assert_scalar(B_rs1)
  assert_scalar(B_rb1)
  assert_scalar(B_f1)
  assert_scalar(narea)
  assert_scalar(narea_0)
  assert_scalar(B_lf1)
  assert_scalar(B_lf2)
  assert_scalar(B_lf3)
  assert_scalar(B_lf4)
  assert_scalar(B_lf5)
  assert_scalar(k_I)
  assert_scalar(latitude)
  assert_scalar(site_prod)

  function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }
    lma       <- with_default("lma")
    rho       <- with_default("rho")
    omega     <- with_default("omega")
    narea     <- with_default("narea", narea)

    ## lma / leaf turnover relationship:
    k_l   <- B_kl1 * (lma / lma_0) ^ (-B_kl2)

    ## rho / mortality relationship:
    d_I  <- B_dI1 * (rho / rho_0) ^ (-B_dI2)

    ## rho / wood turnover relationship:
    k_s  <- B_ks1 *  (rho / rho_0) ^ (-B_ks2)

    ## rho / sapwood respiration relationship:

    ## Respiration rates are per unit mass, so this next line has the
    ## effect of holding constant the respiration rate per unit volume.
    ## So respiration rates per unit mass vary with rho, respiration
    ## rates per unit volume don't.
    r_s <- B_rs1 / rho
    # bark respiration follows from sapwood
    r_b <- B_rb1 / rho

    ## omega / accessory cost relationship
    a_f3 <- B_f1 * omega

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
      I <- plant:::PAR_given_solar_angle(plant:::solar_angle(D, latitude = abs(latitude)))

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

    ## Respiration rates are per unit mass, so convert to mass-based
    ## rate by dividing with lma
    ## So respiration rates per unit mass vary with lma, while
    ## respiration rates per unit area don't.
    r_l  <- B_lf4 * narea / lma

    extra <- cbind(k_l,                # lma
                   d_I, k_s, r_s, r_b, # rho
                   a_f3,               # omega
                   a_p1, a_p2,         # narea
                   r_l)                # lma, narea

    overlap <- intersect(colnames(m), colnames(extra))
    if (length(overlap) > 0L) {
      stop("Attempt to overwrite generated parameters: ",
           paste(overlap, collapse=", "))
    }

    ## Filter extra so that any column where all numbers are with eps
    ## of the default strategy are not replaced:
    if (filter) {
      if (nrow(extra) == 0L) {
        extra <- NULL
      } else {
        pos <- diff(apply(extra, 2, range)) == 0
        if (any(pos)) {
          eps <- sqrt(.Machine$double.eps)
          x1 <- extra[1, pos]
          x2 <- unlist(s[names(x1)])
          drop <- abs(x1 - x2) < eps & abs(1 - x1/x2) < eps
          if (any(drop)) {
            keep <- setdiff(colnames(extra), names(drop)[drop])
            extra <- extra[, keep, drop=FALSE]
          }
        }
      }
    }

    if (!is.null(extra)) {
      m <- cbind(m, extra)
    }
    m
  }
}



##' Hyperparameters for FF16 physiological model
##' @title Hyperparameters for FF16 physiological model
##' @param lma_0 Central (mean) value for leaf mass per area [kg /m2]
##' @param B_kl1_1_log10 Mean Rate of leaf turnover at lma_0 [/yr] in log10
##' @param B_kl1_1_log10 slope of the change Rate of leaf turnover at lma_0 [/yr] with site prod in log10
##' @param B_kl2_1 Mean scaling slope parameter for phi in leaf turnover [dimensionless]
##' @param B_kl2_2 Slope of change scaling slope parameter for phi in leaf turnover [dimensionless] wioth site prod
##' @param rho_0 Central (mean) value for wood density [kg /m3]
##' @param B_dI1 Rate of instantaneous mortality at rho_0 [/yr]
##' @param B_dI2 Scaling slope for wood density in intrinsic mortality [dimensionless]
##' @param B_ks1 Rate of sapwood turnover at rho_0 [/yr]
##' @param B_ks2 Scaling slope for rho in sapwood turnover [dimensionless]
##' @param B_rs1 CO_2 respiration per unit sapwood volume [mol / yr / m3 ]
##' @param B_rb1 CO_2 respiration per unit sapwood volume [mol / yr / m3 ]
##' @param B_f1 Cost of seed accessories per unit seed mass [dimensionless]
##' @param narea nitrogen per leaf area [kg / m2]
##' @param narea_0 central (mean) value for nitrogen per leaf area [kg / m2]
##' @param B_lf1 Potential CO_2 photosynthesis at average leaf nitrogen [mol / d / m2]
##' @param B_lf2 Curvature of leaf photosynthetic light response curve [dimensionless]
##' @param B_lf3 Quantum yield of leaf photosynthetic light response curve [dimensionless]
##' @param B_lf4 CO_2 respiration per unit leaf nitrogen [mol / yr / kg]
##' @param B_lf5 Scaling exponent for leaf nitrogen in maximum leaf photosynthesis [dimensionless]
##' @param k_I light extinction coefficient [dimensionless]
##' @param latitude degrees from equator (0-90), used in solar model [deg]
##' @param site_prod productivity/aridity of the site from -1 to 1.
##' @export
##' @rdname FF16_hyperpar
make_FF16_trait_gradient_slope_elev_hyperpar <- function(
                                lma_0 = 0.1978791,
                                B_kl1_1_log10 = log10(0.4565855),
                                B_kl1_2_log10 = 0,
                                B_kl2_1 = 1.71,
                                B_kl2_2 = 0,
                                rho_0=608.0,
                                B_dI1=0.01,
                                B_dI2=0.0,
                                B_ks1=0.2,
                                B_ks2=0.0,
                                B_rs1=4012.0,
                                B_rb1=2.0*4012.0,
                                B_f1 =3.0,
                                narea=1.87e-3,
                                narea_0=1.87e-3,
                                B_lf1=5120.738 * 1.87e-3 * 24 * 3600 / 1e+06,
                                B_lf2=0.5,
                                B_lf3=0.04,
                                B_lf4=21000,
                                B_lf5=1,
                                k_I=0.5,
                                latitude=0,
                                site_prod=0) {
  assert_scalar <- function(x, name=deparse(substitute(x))) {
    if (length(x) != 1L) {
      stop(sprintf("%s must be a scalar", name), call. = FALSE)
    }
  }
  assert_scalar(lma_0)
  assert_scalar(B_kl1_1_log10)
  assert_scalar(B_kl1_2_log10)
  assert_scalar(B_kl2_1)
  assert_scalar(B_kl2_2)
  assert_scalar(rho_0)
  assert_scalar(B_dI1)
  assert_scalar(B_dI2)
  assert_scalar(B_ks1)
  assert_scalar(B_ks2)
  assert_scalar(B_rs1)
  assert_scalar(B_rb1)
  assert_scalar(B_f1)
  assert_scalar(narea)
  assert_scalar(narea_0)
  assert_scalar(B_lf1)
  assert_scalar(B_lf2)
  assert_scalar(B_lf3)
  assert_scalar(B_lf4)
  assert_scalar(B_lf5)
  assert_scalar(k_I)
  assert_scalar(latitude)
  assert_scalar(site_prod)

  function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }
    lma       <- with_default("lma")
    rho       <- with_default("rho")
    omega     <- with_default("omega")
    narea     <- with_default("narea", narea)

    ## lma / leaf turnover relationship:
    k_l   <- 10^(B_kl1_1_log10 + site_prod*B_kl1_2_log10)* (lma / lma_0) ^ (-(B_kl2_1+B_kl2_2*site_prod))

    ## rho / mortality relationship:
    d_I  <- B_dI1 * (rho / rho_0) ^ (-B_dI2)

    ## rho / wood turnover relationship:
    k_s  <- B_ks1 *  (rho / rho_0) ^ (-B_ks2)

    ## rho / sapwood respiration relationship:

    ## Respiration rates are per unit mass, so this next line has the
    ## effect of holding constant the respiration rate per unit volume.
    ## So respiration rates per unit mass vary with rho, respiration
    ## rates per unit volume don't.
    r_s <- B_rs1 / rho
    # bark respiration follows from sapwood
    r_b <- B_rb1 / rho

    ## omega / accessory cost relationship
    a_f3 <- B_f1 * omega

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

    ## Respiration rates are per unit mass, so convert to mass-based
    ## rate by dividing with lma
    ## So respiration rates per unit mass vary with lma, while
    ## respiration rates per unit area don't.
    r_l  <- B_lf4 * narea / lma

    extra <- cbind(k_l,                # lma
                   d_I, k_s, r_s, r_b, # rho
                   a_f3,               # omega
                   a_p1, a_p2,         # narea
                   r_l)                # lma, narea

    overlap <- intersect(colnames(m), colnames(extra))
    if (length(overlap) > 0L) {
      stop("Attempt to overwrite generated parameters: ",
           paste(overlap, collapse=", "))
    }
    ## Filter extra so that any column where all numbers are with eps
    ## of the default strategy are not replaced:
    if (filter) {
      if (nrow(extra) == 0L) {
        extra <- NULL
      } else {
        pos <- diff(apply(extra, 2, range)) == 0
        if (any(pos)) {
          eps <- sqrt(.Machine$double.eps)
          x1 <- extra[1, pos]
          x2 <- unlist(s[names(x1)])
          drop <- abs(x1 - x2) < eps & abs(1 - x1/x2) < eps
          if (any(drop)) {
            keep <- setdiff(colnames(extra), names(drop)[drop])
            extra <- extra[, keep, drop=FALSE]
          }
        }
      }
    }

    if (!is.null(extra)) {
      m <- cbind(m, extra)
    }
    m
  }
}



