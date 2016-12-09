
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
                             disturbance_mean_interval=10, data_param = NA) {
require(parallel)
mclapply(vec_site_prod,
         FUN,
         disturbance_mean_interval = disturbance_mean_interval,
         data_param = data_param,
         mc.cores = 7)
}


run_assembly <- function(site_prod=1.0, disturbance_mean_interval=10, ...) {

  p <- trait_gradients_base_parameters(site_prod=site_prod)

  p$disturbance_mean_interval <- disturbance_mean_interval
  sys0 <- community(p, bounds_infinite("lma"),
                     fitness_approximate_control=list(type="gp"))

  obj_m0 <- assembler(sys0, list(birth_move_tol=0.5))
  assembler_run(obj_m0, 20)
}

run_assembly_narea<- function(site_prod=1.0, disturbance_mean_interval=10, ...) {

  p <- trait_gradients_base_parameters(site_prod=site_prod)

  p$disturbance_mean_interval <- disturbance_mean_interval

  bounds_narea <- viable_fitness(bounds(narea=c(1E-5, 1E2)), p, x = 0.01)

# now pass in the bounds -- previously we just passed in infinite bounds
  sys0 <- community(p, bounds_narea,
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
obj_m0 <- assembler(sys0, list(birth_move_tol=0.5, compute_viable_fitness = FALSE))

assembler_run(obj_m0, 20)

}

run_assembly_narea_lma<- function(site_prod=1.0, disturbance_mean_interval=10, ...) {

  p <- trait_gradients_base_parameters(site_prod=site_prod)

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


run_assembly_elev <- function(site_prod=1.0, disturbance_mean_interval=10, data_param) {

  p <- trait_gradients_elev_parameters(data_param = data_param, site_prod = site_prod)

  p$disturbance_mean_interval <- disturbance_mean_interval
  sys0 <- community(p, bounds_infinite("lma"),
                     fitness_approximate_control=list(type="gp"))

  obj_m0 <- assembler(sys0, list(birth_move_tol=0.5))
  assembler_run(obj_m0, 20)
}

run_assembly_slope <- function(site_prod=1.0, disturbance_mean_interval=10, data_param) {

  p <- trait_gradients_slope_parameters(data_param = data_param, site_prod = site_prod)

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
  ## ctrl$equilibrium_solver_name <- "hybrid"
                                 #"hybrid" # in default this is "iteration"
  ctrl$equilibrium_verbose <-  TRUE

  FF16_trait_gradient_hyperpar <- make_FF16_trait_gradient_hyperpar(...)
  p <- FF16_Parameters(patch_area=1.0, control=ctrl,
                   hyperpar=FF16_trait_gradient_hyperpar)
      # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0
  p
}

##' Hopefully sensible set of parameters for use with the EBT.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use. Include a variation of LMA LT elevation of the tradeoff changing with MAP
##' @title Sensible, fast (ish) EBT parameters
##' @author Georges Kunstler
##' @export
trait_gradients_elev_parameters<- function(data_param, ...) {
  #plant_log_console()
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.001
  ctrl$equilibrium_eps <- 1e-5

  ctrl$equilibrium_nsteps  <- 80
  ctrl$equilibrium_solver_name <- "iteration" #"hybrid" # in default this is "iteration"
  ctrl$equilibrium_verbose <-  TRUE

  FF16_trait_gradient_hyperpar <-
      make_FF16_trait_gradient_slope_elev_hyperpar(...,
                                B_kl1_1_log10 = log10(0.4565855),
                                B_kl1_2_log10 = data_param[2, 'LMAelev'] *(2/0.6),
                                B_kl2_1 = 1.71,
                                B_kl2_2 = 0 )
    # this slope is for a gradient of precipitation from -1 to 1 but we can not vary that much the productivity (otherwise growth become zero) so I rescale it to have the same slope on a range only between -0.3 and 0.3 (which is a range where growth is non zero)
  p <- FF16_Parameters(patch_area=1.0, control=ctrl,
                   hyperpar=FF16_trait_gradient_hyperpar)
      # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0
  p
}

##' Hopefully sensible set of parameters for use with the EBT.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use. Include a variation of LMA LT tradeoff changing with MAT over MAP for the slope
##' @title Sensible, fast (ish) EBT parameters
##' @author Georges Kunstler
##' @export
trait_gradients_slope_parameters<- function(data_param, ...) {
  #plant_log_console()
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.001
  ctrl$equilibrium_eps <- 1e-5

  ctrl$equilibrium_nsteps  <- 80
  ctrl$equilibrium_solver_name <- "iteration"
                                 #"hybrid" # in default this is "iteration"
  ctrl$equilibrium_verbose <-  TRUE

  FF16_trait_gradient_hyperpar <-
      make_FF16_trait_gradient_slope_elev_hyperpar(...,
                                B_kl1_1_log10 = log10(0.4565855),
                                B_kl1_2_log10 = 0,
                                B_kl2_1 = 1.71,
                                B_kl2_2 = -data_param[2, 'LMAslope']* (2/0.6) )
    # this slope is for a gradient of Temp / precipitation from -1 to 1 but we can not vary that much the productivity (otherwise growth become zero) so I rescale it to have the same slope on a range only between -0.3 and 0.3 (which is a range where growth is non zero)

  p <- FF16_Parameters(patch_area=1.0, control=ctrl,
                   hyperpar=FF16_trait_gradient_hyperpar)
      # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0
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

demo_ver_siteS_lmaS_lightS<- function(fun_param,
                                      site_prodS = c(0.8, 1, 1.2),
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
        xlab = NA, ylab = NA,
        ylim = range(res[ -1 , var, site_prod_n, , ]),
        col = cols[2])
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

plot_demo_site_lma_light <- function(var, res, site_prodS, lmaS, lightS){
site_prodS_n<- paste0('site_prod', site_prodS)
lmaS_n<- paste0('lma', lmaS)
lightS_n<- paste0('light', lightS)

   for (site_prod_n in site_prodS_n){
       plot_demo_var(res, site_prod_n, var, lmaS_n, lightS_n)

   }
}

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



plot_lcp_version <- function(){
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
plot(dhdt_vec_hs, 1 - lcp_vec,  type = "l", xlab = expression(Height~growth~'in'~full~light~low~water~(m~year^{-1})),
     ylab = '1 - light compensation point (%)', cex.lab = 1.2)
}

