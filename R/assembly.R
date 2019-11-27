
run_assembly_FvCB_narea<- function(vpd=0.0, disturbance_mean_interval=10, ...) {
  p <- trait_gradients_FvCB_parameters(vpd=vpd)

  p$disturbance_mean_interval <- disturbance_mean_interval
  # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0

  sys0 <- community(p, rbind(bounds(narea=c(1E-4, 1E-2))),
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
obj_m0 <- assembler(sys0, list(birth_move_tol=0.5, compute_viable_fitness = FALSE))
print('initialized')
assembler_run(obj_m0, 20)
}

run_assembly_FvCB_lma<- function(vpd=0.0, disturbance_mean_interval=10, ...) {
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
  Tleaf <- coefs[1] + coefs[2]*log(vpd + coefs[3])
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



run_assembly_FvCB_narea_lma_Narea_LTR<- function(vpd=0.0,
                                                 disturbance_mean_interval=10,
                                                 ...) {

 p <- trait_gradients_FvCB_LTR_Narea_parameters(vpd=vpd)

  p$disturbance_mean_interval <- disturbance_mean_interval
  # neutralise reproduction
  p$strategy_default$a_f1 <- 0.5
  p$strategy_default$a_f2 <- 0

    sys0 <- community(p, rbind(bounds(narea=c(1E-4, 1E-2)),
                               bounds(lma=c(0.001, 3))),
                   fitness_approximate_control=list(type="gp"))

# and also tell the assembler not to calculate the bounds
    obj_m0 <- assembler(sys0, list(birth_move_tol=0.5,
                                   compute_viable_fitness = FALSE))
print('initialized')
assembler_run(obj_m0, 20)
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


##' Hopefully sensible set of parameters for use with the EBT.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use.
##' @title Sensible, fast (ish) EBT parameters
##' @authorDaniel Falster
##' @export
trait_gradients_FvCB_LTR_Narea_parameters <- function(...) {
  #plant_log_console()
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.001
  ctrl$equilibrium_eps <- 1e-5
  # TO TRY TO SOLVE ERROR MESSAGE
  ctrl$ode_step_size_min <- 1e-8

  ctrl$equilibrium_nsteps  <- 80
  ctrl$equilibrium_solver_name <- "iteration"
  ## ctrl$equilibrium_solver_name <- "hybrid"
                                 #"hybrid" # in default this is "iteration"
  ctrl$equilibrium_verbose <-  TRUE

    FF16_FvCB_trait_gradient_hyperpar <-
        make_FF16_FvCB_trait_gradient_LTRvsNarea_hyperpar(...)
  p <- FF16FvCB_Parameters(patch_area=1.0, control=ctrl,
                   hyperpar=FF16_FvCB_trait_gradient_hyperpar)
  p
}


#### FF16 FvCB with Narea effect on LTR

make_FF16_FvCB_trait_gradient_LTRvsNarea_hyperpar <- function(lma_0=0.1978791,
                                      B_kl1_1 = -0.1420071,
                                      B_kl1_2 =  0.5213681,
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
                                      B_lf1= 31.62 *1000^0.801, # http://doi.wiley.com/10.1111/gcb.12870 conversion of Nare in g m-2
                                      B_lf2=0.7,
                                      B_lf3=0.3,
                                      B_lf4=21000,
                                      B_lf5=0.801, # http://doi.wiley.com/10.1111/gcb.12870
                                      B_lf6=1.67, # Medlyn et al. 2002 / for Walker 1.01
                                      B_lf7 = 1, # Medlyn et al. 2002 / for Walker 0.89
                                      k_I=0.5, latitude=0, vpd = 0, Tleaf= 25) {
  assert_scalar <- function(x, name=deparse(substitute(x))) {
    if (length(x) != 1L) {
      stop(sprintf("%s must be a scalar", name), call. = FALSE)
    }
  }
  assert_scalar(lma_0)
  assert_scalar(B_kl1_1)
  assert_scalar(B_kl1_2)
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
  assert_scalar(B_lf6)
  assert_scalar(B_lf7)
  assert_scalar(k_I)
  assert_scalar(latitude)
  assert_scalar(vpd)
  assert_scalar(Tleaf)
  require(plantecophys)
  require(nlmrt)
  ## TODO: k_I should actually be in default parameter set, so perhaps don't pass into function?

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
   k_l   <- (B_kl1_1 + narea/narea_0*B_kl1_2)* (lma / lma_0) ^ (-(B_kl2_1+B_kl2_2*narea/narea_0))

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

    assimilation_FvCB <- function(I, V, vpd, Tleaf, alpha, theta, lf6, lf7) {
    df_pred <- plantecophys::Photosyn(PPFD=I *1e+06/(24*3600), # conversion of light in mu mol /m2 /s is it ok ?
                        VPD = vpd,
                        Tleaf = Tleaf,
                        Vcmax  = V,
                        Jmax  = lf6*V^lf7, # Jmax ~1.67 Vcmax in Medlyn et al. 20
                        alpha = alpha,
                        theta = theta,
                        gsmodel = "BBLeuning")
    return((df_pred$ALEAF + df_pred$Rd)*24 * 3600 / 1e+06)
    }


    ## Photosynthesis  [mol CO2 / m2 / yr]
    approximate_annual_assimilation <- function(narea, latitude, vpd, Tleaf) {
      E <- seq(0, 1, by=0.02)
      ## Only integrate over half year, as solar path is symmetrical
      D <- seq(0, 365/2, length.out = 10000)
      I <- plant:::PAR_given_solar_angle(
                       plant:::solar_angle(D, latitude = abs(latitude)))

      Vcmax <- B_lf1 * (narea) ^  B_lf5
      theta <- B_lf2
      alpha <- B_lf3
      lf6   <- B_lf6
      lf7   <- B_lf7
      AA <- NA * E

      for (i in seq_len(length(E))) {
        AA[i] <- 2 * plant:::trapezium(D, assimilation_FvCB(
                                    k_I * I * E[i], Vcmax, vpd, Tleaf, alpha, theta, lf6, lf7))
      }
      if(all(diff(AA) < 1E-8)) {
        # line fitting will fail if all have are zero, or potentially same value
        ret <- c(last(AA), 0)
        names(ret) <- c("p1","p2", "p3")
      } else {
        fitxb2 <- nlmrt::nlxb(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3),
                           data = data.frame(E = E, AA = AA),
                           start = list(p1 = 500,
                                        p2 = 600,
                                        p3 = 0.8),
                           lower = 0.01,
                           upper = c(7000, 1200, 1))
        ret <- fitxb2$coefficients
        names(ret) <- c("p1","p2", "p3")
      }
      ret
    }

    # This needed in case narea has length zero, in which case trapezium fails
    a_p1 <- a_p2 <- a_p3 <- 0 * narea
    ## TODO: Remove the 0.5 hardcoded default for k_I here, and deal
    ## with this more nicely.
    if (length(narea) > 0 || k_I != 0.5) {
      i <- match(narea, unique(narea))
      y <- vapply(unique(narea), approximate_annual_assimilation,
                  numeric(3), latitude, vpd, Tleaf)
      a_p1  <- y["p1", i]
      a_p2  <- y["p2", i]
      a_p3  <- y["p3", i]
    }

    ## Respiration rates are per unit mass, so convert to mass-based
    ## rate by dividing with lma
    ## So respiration rates per unit mass vary with lma, while
    ## respiration rates per unit area don't.
    r_l  <- B_lf4 * narea / lma

    extra <- cbind(k_l,                # lma
                   d_I, k_s, r_s, r_b, # rho
                   a_f3,               # omega
                   a_p1, a_p2, a_p2,   # narea
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
