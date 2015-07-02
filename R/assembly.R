
lma_gradient_plot <- function(assembly_lma_1, assembly_lma_2, assembly_lma_3) {

  all_data <- list(assembly_lma_1, assembly_lma_2, assembly_lma_3)
  cols <- c("blue", "orange")
  cols <- rev(brewer.pal(5, "Blues"))

  plot(NA, log="x", xlim=c(0.01, 1), ylim=c(-2,0.5),
    xlab= "Leaf mass per area (kg/m2)",
    ylab="Fitness")
  abline(h=0, col="grey")

  for(i in seq_along(all_data)) {
    data <- all_data[[i]]
    ff <- community_fitness_approximate(last(data$history))
    lma <- sort(c(seq_log_range(data$community$bounds, 400), data$community$traits))
    w <- ff(lma)
    lma_res <- data$community$traits[,"lma"]
    points(lma,w, type='l', col=cols[i])
    points(lma_res,ff(lma_res), type='p', col=cols[i], pch=19)
  }
  legend("topleft", bty="n", legend=c("low", "med", "high"), col=cols, pch=19)
}

run_assembly <- function(disturbance_mean_interval=10, site_prod=1.0) {

  p <- trait_gradients_base_parameters(site_prod=site_prod)

  p$disturbance_mean_interval <- disturbance_mean_interval
  sys0 <- community(p, bounds_infinite("lma"))
# sys0 <- community(p, bounds(lma= c(-Inf, Inf), stc=c(0, 100)))

  obj_m0 <- assembler(sys0, list(birth_move_tol=0.05))
  assembler_run(obj_m0, 20)
}

##' Hopefully sensible set of parameters for use with the EBT.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use.
##' @title Sensible, fast (ish) EBT parameters
##' @authorDaniel Falster
##' @export
trait_gradients_base_parameters <- function(...) {
  plant_log_console()
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.005
  ctrl$equilibrium_eps <- 1e-3

  ctrl$equilibrium_nsteps  <- 20
  ctrl$equilibrium_solver_name <- "hybrid"

  FFW16_trait_gradient_hyperpar <- make_FFW16_trait_gradient_hyperpar(...)
  p <- FFW16_Parameters(patch_area=1.0, control=ctrl,
                   hyperpar=FFW16_trait_gradient_hyperpar)

  # neutralise reproduction
  p$strategy_default$c_r1 <- 0.5
  p$strategy_default$c_r2 <- 0
  p
}


##' Hyperparameterisation of FFW16 model used in this analysis.
##' @title Hyperparameters for plant
##' @param lma_0 LMA value...
##' @param k_l_0 ...
##' @param B4 Slope of lma / leaf turnover log-log relationship
##' @param d0_0 Baseline mortality rate
##' @param narea_0 nitrogen per unit leaf area
##' @param c_PN Photosynthesis per unit leaf N
##' @param c_RN Respiration per unit leaf N
##' @export
##' @rdname FFW16_hyperpar
make_FFW16_trait_gradient_hyperpar <- function(
                              site_prod=1.0,
                              B4=1.71,
                              lma_0=0.2,
                              k_l_0=0.45,
                              narea_0=1.87e-3,
                              c_PN = 80000,
                              c_RN = 21000,
                              d0_0=0.01,
                              c_ds = 0,
                              stc_0=0
                              ) {
  force(site_prod)
  force(B4)
  force(lma_0)
  force(k_l_0)
  force(narea_0)
  force(c_PN)
  force(c_RN)
  force(d0_0)
  force(c_ds)
  force(stc_0)

  function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }

    if(nrow(m)==0L) {
      return(m)
    }

    lma       <- with_default("lma")
    narea     <- with_default("narea", narea_0)
    stc       <- with_default("stc", stc_0)

    ## lma / leaf turnover relationship:
    k_l   <- k_l_0 * (lma / lma_0) ^ (-B4)

    ## narea / photosynthesis / respiration
    ## Photosynthesis per mass leaf N [mol CO2 / kgN / yr]
    ## TODO: this is where we improve photosynthesis model
    c_p1  <- c_PN * narea * site_prod

    ## Respiration per mass leaf N [mol CO2 / kgN / yr]
    ## = (6.66e-4 * (365*24*60*60))
    ## Obtained from global average of ratio of dark respiration rate to
    ## leaf nitrogen content using the GLOPNET dataset
    ## Respiration rates are per unit mass, so convert to mass-based
    ## rate by dividing with lma
    ## So respiration rates per unit mass vary with lma, while
    ## respiration rates per unit area don't.
    c_Rl  <- (c_RN * narea) / lma

    ## Stress tolerance trait -- reduces mortality for given
    ## respiration cost per unit leaf area

    c_Rl  <- c_Rl + (stc / lma)

    ## Baseline mortality rate
    c_d0  <- d0_0
    ## Decreases with investment in stress tolerance
    ## By default c_ds = 0 so trait has no effect
    c_d0  <- d0_0 * exp(-c_ds * stc)

    extra <- cbind(k_l,                   # lma
                   c_d0,                  # stress_tolerance
                   c_p1, c_Rl)            # narea

    overlap <- intersect(colnames(m), colnames(extra))
    if (length(overlap) > 0L) {
      stop("Attempt to overwrite generated parameters: ",
           paste(overlap, collapse=", "))
    }

    ## Filter extra so that any column where all numbers are with eps
    ## of the default strategy are not replaced:
    if (filter) {
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
    cbind(m, extra)
  }
}
