## taken from http://biocycle.atmos.colostate.edu/shiny/photosynthesis/#farquhar.R

farquhar <- function(V.max=50, J.max=100, APAR=500, c.i=30) {

# Model inputs:
  # V.max = maximum rubisco-limited rate in micromoles per (m^2 sec)
  # J.max = maximum light-limited rate in micromoles per (m^2 sec)
  # APAR = absorbed photosynthetically-active radiation in micromoles per (m^2 sec)
  # c.i = intercellular CO2 partial pressure in Pascals (roughly ppm/10)

# Some local parameters we need
p.sfc <- 101325  # surface air pressure (Pascals)
gamma <- 3. # CO2 compensation point (Pascals)
O.i <- 0.209 * p.sfc  # oxygen partial pressure in chloroplast
K.c <- 30 # Michaelis-Menten constant for carboxylation (Pascals)
K.o <- 30000 # Michaelis-Menten constant for oxidation (Pascals)

# Solution of quadratic (Bonan 17.8)
a <- 0.7
b <- -(J.max + 0.385 * APAR)
c <- 0.385 * J.max * APAR
J.1 <- (-b + sqrt(b^2 - 4 * a * c) ) / (2 * a)
J.2 <- (-b - sqrt(b^2 - 4 * a * c) ) / (2 * a)
J <- min(J.1, J.2)

# Rubisco-limited rate of photosynthesis
w.c <- V.max * (c.i - gamma) / (c.i + K.c * (1 + O.i/K.o))  # Bonan 17.6

# Light-limited rate of photosynthesis
w.j <- J * (c.i - gamma) / (4 * (c.i + 2 * gamma))            # Bonan 17.7

# Sink-limited rate of photosynthesis
w.s <- V.max / 2

# Dark respiration
R.d <- 0.015 * V.max

# Net assimilation
A.n <- min(w.c, w.j, w.s)-R.d

return(A.n)
}


 ## Photosynthesis  [mol CO2 / m2 / yr]
    approximate_annual_assimilation <- function(narea, latitude) {
      E <- seq(0, 1, by=0.02)
      ## Only integrate over half year, as solar path is symmetrical
      D <- seq(0, 365/2, length.out = 10000)
      I <- PAR_given_solar_angle(solar_angle(D, latitude = abs(latitude)))

      Amax <- B_lf1 * (narea/narea_0) ^  B_lf5
      theta <- B_lf2
      QY <- B_lf3

      AA <- NA * E

      for (i in seq_len(length(E))) {
        AA[i] <- 2 * trapezium(D, assimilation_rectangular_hyperbolae(
                                    k_I * I * E[i], Amax, theta, QY))
        ## FUNCTION IN Cpp in inst/include/plant/util.h
      }
      if(all(diff(AA) < 1E-8)) {
        # line fitting will fail if all have are zero, or potentially same value
        ret <- c(last(AA), 0)
        names(ret) <- c("p1","p2")
      } else {
        fit <- nls(AA ~ p1 * E/(p2 + E), data.frame(E = E, AA = AA), start = list(p1 = 100, p2 = 0.2))
        ret <- coef(fit)
      }
      ret
    }

# TODO CHANGE 'assimilation_rectangular_hyperbolae' to Farqhuar with water stress

assimilation_rectangular_hyperbolae <- function(I, Amax, theta, QY) {
      x <- QY * I + Amax
      (x - sqrt(x^2 - 4 * theta * QY * I * Amax)) / (2 * theta)
    }
'
