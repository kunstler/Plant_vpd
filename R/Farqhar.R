
#plantecophys
photo_model <- function(){



}

plot_photosyn <-  function(){
require(plantecophys)
plot(1:150, Photosyn(PPFD=1:150)$ALEAF) # pb with light response at low level (check in code with ci set at ca) (OK pb in numerical issue in limitation

##
par(mfrow = c(2,2))
plot(seq(0, 5, length.out = 100),
     Photosyn(VPD = seq(0, 5, length.out = 100), PPFD = 1500,
                          gsmodel = "BBOpti")$ALEAF,
     type = "l")
plot(seq(0, 5, length.out = 100),
     Photosyn(VPD = seq(0, 5, length.out = 100), PPFD = 1500,
                          gsmodel = "BallBerry")$ALEAF,
     type = "l")


plot(seq(0, 10, length.out = 100),
     Photosyn(VPD = seq(0, 10, length.out = 100), PPFD = 1500,
                          gsmodel = "BBOpti")$GS,
     ylim = c(0,1))

plot(seq(0, 10, length.out = 100),
     Photosyn(VPD = seq(0, 10, length.out = 100), PPFD = 1500,
                          gsmodel = "BallBerry")$GS,
     ylim = c(0,1))

#BBOpti
D <- (0:50)/10
V <- 25:125
matA <- matrix(NA, nrow= length(D), ncol = length(V))
for (j in 1:length(V)){
df <- Photosyn(VPD = D,
               Vcmax  = V[j], Jmax = V[j]*1.67,
               gsmodel = "BBOpti")
matA[ ,j] <- df$ALEAF
}
image(D,V,matA, xlab = "VPD", ylab = "Vcmax",
      col = rev(heat.colors(20)))
contour(D,V,matA,add = TRUE)

V_seq <- seq(25, 150, length.out = 5)
cols <- rev(heat.colors(length(V_seq)))
list_df_PPFD_V_VPD1<- vector("list")
list_df_PPFD_V_VPD2.5<- vector("list")
list_df_PPFD_V_VPD5<- vector("list")
for (i in seq_len(length(V_seq))){
list_df_PPFD_V_VPD1[[i]]<- Photosyn(VPD = 1, PPFD = 1:1500,
                          Vcmax  = V_seq[i],
                          Jmax  = V_seq[i]*1.67,
                          gsmodel = "BBOpti")$ALEAF
list_df_PPFD_V_VPD2.5[[i]]<- Photosyn(VPD = 2.5, PPFD = 1:1500,
                          Vcmax  = V_seq[i],
                          Jmax  = V_seq[i]*1.67,
                          gsmodel = "BBOpti")$ALEAF
list_df_PPFD_V_VPD5[[i]]<- Photosyn(VPD = 5, PPFD = 1:1500,
                          Vcmax  = V_seq[i],
                          Jmax  = V_seq[i]*1.67,
                          gsmodel = "BBOpti")$ALEAF
}


df_PPFD_V_VPD1<- as.data.frame(list_df_PPFD_V_VPD1)
df_PPFD_V_VPD2.5<- as.data.frame(list_df_PPFD_V_VPD2.5)
df_PPFD_V_VPD5<- as.data.frame(list_df_PPFD_V_VPD5)
plot(1:1500, df_PPFD_V_VPD1[,1],
     ylim = range(list_df_PPFD_V_VPD1, list_df_PPFD_V_VPD5),
     xlab = expression(paste("PPFD (", mu, "mol", " ",m^-2, " ", s^-1, ")" )),
     ylab = expression(paste(A[n],"(", mu, "mol", " ", m^-2, " ",s^-1, ")" )),
     type = "l", col = cols[1], lty = 3)
for (i in 2:length(V_seq)) lines(1:1500, df_PPFD_V_VPD1[, i], col = cols[i], lty = 3)
for (i in 1:length(V_seq)) lines(1:1500, df_PPFD_V_VPD2.5[, i], col = cols[i], lty = 2)
for (i in 1:length(V_seq)) lines(1:1500, df_PPFD_V_VPD5[, i], col = cols[i], lty = 1)


#BallBerry
D <- (0:30)/10
V <- 25:125
matA <- matrix(NA, nrow= length(D), ncol = length(V))
for (j in 1:length(V)){
df <- Photosyn(VPD = D,
               Vcmax  = V[j], Jmax = V[j]*1.67,
               gsmodel = "BallBerry")
matA[ ,j] <- df$ALEAF
}
image(D,V,matA, xlab = "VPD", ylab = "Vcmax",
      col = rev(heat.colors(40)))
contour(D,V,matA,add = TRUE)

V_seq <- seq(25, 150, length.out = 5)
cols <- rev(heat.colors(length(V_seq)))
list_df_PPFD_V_VPD1<- vector("list")
list_df_PPFD_V_VPD2.5<- vector("list")
list_df_PPFD_V_VPD5<- vector("list")
for (i in seq_len(length(V_seq))){
list_df_PPFD_V_VPD1[[i]]<- Photosyn(VPD = 0.5, PPFD = 1:1500,
                          Vcmax  = V_seq[i],
                          Jmax  = V_seq[i]*1.67,
                          gsmodel = "BallBerry")$ALEAF
list_df_PPFD_V_VPD2.5[[i]]<- Photosyn(VPD = 1, PPFD = 1:1500,
                          Vcmax  = V_seq[i],
                          Jmax  = V_seq[i]*1.67,
                          gsmodel = "BallBerry")$ALEAF
list_df_PPFD_V_VPD5[[i]]<- Photosyn(VPD = 1.5, PPFD = 1:1500,
                          Vcmax  = V_seq[i],
                          Jmax  = V_seq[i]*1.67,
                          gsmodel = "BallBerry")$ALEAF
}


df_PPFD_V_VPD1<- as.data.frame(list_df_PPFD_V_VPD1)
df_PPFD_V_VPD2.5<- as.data.frame(list_df_PPFD_V_VPD2.5)
df_PPFD_V_VPD5<- as.data.frame(list_df_PPFD_V_VPD5)
plot(1:1500, df_PPFD_V_VPD1[,1],
     ylim = range(list_df_PPFD_V_VPD1, list_df_PPFD_V_VPD5),
     xlab = expression(paste("PPFD (", mu, "mol", " ",m^-2, " ", s^-1, ")" )),
     ylab = expression(paste(A[n],"(", mu, "mol", " ", m^-2, " ",s^-1, ")" )),
     type = "l", col = cols[1], lty = 3)
for (i in 2:length(V_seq)) lines(1:1500, df_PPFD_V_VPD1[, i], col = cols[i], lty = 3)
for (i in 1:length(V_seq)) lines(1:1500, df_PPFD_V_VPD2.5[, i], col = cols[i], lty = 2)
for (i in 1:length(V_seq)) lines(1:1500, df_PPFD_V_VPD5[, i], col = cols[i], lty = 1)


###################

V_seq <- seq(25, 150, length.out = 5)
cols <- rev(heat.colors(length(V_seq)))
list_df_V <- vector("list")
for (i in seq_len(length(V_seq))){
list_df_V[[i]]<- Photosyn(VPD = D, Vcmax  = V_seq[i],
                          Jmax  = V_seq[i]*1.67,
                          gsmodel = "BBOpti")$ALEAF
}


df_V <- as.data.frame(list_df_V)
names(df_V) <- paste0("V", V_seq)
plot(D, df_V[,1],
     ylim = range(df_V),
     xlab = "VPD (kPa)",
     ylab = expression(paste(A[n],"(", mu, "mol", m^-2, s^-1, ")" )),
     type = "l", col = cols[1])
for (i in 2:length(V_seq)) lines(D, df_V[, i], col = cols[i])

#TODO plot effect of VPD on function of Daniel with BBopti or BB with EB or not


# why ##' @param k_I light extinction coefficient [dimensionless] in light model ?

                                lma_0=0.1978791
                                B_kl1=0.4565855
                                B_kl2=1.71
                                rho_0=608.0
                                B_dI1=0.01
                                B_dI2=0.0
                                B_ks1=0.2
                                B_ks2=0.0
                                B_rs1=4012.0
                                B_rb1=2.0*4012.0
                                B_f1 =3.0
                                narea=1.87e-3
                                narea_0=1.87e-3
                                B_lf1=5120.738 * 1.87e-3 * 24 * 3600 / 1e+06
                                B_lf2=0.5
                                B_lf3=0.04
                                B_lf4=21000
                                B_lf5=1
                                k_I=0.5
                                latitude=0


library(plantecophys)
df_pred <- Photosyn(PPFD=10:1500,Vcmax  = 75,
                          Jmax  = 75*1.67) # Jmax ~1.67 Vcmax in Medlyn et al. 2002
library(plant)
B_lf1<-5120.738 * 1.87e-3 * 24 * 3600 / 1e+06
B_lf2<-0.5
B_lf3<-0.04
B_lf5<-1
narea<- 1.87e-3
narea_0<-1.87e-3
B_lf1<-5120.738 * 1.87e-3 * 24 * 3600 / 1e+06
Amax <- B_lf1 * (narea/narea_0) ^  B_lf5
theta <- B_lf2
QY <- B_lf3
k_I <-  0.5
assimilation_rectangular_hyperbolae <- function(I, Amax, theta, QY) {
  x <- QY * I + Amax
  (x - sqrt(x^2 - 4 * theta * QY * I * Amax)) / (2 * theta)
}

par(mfrow = c(1,2))
plot(10:1500, df_pred$ALEAF + df_pred$Rd, type = "l",
     xlab = expression(paste("PPFD (", mu, "mol", " ",m^-2, " ", s^-1, ")" )),
     ylab = expression(paste(A,"(", mu, "mol", " ", m^-2, " ",s^-1, ")" )))
plot(seq(0, 216, length.out = 100),
     assimilation_rectangular_hyperbolae(
       seq(0, 216, length.out = 100)*k_I, Amax, theta, QY), type = "l",
     xlab = "mol PAR / m2/ day (ok? and conversion to PPFD?)",
     ylab = expression(paste(A,"(", mu, "mol", " ", m^-2, " ",s^-1, ") or per day ??" )))

## UNIT ??


    ## Photosynthesis  [mol CO2 / m2 / yr]
    approximate_annual_assimilation <- function(narea, latitude) {
      E <- seq(0, 1, by=0.02)
      ## Only integrate over half year, as solar path is symmetrical
      D <- seq(0, 365/2, length.out = 10000)
      I <- plant:::PAR_given_solar_angle(plant:::solar_angle(D, latitude = abs(latitude)))

      Amax <- B_lf1 * (narea/narea_0) ^  B_lf5
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

}
# TODO LINK Vcmax and Jmax (according to Medlyn et al. 2002 Jmax = 1.67 Vcmax)


# Sperry model of stomatal conductance TODO ADD (WHAT full model or just an approximation based on

# cavitation curve
k_xyl <- function(P, b, c, kmax = 5){
kmax * exp(-(P/b)^c)
}

## P <- seq(0, 5, length.out = 100)

## plot(P, k_xyl(P, b = 1, c = 20), type = "l")
## lines(P, k_xyl(P, b = 2, c = 20), lty = 3)
## lines(P, k_xyl(P, b = 2, c = 0.7), lty = 2)

## taken from http://biocycle.atmos.colostate.edu/shiny/photosynthesis/#farquhar.R

## V.max=50
## J.max=100
## APAR=500
## c.i=30
## # Some local parameters we need
## p.sfc <- 101325  # surface air pressure (Pascals)
## gamma <- 3. # CO2 compensation point (Pascals)
## O.i <- 0.209 * p.sfc  # oxygen partial pressure in chloroplast
## K.c <- 30 # Michaelis-Menten constant for carboxylation (Pascals)
## K.o <- 30000 # Michaelis-Menten constant for oxidation (Pascals)

## # Solution of quadratic (Bonan 17.8)
## a <- 0.7
## b <- -(J.max + 0.385 * APAR)
## c <- 0.385 * J.max * APAR
## J.1 <- (-b + sqrt(b^2 - 4 * a * c) ) / (2 * a)
## J.2 <- (-b - sqrt(b^2 - 4 * a * c) ) / (2 * a)
## J <- min(J.1, J.2)
## g <-  0.2


## root2 <- function(V,g=0.20){
## b <- ca-V/g - K
## c <- ca*K+V*gamma/g
## a <- -1
## delta <- b^2 -4*a*c
## return((-b-sqrt(delta))/(2*a))
## }
## V <- 25:125
## plot(V, root2(V), type = "l", ylab = "ci")
## gg <-  (10:70)/100

## plot(gg, root2(V = 50, g = gg), type = "l", ylab = "ci")

## plot(V, V * ( root2(V) - gamma) / ( root2(V)+K) , type = "l", ylab = "A")

## plot(gg, 50 * ( root2(V= 50, g = gg) - gamma) / ( root2(V=50, g = gg)+K) , type = "l", ylab = "A")

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

