#Look at use of FvC model with plantecophys in Plant


##############################
### LOOK AT CONDUCTANCE MODEL

## BBOpti BallBerry BBLeuning but would be good to add one based on soil water potential as in Sperry or in Sterck


## Plot coupled photosynthesis model and stomatal conductance model
plot_FvC_stomatal <- function(gs_type = "BBLeuning"){
require(plantecophys)
mar.default <- c(5,4,4,2) + 0.1
par(mfrow = c(1, 2), mar = mar.default + c(0, 0.3, 0, 0))
D <- (0:50)/10
V <- 25:125
matA <- matrix(NA, nrow= length(D), ncol = length(V))
for (j in 1:length(V)){
  df <- Photosyn(VPD = D,
                 Vcmax  = V[j], Jmax = V[j]*1.67,
                 gsmodel = gs_type)
  matA[ ,j] <- df$ALEAF
}
image(D,V,matA, xlab = "VPD", ylab = "Vcmax",
      col = rev(heat.colors(20)), main = gs_type)
contour(D,V,matA,add = TRUE)

V_seq <- seq(25, 150, length.out = 4)
colfunc <- colorRampPalette(c("grey", "red"))
cols <- colfunc(length(V_seq))

list_df_PPFD_V_VPD0<- vector("list")
list_df_PPFD_V_VPD1.5<- vector("list")
for (i in seq_len(length(V_seq))){
  list_df_PPFD_V_VPD0[[i]]<- Photosyn(VPD = 0, PPFD = 1:1500,
                            Vcmax  = V_seq[i],
                            Jmax  = V_seq[i]*1.67,
                            gsmodel = gs_type)$ALEAF
  list_df_PPFD_V_VPD1.5[[i]]<- Photosyn(VPD = 1.5, PPFD = 1:1500,
                            Vcmax  = V_seq[i],
                            Jmax  = V_seq[i]*1.67,
                            gsmodel = gs_type)$ALEAF
}

df_PPFD_V_VPD0<- as.data.frame(list_df_PPFD_V_VPD0)
df_PPFD_V_VPD1.5<- as.data.frame(list_df_PPFD_V_VPD1.5)
plot(1:1500, df_PPFD_V_VPD0[,1],
     ylim = range(list_df_PPFD_V_VPD0, list_df_PPFD_V_VPD1.5),
     xlab = expression(paste("PPFD (", mu, "mol", " ",m^-2, " ", s^-1, ")" )),
     ylab = expression(paste(A[n],"(", mu, "mol", " ", m^-2, " ",s^-1, ")" )),
     type = "l", col = cols[1], lty = 1)
for (i in 2:length(V_seq)) lines(1:1500, df_PPFD_V_VPD0[, i], col = cols[i], lty = 1, lwd = 2)
for (i in 1:length(V_seq)) lines(1:1500, df_PPFD_V_VPD1.5[, i], col = cols[i], lty = 2, lwd = 2)
legend(x = 0, y = max(df_PPFD_V_VPD0)*0.9, lwd = 2, lty = c(1, 2, rep(1,4)),
       col = c("black", "black", cols),
       legend = c("VPD = 0", "VPD = 1.5", paste("Vcmax", round(V_seq))),
       bty = "n", cex = 0.7)
}

VPD_response <- function(V = 50){
  require(plantecophys)
mar.default <- c(5,4,4,2) + 0.1
par(mfrow = c(1, 2), mar = mar.default + c(0, 0.3, 0, 0))
 plot(seq(0, 5, length.out = 100),
      Photosyn(VPD = seq(0, 5, length.out = 100), PPFD = 1500,
               gsmodel = "BBOpti",
               Vcmax  = V, Jmax  = V*1.67)$ALEAF,
      type = "l", ylim = c(0, 13),
      xlab = expression(paste("VPD (kPa)")),
      ylab = expression(paste(A,"(", mu, "mol", " ", m^-2, " ",s^-1, ")" )))
 lines(seq(0, 5, length.out = 100),
      Photosyn(VPD = seq(0, 5, length.out = 100), PPFD = 1500,
               Vcmax  = V, Jmax  = V*1.67,
               gsmodel = "BallBerry")$ALEAF,
      lty = 2)
 lines(seq(0, 5, length.out = 100),
      Photosyn(VPD = seq(0, 5, length.out = 100), PPFD = 1500,
               Vcmax  = V, Jmax  = V*1.67,
               gsmodel = "BBLeuning")$ALEAF,
      lty = 3)
 #
 plot(seq(0, 10, length.out = 100),
      Photosyn(VPD = seq(0, 10, length.out = 100), PPFD = 1500,
               Vcmax  = V, Jmax  = V*1.67,
               gsmodel = "BBOpti")$GS,
      ylim = c(0,1), type = "l",
      xlab = expression(paste("VPD (kPa)")),
      ylab = expression(paste(Gs," (mol", " ", m^-2, " ",s^-1, ")" )))
 lines(seq(0, 10, length.out = 100),
      Photosyn(VPD = seq(0, 10, length.out = 100), PPFD = 1500,
               Vcmax  = V, Jmax  = V*1.67,
               gsmodel = "BallBerry")$GS,
      lty = 2)
 lines(seq(0, 10, length.out = 100),
      Photosyn(VPD = seq(0, 10, length.out = 100), PPFD = 1500,
               Vcmax  = V, Jmax  = V*1.67,
               gsmodel = "BBLeuning")$GS,
      lty = 3)
legend(x = 0.1, y = 0.8, lty = 1:3,
       legend = c("BBOpti", "BallBerry","BBLeuning"),
       bty = "n", cex = 0.8)
}

## compare original Plant model vs Farquhar
Compare_Photo_Plant_FvC <- function(){
require(plantecophys)
require(plant)
x <- 10:1500
B_lf1<-5120.738 * 1.87e-3 * 24 * 3600 / 1e+06
B_lf2<- 0.5
B_lf3<- 0.04
B_lf5 <- 1
k_I <-  0.5
# Potential CO_2 photosynthesis at average leaf nitrogen [mol / d / m2]
# For comparison, convert Amax to ymol / m2 /s
Amax <- B_lf1 * 1e6 /(24*3600)

assimilation_rectangular_hyperbolae <- function(I, Amax, theta, QY) {
  x <- QY * I + Amax
  (x - sqrt(x^2 - 4 * theta * QY * I * Amax)) / (2 * theta)
}
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 0.3, 0, 0))

plot(x,
     assimilation_rectangular_hyperbolae(
       x*k_I, Amax, theta = B_lf2, QY = B_lf3),
     type = "l", ylim = c(0, 13), lwd = 2,
     xlab = expression(paste("PPFD (", mu, "mol", " ",m^-2, " ", s^-1, ")" )),
     ylab = expression(paste(A,"(", mu, "mol", " ", m^-2, " ",s^-1, ")" )))
V <-  50
df_pred <- Photosyn(PPFD=x,Vcmax  = 50,
                    Jmax  = 100,
                    alpha = 0.24, theta = 0.85)
lines(x,df_pred$ALEAF + df_pred$Rd,
      lty = 1, col = "red", lwd = 2)
V <-  50
df_pred <- Photosyn(PPFD=x,Vcmax  = V,
                    Jmax  = V*1.67,# Jmax ~1.67 Vcmax in Medlyn et al. 2002
                    alpha = 0.24, theta = 0.85)
lines(x,df_pred$ALEAF + df_pred$Rd,
      lty = 2, col = "red", lwd = 2)
V <-  32.5
df_pred <- Photosyn(PPFD=x,Vcmax  = V,
                    Jmax  = V*1.67,# Jmax ~1.67 Vcmax in Medlyn et al. 2002
                    alpha = 0.04*4, theta = 0.5)
lines(x,df_pred$ALEAF + df_pred$Rd,
      lty = 3, col = "red", lwd = 2)
df_pred <- Photosyn(PPFD=x,Vcmax  = V,
                    Jmax  = 1.67*V,#exp(1.197)*V^0.8,
                    alpha = 0.075*4, theta = 0.7) # Jmax ~1.67 Vcmax in Medlyn et al. 2002
lines(x,df_pred$ALEAF + df_pred$Rd,
      lty = 1, col = "blue", lwd = 2)
df_pred <- Photosyn(PPFD=x,Vcmax  = V,
                    Jmax  = 1.67*V,#exp(1.197)*V^0.8,
                    alpha = 0.25, theta = 0.5) # Jmax ~1.67 Vcmax in Medlyn et al. 2002
lines(x,df_pred$ALEAF + df_pred$Rd,
      lty = 1, col = "green", lwd = 2)
legend(x = 400, y = max(assimilation_rectangular_hyperbolae(
       x*k_I, Amax, theta = B_lf2, QY = B_lf3))*0.6,
       lwd = 2, lty = c(1, 1, 2, 3, 1, 1), cex = 0.8,
       col = c("black", "red", "red", "red", "blue", "green"),
       legend = c("Plant", "FvC_plantecophys", "FvC_plantecophys V = 50 and J = 1.67 V",
                  "FvC Plant and V =32.5 and J = 1.67 V", "FvC Troll = von Caemmerer 2000", "Sterck et al. 2011" ),
       bty = "n")
}

Photo_Plant_match_FvC <- function(){
require(plantecophys)
require(plant)
x <- 10:1500
B_lf1<-5120.738 * 1.87e-3 * 24 * 3600 / 1e+06*1.4
B_lf2<-0.85 # too match the shape of the default FvC initialy 0.5
B_lf3<-0.083 # too match the shape of the default FvC initialy 0.04
B_lf5 <- 1
k_I <-  0.5
# Potential CO_2 photosynthesis at average leaf nitrogen [mol / d / m2]
# For comparison, convert Amax to ymol / m2 /s
Amax <- B_lf1 * 1e6 /(24*3600)

assimilation_rectangular_hyperbolae <- function(I, Amax, theta, QY) {
  x <- QY * I + Amax
  (x - sqrt(x^2 - 4 * theta * QY * I * Amax)) / (2 * theta)
}

plot(x,
     assimilation_rectangular_hyperbolae(
       x*k_I, Amax, theta = B_lf2, QY = B_lf3),
     type = "l", ylim = c(0, 13), lwd = 2,
     xlab = expression(paste("PPFD (", mu, "mol", " ",m^-2, " ", s^-1, ")" )),
     ylab = expression(paste(A,"(", mu, "mol", " ", m^-2, " ",s^-1, ")" )))

V <- 50
df_pred <- Photosyn(PPFD=x,Vcmax  = V,
                    Jmax  = 1.67*V,#exp(1.197)*V^0.8,
                    alpha = 0.24, theta = 0.85) # Jmax ~1.67 Vcmax in Medlyn et al. 2002
lines(x,df_pred$ALEAF + df_pred$Rd,
      lty = 1, col = "red", lwd = 2)
legend(x = 1100, y = 2, lwd = c(2,2),
       col = c("black", "red"),
       legend = c("Plant", "FvC_plantecophys"),
       bty = "n")
}



plot_photosyn_annual_plant<-  function(){
require(plantecophys)
require(plant)
#######################
## fit annual A model

B_lf1<-5120.738 * 1.87e-3 * 24 * 3600 / 1e+06
B_lf2<- 0.5
B_lf3<- 0.04
B_lf5 <- 1
k_I <-  0.5

narea<- 1.87e-3
narea_0<-1.87e-3
B_lf1<-5120.738 * 1.87e-3 * 24 * 3600 / 1e+06
Amax <- B_lf1 * (narea/narea_0) ^  B_lf5
theta <- B_lf2
QY <- B_lf3
k_I <-  0.5
latitude <- 0
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
      # TODO need to look at the fit
    fit <- nls(AA ~ p1 * E/(p2 + E), data.frame(E = E, AA = AA),
               start = list(p1 = 100, p2 = 0.2))
    ret <- coef(fit)
  }
  return(list(param = ret, df = data.frame(E = E, AA = AA)))
}

y <- approximate_annual_assimilation(narea, latitude)
a_p1  <- y$param["p1"]
a_p2  <- y$param["p2"]

mar.default <- c(5,4,4,2) + 0.1
par(mfrow = c(1, 2), mar = mar.default + c(0, 0.3, 0, 0))
plot(0:200, assimilation_rectangular_hyperbolae(0:200, Amax, theta, QY),
     type = "l",
     xlab = expression(paste("PPFD (mol", " ",m^-2, " ", d^-1, ")" )),
     ylab = expression(paste(A,"(mol", " ", m^-2, " ",d^-1, ")" )))

E <- seq(0, 1, by=0.02)
plot(y$df$E, y$df$AA,
     xlab = "Percentage full light",
     ylab = expression(paste(A,"(mol", " ", m^-2, " ",y^-1, ")" )))
lines(E, a_p1*E/(a_p2+E), col = "red")
#exponential model
# see https://dl.sciencesocieties.org/publications/aj/pdfs/107/2/786
# https://www.csusm.edu/terl/publications1/Lobo_13_Photo.pdf
# http://www.sciencedirect.com/science/article/pii/S002251938080004X
fit <- nls(AA ~ p1 *(1-exp(-p2* E/p1)), y$df,
           start = list(p1 = 120, p2 = 1000))
lines(E, coef(fit)[["p1"]]*(1-exp(-E*coef(fit)[["p2"]]/coef(fit)[["p1"]])),
      col = "green", lwd = 2)
# see http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2435.2009.01630.x/epdf
# non-rectangular hyperbola
fit <- nls(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3), y$df,
           start = list(p1 = 120, p2 = 400, p3 = 0.9),
           lower=c(10,2, 0.2), upper=c(200,800, 1),
           algorithm = "port")
lines(E, (coef(fit)[["p1"]] +coef(fit)[["p2"]]*E -
          sqrt((coef(fit)[["p1"]]+coef(fit)[["p2"]]*E)^2-
               4*coef(fit)[["p3"]]*coef(fit)[["p2"]]*E*coef(fit)[["p1"]]))/
         (2*coef(fit)[["p3"]]),
      col = "blue", lwd = 2)
legend(x = 0.2, y = max(y$df$AA)*0.3, lwd = c(NA,2, 2, 2), pch = c(1, NA, NA, NA),
       col = c("black", "red", "green", "blue"),
       legend = c("Data Integrated", "Michealis Menten", "Exponential", "Non Rectangular Hyperbola"),
       bty = "n", cex = 0.75)

}

### FvC

plot_photosyn_annual_FvC<-  function(vpd = 0, alpha = 0.3, theta = 0.7, Vcmax = 38){
require(plantecophys)
require(plant)
k_I <-  0.5
latitude <- 0
assimilation_curve <- function(I, V, vpd, alpha, theta) {
df_pred <- Photosyn(PPFD=I *1e+06/(24*3600), # conversion of light in mu mol /m2 /s is it ok ?
                    VPD = vpd,
                    Vcmax  = V,
                    Jmax  = V*1.67,
                    alpha = alpha,
                    theta = theta,
                    gsmodel = "BBLeuning") # Jmax ~1.67 Vcmax in Medlyn et al. 20
return((df_pred$ALEAF + df_pred$Rd)*24 * 3600 / 1e+06)
}

## Photosynthesis  [mol CO2 / m2 / yr]
approximate_annual_assimilation <- function(V, latitude, vpd, alpha, theta) {
  E <- seq(0, 1, by=0.02)
  ## Only integrate over half year, as solar path is symmetrical
  D <- seq(0, 365/2, length.out = 10000)
  I <- plant:::PAR_given_solar_angle(plant:::solar_angle(D, latitude = abs(latitude)))
  AA <- NA * E
# TODO
# add Vcmax dependence on Narea
  for (i in seq_len(length(E))) {
    AA[i] <- 2 * plant:::trapezium(D, assimilation_curve(
                                k_I * I * E[i], V, vpd, alpha, theta))
  }
  if(all(diff(AA) < 1E-8)) {
    # line fitting will fail if all have are zero, or potentially same value
    ret <- c(last(AA), 0)
    names(ret) <- c("p1","p2")
  } else {
    fit <- nls(AA ~ p1 * E/(p2 + E), data.frame(E = E, AA = AA),
               start = list(p1 = 160, p2 = 0.1))
    ret <- coef(fit)
  }
  return(data.frame(E = E, AA = AA))
}

data_a<- approximate_annual_assimilation(V = Vcmax, latitude, vpd, alpha, theta)

fit <- nls(AA ~ p1 * E/(p2 + E), data_a,
           start = list(p1 = 160, p2 = 0.1))
a_p1  <- coef(fit)[["p1"]]
a_p2  <- coef(fit)[["p2"]]

mar.default <- c(5,4,4,2) + 0.1
par(mfrow = c(1, 2), mar = mar.default + c(0, 0.3, 0, 0))
plot(0:200, assimilation_curve(0:200,  V = 38, vpd, alpha , theta),
     type = "l",
     xlab = expression(paste("PPFD (mol", " ",m^-2, " ", d^-1, ")" )),
     ylab = expression(paste(A,"(mol", " ", m^-2, " ",d^-1, ")" )))

E <- seq(0, 1, by=0.02)
plot(data_a$E, data_a$AA,
     xlab = "percentage full light",
     ylab = expression(paste(A,"(mol", " ", m^-2, " ",y^-1, ")" )))

# Michaelis menten
fit <- nls(AA ~ p1 * E/(p2 + E), data_a,
           start = list(p1 = 160, p2 = 0.1))
lines(E, coef(fit)[["p1"]]*E/(coef(fit)[["p2"]]+E), col = "red", lwd = 2)
#exponential model
# see https://dl.sciencesocieties.org/publications/aj/pdfs/107/2/786
# https://www.csusm.edu/terl/publications1/Lobo_13_Photo.pdf
# http://www.sciencedirect.com/science/article/pii/S002251938080004X
fit <- nls(AA ~ p1 *(1-exp(-p2* E/p1)), data_a,
           start = list(p1 = 120, p2 = 1000))
lines(E, coef(fit)[["p1"]]*(1-exp(-E*coef(fit)[["p2"]]/coef(fit)[["p1"]])),
      col = "green", lwd = 2)
# see http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2435.2009.01630.x/epdf
# non-rectangular hyperbola
fit <- nls(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3), data_a,
           start = list(p1 = 120, p2 = 400, p3 = 0.9),
           lower=c(10,2, 0.2), upper=c(200,800, 1), algorithm = "port")
lines(E, (coef(fit)[["p1"]] +coef(fit)[["p2"]]*E -
          sqrt((coef(fit)[["p1"]]+coef(fit)[["p2"]]*E)^2-
               4*coef(fit)[["p3"]]*coef(fit)[["p2"]]*E*coef(fit)[["p1"]]))/
         (2*coef(fit)[["p3"]]),
      col = "blue", lwd = 2)
legend(x = 0.2, y = max(data_a$AA)*0.3, lwd = c(NA,2, 2, 2), pch = c(1, NA, NA, NA),
       col = c("black", "red", "green", "blue"),
       legend = c("Data Integrated", "Michealis Menten", "Exponential", "Non Rectangular Hyperbola"),
       bty = "n", cex = 0.75)
}






plot_photosyn_annual_FvC_Narea<-  function(vpd = 0,
                                              narea = 1.87e-3,
                                              B_lf1= 31.62 *1000^0.801,
                                          # HTTP://DOI.WILEY.COM/10.1111/GCB.12870 CONVERSION OF NARE IN G M-2
                                              B_lf2=0.7,
                                              B_lf3=0.3,
                                              B_lf4=21000,
                                              B_lf5=0.801, # http://doi.wiley.com/10.1111/gcb.12870
                                              B_lf6=1.67,
                                              k_I=0.5,
                                              latitude=0.0){
require(plantecophys)
require(plant)

    assimilation_FvCB <- function(I, V, vpd, alpha, theta, lf6) {
    df_pred <- Photosyn(PPFD=I *1e+06/(24*3600), # conversion of light in mu mol /m2 /s is it ok ?
                        VPD = vpd,
                        Vcmax  = V,
                        Jmax  = V*lf6,
                        alpha = alpha,
                        theta = theta,
                        gsmodel = "BBLeuning") # Jmax ~1.67 Vcmax in Medlyn et al. 20
    return((df_pred$ALEAF + df_pred$Rd)*24 * 3600 / 1e+06)
    }

    ## Photosynthesis  [mol CO2 / m2 / yr]
    approximate_annual_assimilation <- function(narea, latitude, vpd) {
      E <- seq(0, 1, by=0.02)
      ## Only integrate over half year, as solar path is symmetrical
      D <- seq(0, 365/2, length.out = 10000)
      I <- plant:::PAR_given_solar_angle(plant:::solar_angle(D, latitude = abs(latitude)))

      Vcmax <- B_lf1 * (narea) ^  B_lf5
      theta <- B_lf2
      alpha <- B_lf3
      lf6   <- B_lf6
      AA <- NA * E

      for (i in seq_len(length(E))) {
        AA[i] <- 2 * plant:::trapezium(D, assimilation_FvCB(
                                    k_I * I * E[i], Vcmax, vpd, alpha, theta, lf6))
      }
        return(data.frame(E = E, AA = AA))
    }

data_a<- approximate_annual_assimilation(narea, latitude, vpd)

fit <- nls(AA ~ p1 * E/(p2 + E), data_a,
           start = list(p1 = 160, p2 = 0.1))
a_p1  <- coef(fit)[["p1"]]
a_p2  <- coef(fit)[["p2"]]

mar.default <- c(5,4,4,2) + 0.1
par(mfrow = c(1, 2), mar = mar.default + c(0, 0.3, 0, 0))
theta <- B_lf2
alpha <- B_lf3
plot(0:200, assimilation_FvCB(0:200, B_lf1 * (narea) ^  B_lf5, vpd, alpha , theta, B_lf6),
     type = "l",
     xlab = expression(paste("PPFD (mol", " ",m^-2, " ", d^-1, ")" )),
     ylab = expression(paste(A,"(mol", " ", m^-2, " ",d^-1, ")" )))

E <- seq(0, 1, by=0.02)
plot(data_a$E, data_a$AA,
     xlab = "percentage full light",
     ylab = expression(paste(A,"(mol", " ", m^-2, " ",y^-1, ")" )))

# Michaelis menten
fit <- nls(AA ~ p1 * E/(p2 + E), data_a,
           start = list(p1 = 160, p2 = 0.1))
lines(E, coef(fit)[["p1"]]*E/(coef(fit)[["p2"]]+E), col = "red", lwd = 2)
#exponential model
# see https://dl.sciencesocieties.org/publications/aj/pdfs/107/2/786
# https://www.csusm.edu/terl/publications1/Lobo_13_Photo.pdf
# http://www.sciencedirect.com/science/article/pii/S002251938080004X
fit <- nls(AA ~ p1 *(1-exp(-p2* E/p1)), data_a,
           start = list(p1 = 120, p2 = 1000), algorithm = "port")

lines(E, coef(fit)[["p1"]]*(1-exp(-E*coef(fit)[["p2"]]/coef(fit)[["p1"]])),
      col = "green", lwd = 2)
# see http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2435.2009.01630.x/epdf
# non-rectangular hyperbola
fit <- nls(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3), data_a,
           start = list(p1 = 170, p2 = 600, p3 = 0.9))

lines(E, (coef(fit)[["p1"]] +coef(fit)[["p2"]]*E -
          sqrt((coef(fit)[["p1"]]+coef(fit)[["p2"]]*E)^2-
               4*coef(fit)[["p3"]]*coef(fit)[["p2"]]*E*coef(fit)[["p1"]]))/
         (2*coef(fit)[["p3"]]),
      col = "blue", lwd = 2)
legend(x = 0.2, y = max(data_a$AA)*0.3, lwd = c(NA,2, 2, 2), pch = c(1, NA, NA, NA),
       col = c("black", "red", "green", "blue"),
       legend = c("Data Integrated", "Michealis Menten", "Exponential", "Non Rectangular Hyperbola"),
       bty = "n", cex = 0.75)
}

NonRectHyperbola <- function(p, E){
(p[1] +p[2]*E - sqrt((p[1]+p[2]*E)^2-4*p[3]*p[2]*E*p[1]))/(2*p[3])
}

plot_FvCB_param_narea <- function(vpd = 0,
                                    B_lf1= 31.62 *1000^0.801,
                                # HTTP://DOI.WILEY.COM/10.1111/GCB.12870 CONVERSION OF NARE IN G M-2
                                    B_lf2=0.7,
                                    B_lf3=0.3,
                                    B_lf4=21000,
                                    B_lf5=0.801, # http://doi.wiley.com/10.1111/gcb.12870
                                    B_lf6=1.67,
                                    k_I=0.5,
                                    latitude=0.0){
require(plantecophys)
require(plant)
    assimilation_FvCB <- function(I, V, vpd, alpha, theta, lf6) {
      df_pred <- Photosyn(PPFD=I *1e+06/(24*3600), # conversion of light in mu mol /m2 /s is it ok ?
                          VPD = vpd,
                          Vcmax  = V,
                          Jmax  = V*lf6,
                          alpha = alpha,
                          theta = theta,
                          gsmodel = "BBLeuning") # Jmax ~1.67 Vcmax in Medlyn et al. 20
      return((df_pred$ALEAF + df_pred$Rd)*24 * 3600 / 1e+06)
    }







    ## Photosynthesis  [mol CO2 / m2 / yr]
    approximate_annual_assimilation <- function(narea, latitude, vpd) {
        E <- seq(0, 1, by=0.02)
        ## Only integrate over half year, as solar path is symmetrical
        D <- seq(0, 365/2, length.out = 10000)
        I <- plant:::PAR_given_solar_angle(plant:::solar_angle(D, latitude = abs(latitude)))
        Vcmax <- B_lf1 * (narea) ^  B_lf5
        theta <- B_lf2
        alpha <- B_lf3
        lf6   <- B_lf6
        AA <- NA * E
        for (i in seq_len(length(E))) {
          AA[i] <- 2 * plant:::trapezium(D, assimilation_FvCB(
                                      k_I * I * E[i], Vcmax, vpd, alpha, theta, lf6))
        }
          return(data.frame(E = E, AA = AA))
    }

N_seq <-  50
v_narea <- seq(1e-4,1e-1, length.out = N_seq)
v_vpd <- seq(0, 6, length.out = N_seq)
data_coef <- array(NA, dim = c(N_seq, N_seq, 3))
data_coef1 <- data_coef2 <- data_coef3 <- data_coef
    browser()

for (i in 1:N_seq){
  for (j in 1:N_seq){
    print(i)
    print(v_narea[i])
    print(j)
    print(v_vpd[j])
    data_a <- approximate_annual_assimilation(v_narea[i], latitude, vpd = v_vpd[j])
      tryCatch({
      library(nls2)
          print("nls2")
        fitc <- nls2(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3),
                     data = data_a,
                    algorithm = "brute-force",
                   start = expand.grid(p1 = seq(15, 6000, length.out = 30),
                                      p2 = 700,
                                      p3 = c(0.675)))
        print(coef(fitc))
        print("nls")
        fit <- nls(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3),
                   data = data_a,
                   start = coef(fitc))
        data_coef[i, j, ] <- coef(fit)
        }, error=function(e){})

      tryCatch({
        print("nls lin")
        cc <- c(387.451138727026, 55451.4186883312, -219943.812400099,
                1686197.21765355, 44.4400837954403)
        fit <- nls(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3),
                   data = data_a,
                   start = list(p1 = max(1, cc[1]+cc[2]*v_narea[i]+ cc[3]*v_narea[i]^2+
                                      cc[4]*v_narea[i]^3+cc[5]*log(v_narea[i])),
                                p2 = 700,
                                p3 = 0.7))
        data_coef1[i, j, ] <- coef(fit)
        }, error=function(e){})

      tryCatch({
      library(nlrmt)
        cc <- c(387.451138727026, 55451.4186883312, -219943.812400099,
                1686197.21765355, 44.4400837954403)
      fitxb <-     nlxb(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3),
                 data = data_a,
                 start = list(p1 =  max(1, cc[1]+cc[2]*v_narea[i]+ cc[3]*v_narea[i]^2+
                                    cc[4]*v_narea[i]^3+cc[5]*log(v_narea[i])),
                              p2 = 700,
                              p3 = 0.7),
                 lower = 0.1)

      data_coef2[i, j, ] <- fitxb$coefficients
      }, error=function(e){})

      tryCatch({
      library(nlrmt)
      fit <- wrapnls(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3),
                     data = data_a,
                      start = list(p1 =  max(1, cc[1]+cc[2]*v_narea[i]+ cc[3]*v_narea[i]^2+
                                               cc[4]*v_narea[i]^3+cc[5]*log(v_narea[i])),
                                   p2 = 700,
                                   p3 = 0.7),
                     lower = 0.1)

      data_coef3[i, j, ] <- coef(fit)
      }, error=function(e){})
## par(mfrow = c(1, 2))
##       E <- seq(0, 1, length.out = 100)
##       plot(data_a$E, data_a$AA)
##       lines(E,NonRectHyperbola(data_coef2[i,j, 1:3], E))
##       lines(E,NonRectHyperbola(data_coef[i,j, 1:3], E), col = "red")
##       lines(E,NonRectHyperbola(data_coef1[i,j, 1:3], E), col = "green")
##       lines(E,NonRectHyperbola(data_coef3[i,j, 1:3], E), col = "blue")

##     Vcmax <- B_lf1 * (v_narea[i]) ^  B_lf5
##     theta <- B_lf2
##     alpha <- B_lf3
##     lf6   <- B_lf6

##     plot(0:1500, assimilation_FvCB(I= 0:1500, V = Vcmax, vpd = v_vpd[j] , alpha, theta, lf6), type = "l")

 }
}

## colnames(data_coef) <-  c("p1", "p2", "p3", "Vcmax", "narea")
## data_coef <- data.frame(data_coef)

## par(mfrow = c(2,2))
## data_coef$narea <- v_narea
## NN <-  data_coef$narea
## plot(p1~narea, data_coef, type = "p")
## lmpoly3 <- lm(p1 ~  narea +I(narea^2)+I(narea^3)+log(narea), data_coef)
## cc <- coef(lmpoly3)
## lines(NN, cc[1]+cc[2]*NN+ cc[3]*NN^2+
##           cc[4]*NN^3+cc[5]*log(NN),
##       col = "red")
## plot(p2~narea, data_coef, type = "p")
## lmpoly4 <- lm(p2 ~  narea +I(narea^2)+I(narea^3)+I(narea^4), data_coef[data_coef$narea > 0.005, ])
## cc <- coef(lmpoly4)
## lines(NN, cc[1]+cc[2]*NN+ cc[3]*NN^2+
##           cc[4]*NN^3+cc[5]*NN^4,
##       col = "red")
## plot(p3~narea, data_coef, type = "p")
## lmpoly4 <- lm(p3 ~  narea +I(narea^2)+I(narea^3)+I(narea^4), data_coef[data_coef$narea > 0.005, ])
## cc <- coef(lmpoly4)
## lines(NN, cc[1]+cc[2]*NN+ cc[3]*NN^2+
##           cc[4]*NN^3+cc[5]*NN^4,
##       col = "red")

}


## Compare FF16 and FF16FvCB



compare_strategy_growth<- function(light = 1.0, ymax = 19){
require(plant)
derivs <- function(t, y, plant, env) {
 plant$ode_state <- y
 plant$compute_vars_phys(env)
 list(plant$ode_rates)
}

tt <- seq(0, 50, length.out=101)
env <- fixed_environment(light)

p <- scm_base_parameters("FF16")
s1 <- strategy(trait_matrix(2e-3,"narea"), p)
s2 <- strategy(trait_matrix(3e-3,"narea"), p)
s3 <- strategy(trait_matrix(9e-3,"narea"), p)
pl1 <- FF16_PlantPlus(s1)
pl2 <- FF16_PlantPlus(s2)
pl3 <- FF16_PlantPlus(s3)

y01 <- setNames(pl1$ode_state, pl1$ode_names)
yy1 <- deSolve::lsoda(y01, tt, derivs, pl1, env=env)
y02 <- setNames(pl2$ode_state, pl2$ode_names)
yy2 <- deSolve::lsoda(y02, tt, derivs, pl2, env=env)
y03 <- setNames(pl3$ode_state, pl3$ode_names)
yy3 <- deSolve::lsoda(y03, tt, derivs, pl3, env=env)
plot(height ~ time, yy1, type="l", ylim = c(0, ymax), main = paste("light level", light))
lines(yy2[, "time"], yy2[, "height"], lty = 2)
lines(yy3[, "time"], yy3[, "height"], lty = 3)


pF <- scm_base_parameters("FF16FvCB")

sF1 <- strategy(trait_matrix(2e-3,"narea"), pF)
sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
plF1 <- FF16FvCB_PlantPlus(sF1)
plF2 <- FF16FvCB_PlantPlus(sF2)
plF3 <- FF16FvCB_PlantPlus(sF3)

yF01 <- setNames(plF1$ode_state, plF1$ode_names)
yyF1 <- deSolve::lsoda(yF01, tt, derivs, plF1, env=env)
yF02 <- setNames(plF2$ode_state, plF2$ode_names)
yyF2 <- deSolve::lsoda(yF02, tt, derivs, plF2, env=env)
yF03 <- setNames(plF3$ode_state, plF3$ode_names)
yyF3 <- deSolve::lsoda(yF03, tt, derivs, plF3, env=env)
lines(height ~ time, yyF1, type="l", lty = 1, col = "red")
lines(height ~ time, yyF2, type="l", lty = 2, col = "red")
lines(height ~ time, yyF3, type="l", lty = 3, col = "red")
legend("bottomright",
       c("Narea = 2e-3", "Narea = 3e-3", "Narea = 9e-3", "FF16", "FF16FvCB"),
       lty=c(1,2,3,1,1), col=c(rep("black", 3), "red", "black"), bty="n")

}


f_g<- function(x, pl) {
  env <- fixed_environment(x)
  pl$compute_vars_phys(env)
  pl$internals$height_dt
}


compare_strategy_light_growth<- function(){
require(plant)
openness <- seq(0, 1, length.out=101)


p <- scm_base_parameters("FF16")
s1 <- strategy(trait_matrix(2e-3,"narea"), p)
s2 <- strategy(trait_matrix(3e-3,"narea"), p)
s3 <- strategy(trait_matrix(9e-3,"narea"), p)
pl1 <- FF16_PlantPlus(s1)
pl2 <- FF16_PlantPlus(s2)
pl3 <- FF16_PlantPlus(s3)

plot(openness, sapply(openness, f_g, pl1), type="l", xlim=c(0, 1),
     las=1, xlab="Canopy openness", ylab="Height growth rate (m / yr)",
     ylim = c(0, 1.4))
lines(openness, sapply(openness, f_g, pl2), lty = 2)
lines(openness, sapply(openness, f_g, pl3), lty = 3)


pF <- scm_base_parameters("FF16FvCB")

sF1 <- strategy(trait_matrix(2e-3,"narea"), pF)
sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
plF1 <- FF16FvCB_PlantPlus(sF1)
plF2 <- FF16FvCB_PlantPlus(sF2)
plF3 <- FF16FvCB_PlantPlus(sF3)

lines(openness, sapply(openness, f_g, plF1), lty = 1, col = "red")
lines(openness, sapply(openness, f_g, plF2), lty = 2, col = "red")
lines(openness, sapply(openness, f_g, plF3), lty = 3, col = "red")
legend("topleft",
       c("Narea = 2e-3", "Narea = 3e-3", "Narea = 9e-3", "FF16", "FF16FvCB"),
       lty=c(1,2,3,1,1), col=c(rep("black", 3), "red", "black"), bty="n")

}


compare_strategy_growth_narea<- function(n_length = 50){
require(plant)

f1 <- function(narea, lights){
    browser()
  p <- scm_base_parameters("FF16")
  s <- strategy(trait_matrix(narea,"narea"), p)
  pl <- FF16_PlantPlus(s)
  pl$height <- 0.5
  sapply(lights, f_g, pl)
}

f2<- function(narea, lights){
    print(narea)
    browser()
  p <- scm_base_parameters("FF16FvCB")
  s <- strategy(trait_matrix(narea,"narea"), p)
  pl <- FF16FvCB_PlantPlus(s)
  pl$height <- 0.5
  sapply(lights, f_g, pl)
}
f2b <- function(n, lights){
tryCatch(f2(n, lights), error=function(err) NA)
}

v_narea <- 10^seq(log10(1e-4),log10(1e-1), length.out = n_length)
v_g <- t(sapply(v_narea, f1, lights = c(0.25, 1.0)))
v_g_F <- t(sapply(v_narea, f2, lights = c(0.25, 1.0)))

plot(v_narea, v_g[, 2], type="l", ylim = range(v_g, v_g_F),
     las=1, xlab="Nitrogen per area", ylab="Height growth rate", log = "x")
lines(v_narea, v_g[, 1], lty = 2)
lines(v_narea, v_g_F[, 2], col = "red")
lines(v_narea, v_g_F[, 1], col = "red", lty = 2)
}




# much less accurate than what is normaly done in Plant but just while lcp_whole_plant doesn't work ...
lcp_georges <- function(pl){
  openness <- seq(0, 1, length.out=501)
  v <- sapply(openness, f_g, pl)
  lcp <- mean(openness[c(max(which(v < 1e-10)), min(which(v > 1e-10)))])
  return(lcp)
}


compare_strategy_lcp<- function(){
require(plant)

p <- scm_base_parameters("FF16")
s1 <- strategy(trait_matrix(2e-3,"narea"), p)
s2 <- strategy(trait_matrix(3e-3,"narea"), p)
s3 <- strategy(trait_matrix(9e-3,"narea"), p)
pl1 <- FF16_PlantPlus(s1)
pl2 <- FF16_PlantPlus(s2)
pl3 <- FF16_PlantPlus(s3)

openness <- seq(0, 1, length.out=51)

lcp1 <- lcp_georges(pl1)
lcp2 <- lcp_georges(pl2)
lcp3 <- lcp_georges(pl3)

x <- c(lcp1, openness[openness > lcp1])
plot(x, sapply(x, f_g, pl1), type="l", xlim=c(0, 1), ylim = c(0, 1.5),
     las=1, xlab="Canopy openness", ylab="Height growth rate (m / yr)")
points(lcp1, 0.0, pch=19)
x <- c(lcp2, openness[openness > lcp2])
lines(x, sapply(x, f_g, pl2), lty = 2)
points(lcp2, 0.0, pch=19)
x <- c(lcp3, openness[openness > lcp3])
lines(x, sapply(x, f_g, pl3), lty = 3)
points(lcp3, 0.0, pch=19)

## FF16FvCB
pF <- scm_base_parameters("FF16FvCB")

sF1 <- strategy(trait_matrix(1e-3,"narea"), pF)
sF2 <- strategy(trait_matrix(3e-3,"narea"), pF)
sF3 <- strategy(trait_matrix(9e-3,"narea"), pF)
plF1 <- FF16FvCB_PlantPlus(sF1)
plF2 <- FF16FvCB_PlantPlus(sF2)
plF3 <- FF16FvCB_PlantPlus(sF3)

lcpF1 <- lcp_georges(plF1)
lcpF2 <- lcp_georges(plF2)
lcpF3 <- lcp_georges(plF3)

x <- c(lcpF1, openness[openness > lcpF1])
lines(x, sapply(x, f_g, plF1), col = 'red')
points(lcpF1, 0.0, pch=19, col = "red")
x <- c(lcpF2, openness[openness > lcpF2])
lines(x, sapply(x, f_g, plF2), lty = 2, col = 'red')
points(lcpF2, 0.0, pch=19, col = "red")
x <- c(lcpF3, openness[openness > lcpF3])
lines(x, sapply(x, f_g, plF3), lty = 3, col = 'red')
points(lcpF3, 0.0, pch=19, col = "red")
}


## lcp in function of Narea
compare_strategy_lcp_narea<- function(n_length = 50){
require(plant)

f1 <- function(n){
  p <- scm_base_parameters("FF16")
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16_PlantPlus(s)
  pl$height <- 0.5
  lcp_georges(pl)
}

f2<- function(n){
  p <- scm_base_parameters("FF16FvCB")
  s <- strategy(trait_matrix(n,"narea"), p)
  pl <- FF16FvCB_PlantPlus(s)
  pl$height <- 0.5
  lcp_georges(pl)
}
f2b <- function(n){
tryCatch(f2(n), error=function(err) NA)
}

v_narea <- 10^seq(log10(1e-4),log10(1e-1), length.out = n_length)
v_lcp <- sapply(v_narea, f1)
v_lcp_F <- sapply(v_narea, f2b)

plot(v_narea, v_lcp, type="l", ylim = c(0, 1),
     las=1, xlab="Nitrogen per area", ylab="WP LCP",
     log = "x")
lines(v_narea, v_lcp_F, col = "red")
}




## look at patch dynamcis for two different N value

compare_strategy_patch<- function(){
require(plant)
p0 <- scm_base_parameters("FF16")
p0$control$equilibrium_nsteps <- 30
p0$control$equilibrium_solver_name <- "hybrid"
p0$disturbance_mean_interval <- 30.0

## First, with a single species:
p <- expand_parameters(trait_matrix(c(3e-3, 5e-3, 9e-3), "narea"), p0, FALSE)


p_eq <- equilibrium_seed_rain(p)

## Then collect the patch-level dynamics:
data <- run_scm_collect(p_eq)

tt <- data$time
h1 <- data$species[[1]]["height", , ]
h2 <- data$species[[2]]["height", , ]
h3 <- data$species[[3]]["height", , ]

cols <- c("#e41a1c", "#377eb8", "#4daf4a")

## Relativise the log densities onto (-4, max)
d1 <- data$species[[1]]["log_density", , ]
d2 <- data$species[[2]]["log_density", , ]
d3 <- data$species[[3]]["log_density", , ]
rel <- function(x, xmin) {
  x[x < xmin] <- xmin
  xmax <- max(x, na.rm=TRUE)
  (x - xmin) / (xmax - xmin)
}



rd1 <- rel(d1, -4)
rd2 <- rel(d2, -4)
rd3 <- rel(d3, -4)

n <- length(tt)
x1 <- matrix(rep(tt, ncol(h1)), nrow(h1))
x2 <- matrix(rep(tt, ncol(h2)), nrow(h2))
x3 <- matrix(rep(tt, ncol(h3)), nrow(h3))
col1 <- matrix(make_transparent(cols[[1]], rd1), nrow(d1))
col2 <- matrix(make_transparent(cols[[2]], rd2), nrow(d2))
col3 <- matrix(make_transparent(cols[[3]], rd3), nrow(d3))
par(mfrow = c(1,2))
plot(NA, xlim=range(tt),
     las=1, xlab="Time (years)", ylab="Cohort height (m)",
     ylim = range(h1,h2,h3, na.rm = TRUE))
segments(x1[-1, ], h1[-1, ], x1[-n, ], h1[-n, ], col=col1[-n, ], lend="butt")
segments(x2[-1, ], h2[-1, ], x2[-n, ], h2[-n, ], col=col2[-n, ], lend="butt")
segments(x3[-1, ], h3[-1, ], x3[-n, ], h3[-n, ], col=col3[-n, ], lend="butt")

## FvCB
pF0 <- scm_base_parameters("FF16FvCB")
pF0$control$equilibrium_nsteps <- 30
pF0$control$equilibrium_solver_name <- "hybrid"
pF0$disturbance_mean_interval <- 30.0

## First, with a single species:
pF <- expand_parameters(trait_matrix(c(3e-3, 5e-3, 9e-3), "narea"), pF0, FALSE)

pF_eq <- equilibrium_seed_rain(pF)

## Then collect the patch-level dynamics:
dataF <- run_scm_collect(pF_eq)

tt <- dataF$time
h1 <- dataF$species[[1]]["height", , ]
h2 <- dataF$species[[2]]["height", , ]
h3 <- dataF$species[[3]]["height", , ]

cols <- c("#e41a1c", "#377eb8", "#4daf4a")

## Relativise the log densities onto (-4, max)
d1 <- dataF$species[[1]]["log_density", , ]
d2 <- dataF$species[[2]]["log_density", , ]
d3 <- dataF$species[[3]]["log_density", , ]

rd1 <- rel(d1, -4)
rd2 <- rel(d2, -4)
rd3 <- rel(d3, -4)
n <- length(tt)
x1 <- matrix(rep(tt, ncol(h1)), nrow(h1))
x2 <- matrix(rep(tt, ncol(h2)), nrow(h2))
x3 <- matrix(rep(tt, ncol(h3)), nrow(h3))
col1 <- matrix(make_transparent(cols[[1]], rd1), nrow(d1))
col2 <- matrix(make_transparent(cols[[2]], rd2), nrow(d2))
col3 <- matrix(make_transparent(cols[[3]], rd3), nrow(d3))
plot(NA, xlim=range(tt),
     las=1, xlab="Time (years)", ylab="Cohort height (m)",
     ylim = range(h1,h2,h3, na.rm = TRUE))
segments(x1[-1,], h1[-1, ], x1[-n, ], h1[-n, ], col=col1[-n, ], lend="butt")
segments(x2[-1, ], h2[-1, ], x2[-n, ], h2[-n, ], col=col2[-n, ], lend="butt")
segments(x3[-1, ], h3[-1, ], x3[-n, ], h3[-n, ], col=col3[-n, ], lend="butt")
legend("right", legend = c("Narea = 3e-3", "Narea = 5e-3", "Narea = 9e-3"), lty=1, col = cols, bty="n")
}


compare_fit_Narea_Vcmax <- function(){
a1 <- c(1.99, 6.35, 5.4, 34.05, 4.61)
b1 <- c(10.71, 25.88, 30.38, 9.71, 30.20)
# when estimation use also Amax not only Vcmax
a2 <- c(8.9, 4.19, 5.73, 6.32, 14.71)
b2 <-  c(9.30, 26.19, 26.81, 18.15, 23.15)
names(a1) <- names(b1) <- names(a2) <- names(b2) <-c("Tropical_Trees_Oxisols",
                                                    "Tropical_Trees_NonOxisols",
                                                     "Temperate_Trees",
                                                     "Coniferous_Trees",
                                                     "Shrubs")
Narea <- seq(0.1, 10, length.out = 100)
par(mfrow = c(1,2))
plot(Narea, 31.62*Narea^0.801, type = "l", xlab = "Narea (g m-2)",
     ylab = "Vcmax")
lines(Narea, exp(3.712)*Narea^0.650, lty = 2)
lines(Narea, 10^1.57*Narea^0.55, lty = 3)
lines(Narea, mean(a1) + Narea*mean(b1), lty = 4 , col = 2)
lines(Narea, mean(a2) + Narea*mean(b2), lty = 4 , col = 3)
legend("topleft", legend = c("Sakschewski-TRY", "Walker2017",
                             "Domingues2011", "Kattge2009_V_meanPFT",
                             "Kattge2009_VA_meanPFT"),
       lty = c(1:4,4), col = c(1,1,1,2,3))
plot(Narea, 31.62*Narea^0.801, type = "l", lty = 0, xlab = "Narea (g m-2)",
     ylab = "Vcmax")
for (i in 1:5){
    lines(Narea, a1[i] + Narea*b1[i], lty = 1 , col = i+1)
    lines(Narea, a2[i] + Narea*b2[i], lty = 2 , col = i+1)
}
legend("topleft", legend = names(a1),
       lty = 1, col = 2:6)
}


compare_fit_Vcmax_Jmax <-  function(){
V <- 5:200
plot(V, V*1.67, type = "l", xlab = "Vcmax", ylab = "Jmax")
lines(V, exp(1.01)*V^0.89, lty = 2)
lines(V, exp(1.669)*V^0.75, lty = 3)
lines(V, exp(1.425)*V^0.837, lty = 4)
legend("topleft", legend = c("Medlyn", "Walker", "Kattge", "Wullschleger"),
       lty = 1:4)
}



