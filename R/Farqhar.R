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

N_seq <-  75
v_narea <- seq(1e-4,1e-1, length.out = N_seq)
v_vpd <- seq(0, 6, length.out = N_seq)
data_coef <- array(NA, dim = c(N_seq, N_seq, 3))

for (i in 1:N_seq){
  for (j in 1:N_seq){
    print(i)
    print(v_narea[i])
    print(j)
    print(v_vpd[j])
    data_a <- approximate_annual_assimilation(v_narea[i], latitude, vpd = v_vpd[j])

      tryCatch({
      library(nlmrt)
      fitxb2 <- nlxb(AA ~ (p1 +p2*E - sqrt((p1+p2*E)^2-4*p3*p2*E*p1))/(2*p3),
                           data = data_a,
                           start = list(p1 = 500,
                                        p2 = 600,
                                        p3 = 0.8),
                           lower = 0.01,
                           upper = c(7000, 1000, 1))
      data_coef[i, j, ] <- fitxb2$coefficients
      rm(fitw)
      }, error=function(e){})

 }
}

df <- data_coef
df_p<- data.frame(narea = rep(v_narea, N_seq),
                 vpd = rep(v_vpd, each = N_seq),
                 p1 = as.vector(df[,,1]),
                 p2 = as.vector(df[,,2]),
                 p3 = as.vector(df[,,3]))
library(ggplot2)
p1 <- ggplot(df_p, aes(narea, vpd)) +
  geom_raster(aes(fill = p1)) +
  scale_fill_gradientn(colours = terrain.colors(10))
p2 <- ggplot(df_p, aes(narea, vpd)) +
  geom_raster(aes(fill = p2)) +
  scale_fill_gradientn(colours = terrain.colors(10))
p3 <- ggplot(df_p, aes(narea, vpd)) +
  geom_raster(aes(fill = p3)) +
  scale_fill_gradientn(colours = terrain.colors(10))

lm_p1 <- lm(p1~narea*vpd+ I(narea^2) + I(vpd^2) + I(narea^2):I(vpd^2)+
                I(narea^3) + I(vpd^3) + I(narea^3):I(vpd^3) , df_p)
ppb1 <- ggplot(data.frame(obs = df_p$p1, pred = predict(lm_p1)), aes(obs, pred)) + geom_point()+
        geom_abline(slope = 1, intercept = 0)
df_p$p1_pred <- predict(lm_p1)
pp1 <- ggplot(df_p, aes(narea, vpd)) +
  geom_raster(aes(fill = p1_pred)) +
  scale_fill_gradientn(colours = terrain.colors(10))
lm_p2 <- lm(p2~narea*vpd+ I(narea^2) + I(vpd^2) + I(narea^2):I(vpd^2)+
                I(narea^3) + I(vpd^3) + I(narea^3):I(vpd^3), df_p)
ppb2 <- ggplot(data.frame(obs = df_p$p2, pred = predict(lm_p2)), aes(obs, pred)) + geom_point()+
    geom_abline(slope = 1, intercept = 0)
df_p$p2_pred <- predict(lm_p2)
pp2 <- ggplot(df_p, aes(narea, vpd)) +
  geom_raster(aes(fill = p2_pred)) +
  scale_fill_gradientn(colours = terrain.colors(10))
lm_p3 <- lm(p3~narea*vpd+ I(narea^2) + I(vpd^2) + I(narea^2):I(vpd^2)+
                I(narea^3) + I(vpd^3) + I(narea^3):I(vpd^3), df_p)
ppb3 <- ggplot(data.frame(obs = df_p$p3, pred = predict(lm_p3)), aes(obs, pred)) + geom_point()+
    geom_abline(slope = 1, intercept = 0)
df_p$p3_pred <- predict(lm_p3)
pp3 <- ggplot(df_p, aes(narea, vpd)) +
  geom_raster(aes(fill = p3_pred)) +
  scale_fill_gradientn(colours = terrain.colors(10))

## pp <- multiplot(p1, p2, p3, pp1, pp2, pp3, ppb1, ppb2, ppb3,cols = 3)

pdf("../figures/NRH_params.pdf", width = 12, height = 10)
pp <- multiplot(p1, p2, p3, cols = 2)
dev.off()

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



