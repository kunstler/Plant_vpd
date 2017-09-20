#Look at use of FcV model with plantecophys in Plant

# Default param for FcV for forest tree ??
# In plantecophys  alpha = 0.24, theta = 0.85, Jmax = 100, Vcmax = 50
# in Medlyn et al. 2002 alpha = 0.3 theta 0.9
# In Plant theta = 0.5 alpha is not on the same unite (mol of CO2 per mol of photon wherease in FcV in mol electron per mol photon) = 0.04 (conversion see Mercado et al 2009 , the conversion of the apparent quantum yield in micromolCO2/micromol quantum into micromol e-/micxromol quantum is done by multipliyng by 4, since four electrons are needed to regenerate RuBP 4*0.04 = 0.16 but we also probably need to correct for leaf absorbance 0.8 in Medlyn et al. 2002)
# Sterck et al. 2011 theta = 0.5 alpha = 0.25 (not exactly the same alpha ...)
# Valladares et al. 1997 tropicla shrub theta between 0.51 and 0.8 and alpha (in mol CO2) between 0.048 and 0.066
# In TROLL Marechau & Chave Ecol. Monogr. in press alpha = 0.075*4 = 0.3 and theta = 0.7 based on von Caemmerer 2000

# In von  Caemmerer 2000 the non rectangular hyperbola light response curve have a initial slope $\alpha_1 \times \PHI_{PSII} \times \beta$ with $\alpha_1$ the leaf absorptance, $\PHI_{PSII}$ the maximum quantum yield of the photosystem II and $\beta$ the fraction of absorbed light that reach PSII.

# LINK Vcmax and Jmax (according to Medlyn et al. 2002 Jmax = 1.67 Vcmax)?
# in Walker, A. P. et al. 2014. The relationship of leaf photosynthetic traits - Vcmax and Jmax - to leaf nitrogen, leaf phosphorus, and specific leaf area: a meta-analysis and modeling study. - Ecology and Evolution 4: 3218–3235. It is a power function exp(1.01)*V^0.89 in Kattge et al. 2011 exp(1.669)*V^0.75, in Wullschleger exp(1.425)*V^0.837

#Vcmax vs Narea
# Sakschewski et al. 2014 Vcmax = 31.62*Narea^0.801
# Or Vcmax and Jmax functiuon of N, P and LMA see Domingues, T. F. et al. 2010. Co-limitation of photosynthetic capacity by nitrogen and phosphorus in West Africa woodlands: Nutrients constraints on photosynthesis in African woodlands. - Plant, Cell & Environment 33: 959–980.
#

## V <- 5:200
## plot(V, V*1.67, type = "l", xlab = "Vcmax", ylab = "Jmax")
## lines(V, exp(1.01)*V^0.89, lty = 2)
## lines(V, exp(1.669)*V^0.75, lty = 3)
## lines(V, exp(1.425)*V^0.837, lty = 3)



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
list_df_PPFD_V_VPD1.9<- vector("list")
for (i in seq_len(length(V_seq))){
  list_df_PPFD_V_VPD0[[i]]<- Photosyn(VPD = 0, PPFD = 1:1500,
                            Vcmax  = V_seq[i],
                            Jmax  = V_seq[i]*1.67,
                            gsmodel = gs_type)$ALEAF
  list_df_PPFD_V_VPD1.9[[i]]<- Photosyn(VPD = 1.9, PPFD = 1:1500,
                            Vcmax  = V_seq[i],
                            Jmax  = V_seq[i]*1.67,
                            gsmodel = gs_type)$ALEAF
}

df_PPFD_V_VPD0<- as.data.frame(list_df_PPFD_V_VPD0)
df_PPFD_V_VPD1.9<- as.data.frame(list_df_PPFD_V_VPD1.9)
plot(1:1500, df_PPFD_V_VPD0[,1],
     ylim = range(list_df_PPFD_V_VPD0, list_df_PPFD_V_VPD1.9),
     xlab = expression(paste("PPFD (", mu, "mol", " ",m^-2, " ", s^-1, ")" )),
     ylab = expression(paste(A[n],"(", mu, "mol", " ", m^-2, " ",s^-1, ")" )),
     type = "l", col = cols[1], lty = 1)
for (i in 2:length(V_seq)) lines(1:1500, df_PPFD_V_VPD0[, i], col = cols[i], lty = 1, lwd = 2)
for (i in 1:length(V_seq)) lines(1:1500, df_PPFD_V_VPD1.9[, i], col = cols[i], lty = 2, lwd = 2)
legend(x = 0, y = max(df_PPFD_V_VPD0)*0.9, lwd = 2, lty = c(1, 2, rep(1,4)),
       col = c("black", "black", cols),
       legend = c("VPD = 0", "VPD = 1.9", paste("Vcmax", round(V_seq))),
       bty = "n", cex = 0.5)
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
Compare_Photo_Plant_FcV <- function(){
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
       legend = c("Plant", "FcV_plantecophys", "FcV_plantecophys V = 50 J = 1.67 V",
                  "FcV Plant and V =32.5 J = 1.67 J", "FcV Troll = von Caemmerer 2000", "Sterck et al. 2011" ),
       bty = "n")
}

Photo_Plant_match_FcV <- function(){
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
       legend = c("Plant", "FcV_plantecophys"),
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
lines(E, a_p1*E/(a_p2+E))
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
           start = list(p1 = 120, p2 = 400, p3 = 0.9))
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



## plot_photosyn_annual_FvC(vpd = 0, alpha = 0.3, theta = 0.7)

# Sperry model of stomatal conductance TODO ADD (WHAT full model or just an approximation based on only on cavitation curve ? see Sterck for similar idea)

# cavitation curve
k_xyl <- function(P, b, c, kmax = 5){
kmax * exp(-(P/b)^c)
}

