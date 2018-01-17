
#looks up MAT and MAP for given lat/lon values
GetClimate <-function(lats,lons, wc) {

  plot.mat <- numeric(length(lats))
  plot.tmax <- numeric(length(lats))
  plot.map <- numeric(length(lats))
  plot.pdry <- numeric(length(lats))

  plot.mat <- extract(wc,
                      cbind(lons,
                            lats),
                      layer=1,nl=1)/10
  #retrieve MAT, divide by 10 for deg C
  plot.tmax <- extract(wc,
                       cbind(lons,
                             lats),
                       layer=5,nl=1)/10
  #retrieve T max of warmest month, divide by 10 for deg C
  plot.map <- extract(wc,
                      cbind(lons,
                            lats),
                      layer=12,nl=1)
  #retrieve MAP
  plot.pdry <- extract(wc,
                       cbind(lons,
                             lats),
                       layer=14,nl=1)
  #retrieve P dryest month
  return(data.frame(MAT=plot.mat,MAP=plot.map, TMAX = plot.tmax, PDRY = plot.pdry))

} #GetClimate


#get aridity index

download_aridity_index<- function(dir_temp = "download"){
require(raster)
require(R.utils)
# download raster

url_clim <- "https://www.dropbox.com/sh/e5is592zafvovwf/AACSS163OQ2nm5m1jmlZk4Gva/Global%20PET%20and%20Aridity%20Index/Global%20Aridity%20-%20Annual.zip?dl=1"
raster_name_zip <- "Global Aridity - Annual.zip"
if(!dir.exists(file.path(dir_temp, "AI_annual"))){
    dir.create(dir_temp)
    download.file(url_clim, file.path(dir_temp, raster_name_zip))
    unzip(zipfile = file.path(dir_temp, raster_name_zip),
          exdir = dir_temp)
}

raster_aridity <- raster(file.path(dir_temp, "AI_annual",
                                   "ai_yr",  "w001001x.adf"))

return(raster_aridity)
}

extract_aridty <- function(lats,lons,raster_aridity){
plot_aridity <- raster::extract(raster_aridity, cbind(lons, lats))/10000
nas <- which(is.na(plot_aridity))
    if (length(nas)>0) {
      good <- which(!is.na(plot_aridity))
      near <- NearestPlot(lons[nas],lats[nas],lons[good],lats[good])
      plot_aridity[nas] <- plot_aridity[good[near]]
    }
return(plot_aridity)
}

#Returns the index of the nearest plot(x2,y2) to each of the target plots (x1,y1)
NearestPlot <- function(x1,y1,x2,y2) {

  res = numeric(length(x1))

  for (i in 1:length(x1)) {
    dists = sqrt((x1[i] - x2)^2 + (y1[i] - y2)^2)
    res[i] = which.min(dists)
    }
  return(res)

} #NearestPlot


############
#### COMPUTE VPD AND ITS CORRELATION GLOBALY

# STATURATED VAPOR PRESSURE FUNCTIONS


#Teten formula in FAO 56 in
svp <- function(T){
0.6108 * exp(17.27 * T / (T + 237.3))
}

fun_unzip <- function(name, down_path = "download"){
unzip(name, junkpaths = TRUE, exdir = down_path)
}

read_all_month_stack <- function(name){
months <- c(paste0("0", 1:9), 10:12)
f <- function(n, name) raster(file.path("download", paste0(name,n,".tif")))
stack(lapply(months, f, name = name))
}

vpd_tavg_corr <- function(){
stack_tmin <- read_all_month_stack("wc2.0_2.5m_tmin_")
stack_tmax <- read_all_month_stack("wc2.0_2.5m_tmax_")
stack_tavg <- read_all_month_stack("wc2.0_2.5m_tavg_")
stack_vapr <- read_all_month_stack("wc2.0_2.5m_vapr_")
stack_t <- stack(stack_tmin, stack_tmax)
f <- function(x) (svp(x[1:12]) +svp(x[13:24]))/2
stack_svp <- calc(stack_t, f)
stack_vp <- stack(stack_vapr, stack_svp)
f <- function(x) (x[1:12] -x[13:24])
stack_vpd <- calc(stack_vp,f)
return(list(vpd = stack_vpd, tavg = stack_tavg))
}

mean_vpd_t_cor <- function(list_clim){
vpd_m <- calc(list_clim$vpd, function(x) mean(x[1:12]))
t_m <- calc(list_clim$tavg, function(x) mean(x[1:12]))
df2 <- data.frame(T = as.vector(t_m), vpd = as.vector(vpd_m))
list_vpd_t<- vector("list")
for (i in 1:12){
list_vpd_t[[i]] <- data.frame(T = as.vector(list_clim$tavg[[i]]), vpd = as.vector(list_clim$vpd[[i]]),
                              months = rep(i, length.out = dim(list_clim$vpd[[i]])[1]*
                                                          dim(list_clim$vpd[[i]])[2]))
}

df1 <- do.call("rbind",list_vpd_t)
# restrict data to be above -15C
df2 <- df2[df2$T> -15 & !is.na(df2$T), ]

df1 <- df1[df1$T> -40 & !is.na(df1$T), ]
return(list(df_months= df1, df_avg= df2))
}

plot_reg_vpd_t <- function(ll_d){
df1 <- ll_d$months
df2 <- ll_d$avg

df1$vpd2 <- df1$vpd - max(df1$vpd)*1.01
df2$vpd2 <- df2$vpd - max(df2$vpd)*1.01
browser()

res2 <- lm(T~log(-vpd2), data = df2)
res1 <- lm(T~log(-vpd2), data = df2)

coef1 <- coefficients(res1)
coef2 <- coefficients(res2)
seq_vpd <- seq(-4, 0, length.out = 100)
pred1 <- coef1[1] +coef1[2]*log(seq_vpd - max(df1$vpd)*1.01)
pred2 <- coef2[1] +coef2[2]*log(seq_vpd - max(df2$vpd)*1.01)

pdf("figures/vpd_t.pdf")
x11()
plot(df1$vpd, df1$T, col = df1$months, xlab = "vpd", ylab ="Tavg months", main = "monthly data")
lines(seq_vpd, pred1, lwd =2)
print(coef1)
x11()
plot(df2$vpd, df2$T, xlab = "vpd", ylab ="Tavg months", main = "annual data")
lines(seq_vpd, pred2, lwd =2)
print(coef2)
dev.off()
}

