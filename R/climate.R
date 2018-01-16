
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
browser()

}
## T <- seq(-10, 30, length.out = 20)
## plot(T, svp_a(T)/10, type = "l", ylim = range(svp_a(T)/10, svp_b(T)))
## lines(T, svp_b(T), lty = 2)
## lines(T, svp_c(T)/10, lty = 3, col = "red")
## lines(T, svp_d(T)/10, lty = 4)


