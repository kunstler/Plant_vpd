
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


