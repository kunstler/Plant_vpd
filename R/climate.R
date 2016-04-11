
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



