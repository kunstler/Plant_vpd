process_wright_2004 <- function(filename, sitevars_file, AridityIndex) {

  ## There are several strategies for reading in an excel file, but
  ## this one works quite well.
  d <-read.xlsx2(filename, sheetIndex=1, startRow=11,
                  stringsAsFactors=FALSE, check.names=FALSE)

  ## Do some name translations:
  tr <- c("Code"="Code",
          "Dataset"="Dataset",
          "BIOME"="Biome",
          "Species"="Species",
          "GF"="GrowthForm",
          "Decid/E'green"="Deciduous",
          "Needle/Broad lf"="Needle",
          "C3C4"="C3",
          "N2-fixer"="N2fixer",
          "log LL"="LogLeafLifespan",
          "log LMA"="LogLMA",
          "log Nmass"="Log.N.mass",
          "log Narea"="Log.N.area",
          "log Pmass"="Log.P.mass",
          "log Parea"="Log.P.area",
          "log Amass"="Log.A.mass",
          "log Aarea"="Log.A.area",
          "log Gs"="Log.Gs",
          "log Rdmass"="Log.Rd.mass",
          "log Rdarea"="Log.Rd.area",
          "Ca - Ci"="CaCi")
  names(d)[match(names(tr), names(d))] <- tr

  ## Drop blank columns
  d <- d[names(d) != " "]

  ## Data tweaking:
  d[["Code"]] <- as.integer(d[["Code"]])
  d[["CaCi"]] <- as.numeric(d[["CaCi"]])

  d[["Deciduous"]] <- category_to_logical(d[["Deciduous"]], "D")
  d[["Needle"]]    <- category_to_logical(d[["Needle"]],    "N")
  d[["C3"]]        <- category_to_logical(d[["C3"]],        "C3")
  d[["N2fixer"]]   <- category_to_logical(d[["N2fixer"]],   "Y")

  names(d) <- gsub("Log\\.", "Log", names(d))
  re <- "Log"
  i_log <- grep(re, names(d))
  d[i_log] <- lapply(d[i_log], as.numeric)
  d_unlogged <- as.data.frame(10^d[i_log])
  names(d_unlogged) <- sub(re, "", names(d_unlogged))

  d <- cbind(d[-c(i_log)], d_unlogged)

  # add location info
  sitevars <- read.csv(sitevars_file, stringsAsFactors = FALSE,
                       encoding = "UTF-8")
  ## add aridity index
  sitevars$AI <- extract_aridty(sitevars$latitude, sitevars$longitude,AridityIndex)

  data <- merge(d, sitevars, by.x = 'Dataset', by.y = 'dataset_location',
                all.x = TRUE, all.y = FALSE, sort = FALSE)

  #lowercase names
  names(data) <- tolower(names(data))

  # unit conversions
  data$lma <- data$lma/1000 # Converts to kg
  data$n.area <- data$n.area/1000 # Converts to kg
  data$a.area <- (data$a.area * 31557.6)*10^-6
       # converts to mol/kg/yr from micro-mol/g/s
  data$rd.area <- (data$rd.area * 31557.6)*10^-6 # converts to mol/kg/yr from micro-mol/g/s

  data$leaflifespan <- data$leaflifespan/12 ## convert LL from months to years
  data$leaf_turnover <- 1/data$leaflifespan ## per year

  data$mat_o_map<-  (data$mat+18)/data$map

  names(data)[names(data) %in%
              c('dataset', 'mat_degc', 'map_mm')]<- c('location','mat', 'map')

  ## return(subset(data,   data[["map"]] < 4000 & #why did you subset the data to have more than 400m of MAP ??
  ##    data[["growthform"]] %in% c("S","T")))
  return(data)
}

figure_lma_climate <- function(data) {
  lma <- data[["lma"]]
  narea <- data[["n.area"]]
  mat_o_map <- data[["mat_o_map"]]
  map <- data[["map"]]
  ai <- data[["ai"]]


  myplot <- function(...,  xlabel=FALSE,  ylabel=FALSE) {
    plot(..., ann=FALSE, xaxt='n', yaxt='n', pch=16, cex=0.9)
    axis(1, labels=xlabel)
    axis(2, labels=ylabel, las=1)
  }

  mylabel <- function(lab, ...) label(0.05, 0.95, lab, log.y=TRUE, ...)

  par(mfrow=c(2,3), mar=c(1,1,1,1), oma=c(4,5,0,0))

  myplot(map, lma, log ="xy",  xlabel=FALSE, ylabel=TRUE)
  mtext(expression(paste("Leaf-mass per area (kg ", m^-2,")")), 2, line=4)
  mylabel("a", log.x=TRUE)

  myplot(mat_o_map, lma, log ="xy", xlabel=FALSE)
  mylabel("b", log.x=TRUE)

  myplot(ai, lma, log ="xy", xlabel=FALSE)
  mylabel("c", log.x=TRUE)

  myplot(map, narea, log ="xy",xlabel=TRUE, ylabel=TRUE)
  mtext("Annual precipitation (mm)",1, line=4)
  mtext(expression(paste("Leaf-nitrogen per area (kg ", m^-2,")")), 2, line=4)
  mylabel("d", log.x=TRUE)

  myplot(mat_o_map, narea, log ="xy", xlabel=TRUE)
  mtext("Temperature (C) over Precipitation (mm)", 1, line=4)
  mylabel("e", log.x=TRUE)

  myplot(ai, narea, log ="xy", xlabel=TRUE)
  mtext("Aridity Index MAP/PET", 1, line=4)
  mylabel("d", log.x=TRUE)
}

figure_lma_tradeoff <- function(data) {
  data <- subset(data, !is.na(data[["lma"]] * data[["leaf_turnover"]])
    & table(data[["location"]])[data[["location"]]] > 9)
  location <- data[["location"]]
  lma <- data[["lma"]]
  leaf_turnover <- data[["leaf_turnover"]]

  sm1 <- sma(leaf_turnover ~ lma * location, log="xy")

  par(mar=c(4.6, 4.6, .5, .5))
  plot(NA, type="n", log="xy", xlim=c(0.01, 1.28), ylim=c(0.03, 32),
       xlab="", ylab="", las=1)
  mtext(expression(paste("Leaf-construction cost (kg ", m^-2,")")), line=3, side = 1)
  mtext(expression(paste("Leaf turnover rate (",yr^-1,")")), line=3, side = 2)

  points(lma, leaf_turnover, col=make_transparent("grey", 0.5), pch=16)
  plot(sm1, add=TRUE, type="l", lwd=1.25, p.lines.transparent=0.15)

  x <- seq_log_range(c(0.001,3), 40)
  points(x, 0.0286*x^-1.71, type='l', col='black', lwd=2)

  title <- sprintf("%d sites, %d species",
                   length(unique(location)),
                   sum(!is.na(leaf_turnover)))

}


figure_lma_tradeoff_climate <- function(data,
                                        var_clim = 'mat_o_map') {
  data <- subset(data, !is.na(data[["lma"]] * data[["leaf_turnover"]]))
  lma <- data[["lma"]]
  leaf_turnover <- data[["leaf_turnover"]]
  if(var_clim == 'map'){
    data$levels <- cut(data[[var_clim]],
                     breaks = c(400, 600, 900, 1200, 2400, 4800) ,
                     labels = FALSE, ordered_result = TRUE)
  }
  if(var_clim == 'ai'){
   data$levels <- cut(data[[var_clim]],
                    breaks = c(0, 0.6, 1, 1.4, 2, 4) ,
                    labels = FALSE, ordered_result = TRUE)
   }
  if(var_clim == "mat_o_map"){
  data$levels <- cut(data[[var_clim]],
                     breaks = c(-0.03, 0, 0.01, 0.02, 0.04, 0.06),
                     labels = FALSE,
                     ordered_result = TRUE)
  }


  groups<- data[['levels']]
  sm1 <- sma(leaf_turnover ~ lma * groups, log="xy")

  colfunc <- colorRampPalette(c("red", "blue"))
  cols <- colfunc(length(unique(data[['levels']])))
  if(var_clim == "mat_o_map") rev(cols)
  col_sm1 <- cols[data[["levels"]][match(sm1[["groups"]], groups)]]

  par(mar=c(4.6, 4.6, .5, .5))
  plot(NA, type="n", log="xy", xlim=c(0.01, 1.28), ylim=c(0.03, 32),
       xlab="", ylab="", las=1)
  mtext(expression(paste("Leaf-construction cost (kg ", m^-2,")")),
        line=3, side = 1)
  mtext(expression(paste("Leaf turnover rate (",yr^-1,")")), line=3, side = 2)

  points(lma, leaf_turnover, col=cols[data$levels],
         pch=16)
  plot(sm1, add=TRUE, col=col_sm1, type="l", lwd=2)

  x <- seq_log_range(c(0.001,3), 40)
  points(x, 0.0286*x^-1.71, type='l', col='black', lwd=1, lty = 2)

  title <- sprintf("%d sites, %d species",
                   length(unique(data$location)),
                   sum(!is.na(leaf_turnover)))
  legend("topright", legend=paste(var_clim, " class",
                                  1:length(unique(data[['levels']]))), bty="n",
         pch=16, col=cols, cex=1, title=title)

}

plot_coef_sma <- function(df, var){
  plot(df[[var]], df$elevationm,
       ylim = range(df$elevationch, df$elevationcl),
       pch = 16, xlab = var, ylab = "SMA LTR-LCC-tradeoff intercept")
  segments(df[[var]], df$elevationcl, df[[var]], df$elevationch)
  obj_lm <- lm(formula(paste("elevationm ~ ", var)), data = df ,
            weights = df$elevw)
  abline(obj_lm,
         col = 'gray')
  if(summary(obj_lm)$coefficients[2,4] < 0.01) {
   text(x=min(df[[var]]), y=max(df$elevationch)*0.85,
        paste0("* R2 =",round(summary(obj_lm)$r.squared, 3)) , cex=1.2, col = "red", pos = 4)
    }else{
      text(x=min(df[[var]])+0.1, y=max(df$elevationch)*0.85,
           "NS",
           cex=1.2, col = "red", pos = 4)
  }
  plot(df[[var]], df$slopem,
       ylim = range(df$slopech, df$slopecl),
       pch = 16, xlab = var, ylab = "SMA LTR-LCC-tradeoff slope")
  segments(df[[var]], df$slopecl, df[[var]], df$slopech)
  obj_lm <- lm(formula(paste("slopem ~ ", var)), data = df , weights = df$slopw)
  abline(obj_lm,
         col = 'gray')
  if(summary(obj_lm)$coefficients[2,4] < 0.01){
      text(x=min(df[[var]]), y=max(df$slopech)-0.2,
           paste0("* R2 =",round(summary(obj_lm)$r.squared, 3)),
           cex=1.2, col = "red", pos = 4)
  }else{
      text(x=min(df[[var]])+0.1, y=max(df$slopech)-0.2,
           "NS",
           cex=1.2, col = "red", pos = 4)
  }
}

figure_lma_tradeoff_climate_slope_elev<- function(data) {
  require(dplyr)
  data <- subset(data, !is.na(data[["lma"]] * data[["leaf_turnover"]])
    & table(data[["location"]])[data[["location"]]] > 5)
  location <- data[["location"]]
  lma <- data[["lma"]]/10^(log10(mean(data[["lma"]])))
  leaf_turnover <- data[["leaf_turnover"]]
  sm <- sma(leaf_turnover ~ lma, log="xy")
  print(summary(sm))
  sm1 <- sma(leaf_turnover ~ lma * location, log="xy")
  table_coef <- do.call("rbind",lapply(sm1$coef,
                         function(x){ d <- (cbind(x[1,], x[2, ]));
                                      colnames(d) <-  as.vector(t(outer(c('elevation', 'slope'),
                                                                      c('m', 'cl', 'ch'),
                                                                      paste0)));
                                      return(data.frame(d))}))
  df <- left_join(data.frame(location = rownames(table_coef),
                             table_coef,
                             pval = unlist(sm1$pval)),
                  data[!duplicated(data$location), c("location", "mat", "map", "mat_o_map", "ai")],
                  by = "location")
  df$elevw <- 1/(df$elevationch - df$elevationcl)
  df$slopw <- 1/(df$slopech - df$slopecl)
  df <-  df[df$pval <= 0.05, ]
  df$mat_o_map<-  scale(df$mat_o_map)
  df$map <-  scale(df$map)
  df$mat <-  scale(df$mat)

  par(mfrow = c(3,2), mar=c(2.5, 2.5, .5, .5), mgp = c(1.5, 0.5, 0))
  #MAP
  plot_coef_sma(df, "map")
  # MAP
  print("elevation vs map")
  print(summary(lm(elevationm~scale(map),
                   data = df ,
                   weights = df$slopw)))
  print("slope vs map over mat")
  print(summary(lm(slopem~scale(map),
                   data = df ,
                   weights = df$slopw)))

  #MAT/MAP
  plot_coef_sma(df, "mat_o_map")
  #MAT/MAP

  print("elevation vs map over mat")
  print(summary(lm(elevationm~scale(mat_o_map),
                   data = df ,
                   weights = df$slopw)))
  print("slope vs map over mat")
  print(summary(lm(slopem~scale(mat_o_map),
                   data = df ,
                   weights = df$slopw)))

  #ai
  plot_coef_sma(df, "ai")
  #ai

  print("elevation vs ai")
  print(summary(lm(elevationm~scale(ai),
                   data = df ,
                   weights = df$slopw)))
  print("slope vs ai")
  print(summary(lm(slopem~scale(ai),
                   data = df ,
                   weights = df$slopw)))


}



figure_B_kl_climate<- function(data) {
  require(dplyr)

  data <- subset(data, !is.na(data[["lma"]] * data[["leaf_turnover"]])
                       & table(data[["location"]])[data[["location"]]] > 5)
  location <- data[["location"]]
  lma <- data[["lma"]]/0.1978791
  leaf_turnover <- data[["leaf_turnover"]]

  sm <- sma(leaf_turnover ~ lma, log="xy")
  print(summary(sm))
  sm1 <- sma(leaf_turnover ~ lma * location, log="xy")
  table_coef <- do.call("rbind",lapply(sm1$coef,
                 function(x){ d <- (cbind(x[1,], x[2, ]));
                              colnames(d) <-  as.vector(t(outer(c('elevation',
                                                                  'slope'),
                                                                 c('m', 'cl',
                                                                   'ch'),
                                                                 paste0)));
                                      return(data.frame(d))}))
  df <- left_join(data.frame(location = rownames(table_coef),
                             table_coef,
                             pval = unlist(sm1$pval)),
                  data[!duplicated(data$location), c("location", "mat", "map")],
                  by = "location")
  df <-  df[df$pval <= 0.05, ]
  df$mat_o_map<-  scale(df$mat/df$map)
  df$mat<-  scale(df$mat)
  df$map<-  scale(df$map)
  df$elevw <- 1/(df$elevationch - df$elevationcl)
  df$slopw <- 1/(df$slopech - df$slopecl)

  param_P<- data.frame(coef = c("inter", "slope"),
                       elev = coef(lm(elevationm~map,
                                      data = df ,
                                      weights = df$elevw)),
                       slop = coef(lm(slopem~map,
                                     data = df ,
                                     weights = df$slopw))
                      )

pval_elev_P<- as.integer(summary(lm(elevationm~map,data = df,
                                    weights = df$elevw))$coefficients[2,4] > 0.01)
pval_slop_P<- as.integer(summary(lm(slopem~map,data = df,
                                    weights = df$slopw))$coefficients[2,4] > 0.01)


  param_TP<- data.frame(coef = c("inter", "slope"),
                      elev = coef(lm(elevationm~mat_o_map,
                                     data = df ,
                                     weights = df$elevw)),
                      slop = coef(lm(slopem~mat_o_map,
                                     data = df ,
                                     weights = df$slopw)))

pval_elev_TP<- as.integer(summary(lm(elevationm~mat_o_map,data = df,
                                    weights = df$elevw))$coefficients[2,4] > 0.01)
pval_slop_TP<- as.integer(summary(lm(slopem~mat_o_map,data = df,
                                    weights = df$slopw))$coefficients[2,4] > 0.01)

  seq_stress <- seq(from = -1.5, to = 1.5, length.out = 100)
  par(mfrow = c(2, 2))
  plot(seq_stress, 10^(param_P[1, 2] + seq_stress * param_P[2,2]),
       xlab = "MAP", ylab = "B_kl1", type = "l",
       xlim = range(df$map), ylim = range(10^df$elevationm),
       lty = pval_elev_P+1)
  points(df$map, 10^df$elevationm)
  segments(df$map, 10^df$elevationcl, df$map, 10^df$elevationch)
  abline(v=0, col = 'red')
  plot(seq_stress, param_P[1, 3] + seq_stress * param_P[2,3],
       xlab = "MAP", ylab = "B_kl2", type = "l",
       xlim = range(df$map), ylim = range(df$slopem),
       lty = pval_slop_P+1)
  points(df$map, df$slopem)
  segments(df$map, df$slopecl, df$map, df$slopech)
  abline(v=0, col = 'red')
  plot(seq_stress, 10^(param_TP[1, 2] + seq_stress * param_TP[2,2]),
       xlab = "MAT over MAP", ylab = "B_kl1", type = "l",
       xlim = range(df$mat_o_map), ylim = range(10^df$elevationm),
       lty = pval_elev_TP+1)
  points(df$mat_o_map, 10^df$elevationm)
  segments(df$mat_o_map, 10^df$elevationcl, df$mat_o_map, 10^df$elevationch)
  abline(v=0, col = 'red')
  plot(seq_stress, param_TP[1, 3] + seq_stress * param_TP[2,3],
       xlab = "MAT over MAP", ylab = "B_kl2", type = "l",
       xlim = range(df$mat_o_map), ylim = range(df$slopem),
       lty = pval_slop_TP+1)
  points(df$mat_o_map, df$slopem)
  segments(df$mat_o_map, df$slopecl, df$mat_o_map, df$slopech)
  abline(v=0, col = 'red')
}


param_B_kl_climate_P<- function(data) {
  require(dplyr)

  data <- subset(data, !is.na(data[["lma"]] * data[["leaf_turnover"]])
                       & table(data[["location"]])[data[["location"]]] > 5)
  location <- data[["location"]]
  lma <- data[["lma"]]/0.1978791
  leaf_turnover <- data[["leaf_turnover"]]

  sm <- sma(leaf_turnover ~ lma, log="xy")
  print(summary(sm))
  sm1 <- sma(leaf_turnover ~ lma * location, log="xy")
  table_coef <- do.call("rbind",lapply(sm1$coef,
                 function(x){ d <- (cbind(x[1,], x[2, ]));
                              colnames(d) <-  as.vector(t(outer(c('elevation',
                                                                  'slope'),
                                                                 c('m', 'cl',
                                                                   'ch'),
                                                                 paste0)));
                                      return(data.frame(d))}))
  df <- left_join(data.frame(location = rownames(table_coef),
                             table_coef,
                             pval = unlist(sm1$pval)),
                  data[!duplicated(data$location), c("location", "mat", "map")],
                  by = "location")
  df <-  df[df$pval <= 0.05, ]
  df$mat_o_map<-  scale(df$mat/df$map)
  df$mat<-  scale(df$mat)
  df$map<-  scale(df$map)
  df$elevw <- 1/(df$elevationch - df$elevationcl)
  df$slopw <- 1/(df$slopech - df$slopecl)

  param <- data.frame(coef = c("a", "b"),
                      LMAelev = coef(lm(elevationm~map,
                                     data = df, weights = df$elevm )),
                      LMAslope = coef(lm(slopem~map,
                                     data = df, weights = df$elevm )))

pval_elev_P<- summary(lm(elevationm~map,data = df,
                         weights = df$elevw))$coefficients[2,4] > 0.01
pval_slop_P<- summary(lm(slopem~map,data = df,
                         weights = df$slopw))$coefficients[2,4] > 0.01
if(pval_elev_P) param$LMAelev <-  NA
if(pval_slop_P) param$LMAslope <-  NA

  write.csv(param, file = "output/data_slope_P.csv", row.names = FALSE)
}

param_B_kl_climate_TP<- function(data) {
  require(dplyr)

  data <- subset(data, !is.na(data[["lma"]] * data[["leaf_turnover"]])
                       & table(data[["location"]])[data[["location"]]] > 9)
  location <- data[["location"]]
  lma <- data[["lma"]]/0.1978791
  leaf_turnover <- data[["leaf_turnover"]]

  sm <- sma(leaf_turnover ~ lma, log="xy")
  print(summary(sm))
  sm1 <- sma(leaf_turnover ~ lma * location, log="xy")
  table_coef <- do.call("rbind",lapply(sm1$coef,
                 function(x){ d <- (cbind(x[1,], x[2, ]));
                              colnames(d) <-  as.vector(t(outer(c('elevation',
                                                                  'slope'),
                                                                 c('m', 'cl',
                                                                   'ch'),
                                                                 paste0)));
                                      return(data.frame(d))}))
  df <- left_join(data.frame(location = rownames(table_coef),
                             table_coef,
                             pval = unlist(sm1$pval)),
                  data[!duplicated(data$location), c("location", "mat", "map", "mat_o_map")],
                  by = "location")
  df <-  df[df$pval <= 0.05, ]
  df$mat_o_map<-  scale(df$mat_o_map)
  df$mat<-  scale(df$mat)
  df$map<-  scale(df$map)
  df$elevw <- 1/(df$elevationch - df$elevationcl)
  df$slopw <- 1/(df$slopech - df$slopecl)

  param <- data.frame(coef = c("a", "b"),
                      LMAelev = coef(lm(elevationm~mat_o_map,
                                     data = df, , weights = df$slopm )),
                      LMAslope = coef(lm(slopem~mat_o_map,
                                     data = df , weights = df$slopm)))
pval_elev_TP<- summary(lm(elevationm~mat_o_map,data = df,
                         weights = df$elevw))$coefficients[2,4] > 0.01
pval_slop_TP<- summary(lm(slopem~mat_o_map,data = df,
                         weights = df$slopw))$coefficients[2,4] > 0.01
if(pval_elev_TP) param$LMAelev <-  NA
if(pval_slop_TP) param$LMAslope <-  NA

  write.csv(param, file = "output/data_slope_TP.csv", row.names = FALSE)
}
