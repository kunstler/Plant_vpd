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
  ## add VPD
  sitevars$vpd <- extract_vpd(sitevars$latitude, sitevars$longitude)


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
  names(data)[names(data) %in%
              c('dataset', 'mat_degc', 'map_mm')]<- c('location','mat', 'map')

  data$mat_o_map<-  (data$mat+18)/data$map

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
  vpd <- data[["vpd"]]

  myplot <- function(...,  xlabel=FALSE,  ylabel=FALSE) {
    plot(..., ann=FALSE, xaxt='n', yaxt='n', pch=16, cex=0.9)
    axis(1, labels=xlabel)
    axis(2, labels=ylabel, las=1)
  }

  mylabel <- function(lab, ...) label(0.05, 0.95, lab, log.y=TRUE, ...)

  par(mfrow=c(2,3), mar=c(1,1,1,1), oma=c(4,5,0,0))

  myplot(mat_o_map, lma, log ="xy",  xlabel=FALSE, ylabel=TRUE)
  mtext(expression(paste("Leaf-mass per area (kg ", m^-2,")")), 2, line=4)
  mylabel("a", log.x=TRUE)

  myplot(ai, lma, log ="xy", xlabel=FALSE)
  mylabel("b", log.x=TRUE)

  myplot(vpd, lma, log ="y", xlabel=FALSE)
  mylabel("c")

  myplot(mat_o_map, narea, log ="xy",xlabel=TRUE, ylabel=TRUE)
  mtext("Temp (C) over Precip (mm)", 1, line=4)
  mtext(expression(paste("Leaf-nitrogen per area (kg ", m^-2,")")), 2, line=4)
  mylabel("d", log.x=TRUE)

  myplot(ai, narea, log ="xy", xlabel=TRUE)
  mtext("Aridity Index MAP/PET", 1, line=4)
  mylabel("e", log.x=TRUE)

  myplot(vpd, narea, log ="y", xlabel=TRUE)
  mtext("Vapour pressure deficit", 1, line=4)
  mylabel("f")
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


figure_lma_tradeoff_narea <- function(data) {
    data <- subset(data,
                   !is.na(data[["n.area"]] * data[["lma"]] *
                          data[["leaf_turnover"]]))
  lma <- data[["lma"]]
  leaf_turnover <- data[["leaf_turnover"]]
  data$n.area<-  data$n.area/1.87e-3
  data$levels <- cut(data$n.area,
                     breaks = quantile(data$n.area,
                                       probs = seq(0,1, length.out = 7),
                                       na.rm = TRUE),
                     labels = FALSE, ordered_result = TRUE,
                     include.lowest = TRUE)
  groups<- data[['levels']]
  sm1 <- sma(leaf_turnover ~ lma * groups, log="xy")

  colfunc <- colorRampPalette(c("red", "blue"))
  cols <- colfunc(length(unique(data[['levels']])))
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
  legend("topright", legend=paste("Narea", " class",
                                  1:length(unique(data[['levels']]))), bty="n",
         pch=16, col=cols, cex=1, title=title)

}


figure_B_kl_narea<- function(data) {
  require(dplyr)
  data <- subset(data, !is.na(data[["n.area"]] * data[["lma"]] *
                              data[["leaf_turnover"]]))
  data$naream<-  data$n.area/1.87e-3
  data$levels <- cut(data$naream,
                   breaks = quantile(data$naream,
                                     probs = seq(0,1, length.out = 15),
                                     na.rm = TRUE),
                   labels = FALSE, ordered_result = TRUE,
                   include.lowest = TRUE)
  levels <- data$levels
  data_summ <- data %>% group_by(levels) %>%
      summarise(n.area.mean = mean(naream)) %>% ungroup()
  lma <- data[["lma"]]/0.1978791
  leaf_turnover <- data[["leaf_turnover"]]
  sm <- sma(leaf_turnover ~ lma, log="xy")
  print(summary(sm))
  sm1 <- sma(leaf_turnover ~ lma * levels, log="xy")
  table_coef <- do.call("rbind",lapply(sm1$coef,
                         function(x){ d <- (cbind(x[1,], x[2, ]));
                             colnames(d) <-  as.vector(t(outer(c('elevation',
                                                                 'slope'),
                                                               c('m', 'cl',
                                                                 'ch'),
                                                                      paste0)));
                                      return(data.frame(d))}))
  df <- left_join(data.frame(levels = as.integer(rownames(table_coef)),
                             table_coef,
                             pval = unlist(sm1$pval)),
                  data_summ,
                  by = "levels")
  df$elevationm <- 10^df$elevationm
  df$elevationch <- 10^df$elevationch
  df$elevationcl <- 10^df$elevationcl
  df$elevw <- 1/(df$elevationch - df$elevationcl)
  df$slopw <- 1/(df$slopech - df$slopecl)
  df <-  df[df$pval <= 0.05, ]
  df$naream<-  df$n.area.mean
  param_N<- data.frame(coef = c("inter", "slope"),
                       elev = coef(lm(elevationm~naream,
                                      data = df ,
                                      weights = df$elevw)),
                       slop = coef(lm(slopem~naream,
                                     data = df ,
                                     weights = df$slopw))
                       )
  print(param_N)
  pval_elev_N<- as.integer(summary(lm(elevationm~naream, data = df,
                                weights = df$elevw))$coefficients[2,4] > 0.01)
  pval_slop_N<- as.integer(summary(lm(slopem~naream, data = df,
                                weights = df$slopw))$coefficients[2,4] > 0.01)
  seq_stress <- seq(from = min(df$naream), to = max(df$naream),
                    length.out = 100)
  par(mfrow = c(1, 2))
  plot(seq_stress, param_N[1, 2] + seq_stress * param_N[2,2],
       xlab = "Narea over Narea mean", ylab = "B_kl1", type = "l",
       xlim = range(df$naream),
       ylim = range(df$elevationm, df$elevationch, df$elevationcl),
       lty = pval_elev_N+1)
  points(df$naream, df$elevationm)
  segments(df$naream, df$elevationcl, df$naream, df$elevationch)
  abline(v=1, col = 'red')
  plot(seq_stress, param_N[1, 3] + seq_stress * param_N[2,3],
       xlab = "Narea over Narea mean", ylab = "B_kl2", type = "l",
       xlim = range(df$naream),
       ylim = range(df$slopem, df$slopech, df$slopecl),
       lty = pval_slop_N+1)
  points(df$naream, df$slopem)
  segments(df$naream, df$slopecl, df$naream, df$slopech)
  abline(v=1, col = 'red')
}




param_B_kl_narea<- function(data) {
  require(dplyr)

  data <- subset(data, !is.na(data[["lma"]] * data[["leaf_turnover"]] * data[["n.area"]]))
  data$naream<-  data$n.area/1.87e-3
  data$levels <- cut(data$naream,
                   breaks = quantile(data$naream,
                                     probs = seq(0,1, length.out = 15),
                                     na.rm = TRUE),
                   labels = FALSE, ordered_result = TRUE,
                   include.lowest = TRUE)
  data_summ <- data %>% group_by(levels) %>% summarise(n.area.mean = mean(naream)) %>% ungroup()
  levels <- data$levels
  lma <- data[["lma"]]/0.1978791
  leaf_turnover <- data[["leaf_turnover"]]

  sm <- sma(leaf_turnover ~ lma, log="xy")
  print(summary(sm))
  sm1 <- sma(leaf_turnover ~ lma * levels, log="xy")
  table_coef <- do.call("rbind",lapply(sm1$coef,
                 function(x){ d <- (cbind(x[1,], x[2, ]));
                              colnames(d) <-  as.vector(t(outer(c('elevation',
                                                                  'slope'),
                                                                 c('m', 'cl',
                                                                   'ch'),
                                                                 paste0)));
                                      return(data.frame(d))}))
  df <- left_join(data.frame(levels = as.integer(rownames(table_coef)),
                             table_coef,
                             pval = unlist(sm1$pval)),
                  data_summ,
                  by = "levels")
  
  df$elevationm <- 10^df$elevationm
  df$elevationcl <- 10^df$elevationcl
  df$elevationch <- 10^df$elevationch
  df$elevw <- 1/(df$elevationch - df$elevationcl)
  df$slopw <- 1/(df$slopech - df$slopecl)
  df <-  df[df$pval <= 0.05, ]
  df$naream<-  df$n.area.mean
  
  param <- data.frame(coef = c("a", "b"),
                      LMAelev = coef(lm(elevationm~naream,
                                     data = df, , weights = df$elevw )),
                      LMAslope = coef(lm(slopem~naream,
                                     data = df , weights = df$slopw)))
pval_elev_TP<- summary(lm(elevationm~naream,data = df,
                         weights = df$elevw))$coefficients[2,4] > 0.01
pval_slop_TP<- summary(lm(slopem~naream,data = df,
                         weights = df$slopw))$coefficients[2,4] > 0.01
if(pval_elev_TP) param$LMAelev <-  NA
if(pval_slop_TP) param$LMAslope <-  NA

  write.csv(param, file = "output/data_slope_narea.csv", row.names = FALSE)
}

