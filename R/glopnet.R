process_wright_2004 <- function(filename, sitevars_file) {

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
  sitevars <- read.csv(sitevars_file, stringsAsFactors=FALSE)
  data <- merge(d, sitevars, by.x = 'Dataset', by.y = 'dataset_location')

  #lowercase names
  names(data) <- tolower(names(data))

  # unit conversions
  data$lma <- data$lma/1000 # Converts to kg
  data$n.area <- data$n.area/1000 # Converts to kg
  data$a.area <- (data$a.area * 31557.6)*10^-6 # converts to mol/kg/yr from micro-mol/g/s
  data$rd.area <- (data$rd.area * 31557.6)*10^-6 # converts to mol/kg/yr from micro-mol/g/s

  data$leaflifespan <- data$leaflifespan/12 ## convert LL from months to years
  data$leaf_turnover <- 1/data$leaflifespan ## per year
  names(data)[names(data) %in% c('dataset', 'mat_degc', 'map_mm')]<- c('location','mat', 'map')

  subset(data, data[["map"]] > 400
    & data[["growthform"]] %in% c("S","T"))
}

figure_lma_climate <- function(data) {

  lma <- data[["lma"]]
  narea <- data[["n.area"]]
  map <- data[["map"]]
  mat <- data[["mat"]]

  grey <- make_transparent("grey", 0.3)
  blues <- rev(brewer.pal(8, "Blues"))
  greens <- rev(brewer.pal(8, "Greens"))

  bins <- seq_log_range(range(map), 6)
  map_bin <- as.numeric(
              cut(map, bins, seq_len(length(bins)-1), include.lowest = TRUE))

  bins <- seq_range(range(mat), 6)
  mat_bin <- as.numeric(
              cut(mat, bins, seq_len(length(bins)-1), include.lowest = TRUE))


  myplot <- function(...,  xlabel=FALSE,  ylabel=FALSE) {
    plot(..., ann=FALSE, xaxt='n', yaxt='n', pch=16, cex=0.9)
    axis(1, labels=xlabel)
    axis(2, labels=ylabel, las=1)
  }

  mylabel <- function(lab, ...) label(0.05, 0.95, lab, log.y=TRUE, ...)

  par(mfrow=c(2,2), mar=c(1,1,1,1), oma=c(4,5,0,0))

  myplot(mat, lma, log ="y", col=greens[map_bin], xlabel=FALSE, ylabel=TRUE)
  mtext(expression(paste("Leaf-mass per area (kg ", m^-2,")")), 2, line=4)
  mylabel("a")

  myplot(map, lma, log ="xy", col=blues[mat_bin], xlabel=FALSE)
  mylabel("b", log.x=TRUE)

  myplot(mat, narea, log ="y", col=greens[map_bin], xlabel=TRUE, ylabel=TRUE)
  mtext("Av. Temperature (deg C)",1, line=4)
  mtext(expression(paste("Leaf-nitrogen per area (kg ", m^-2,")")), 2, line=4)
  mylabel("c")

  myplot(map, narea, log ="xy", col=blues[mat_bin], xlabel=TRUE)
  mtext("Precipitation (mm)", 1, line=4)
  mylabel("d", log.x=TRUE)
}

figure_lma_tradeoff <- function(data) {

  data <- subset(data, !is.na(data[["lma"]] * data[["leaf_turnover"]])
    & table(data[["location"]])[data[["location"]]] > 5)

  location <- data[["location"]]
  lma <- data[["lma"]]
  leaf_turnover <- data[["leaf_turnover"]]
  biomes <- sort(unique(data[["biome"]]))

  col_table <- nice_colors(length(biomes))
  names(col_table) <- biomes

  sm1 <- sma(leaf_turnover ~ lma * location, log="xy")
  col_sm1 <- col_table[data[["biome"]][match(sm1[["groups"]], location)]]

  par(mar=c(4.6, 4.6, .5, .5))
  plot(NA, type="n", log="xy", xlim=c(0.01, 1.28), ylim=c(0.03, 32),
       xlab="", ylab="", las=1)
  mtext(expression(paste("Leaf-construction cost (kg ", m^-2,")")), line=3, side = 1)
  mtext(expression(paste("Leaf turnover rate (",yr^-1,")")), line=3, side = 2)

  points(lma, leaf_turnover, col=make_transparent("grey", 0.5), pch=16)
  plot(sm1, add=TRUE, col=col_sm1, type="l", lwd=1.25, p.lines.transparent=0.15)

  x <- seq_log_range(c(0.001,3), 40)
  points(x, 0.0286*x^-1.71, type='l', col='black', lwd=2)

  title <- sprintf("%d sites, %d species",
                   length(unique(location)),
                   sum(!is.na(leaf_turnover)))
  legend("topright", legend=tolower(names(col_table)), bty="n",
         pch=16, col=col_table, cex=1, title=title)
}
