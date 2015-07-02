
figure_lma_climate <- function(data) {
  data <- subset(data, data[["map"]] > 400
    & data[["growthform"]] %in% c("S","T"))

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
  mtext(expression(paste("Leaf-construction cost (kg ", m^-2,")")), 2, line=4)
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

  points(lma, leaf_turnover, col=make_transparent("grey", 0.3))
  plot(sm1, add=TRUE, col=col_sm1, type="l", lwd=1.25, p.lines.transparent=0.15)

  x <- seq_log_range(c(0.001,3), 40)
  points(x, 0.0286*x^-1.71, type='l', col='black', lwd=2)

  title <- sprintf("%d sites, %d species",
                   length(unique(location)),
                   sum(!is.na(leaf_turnover)))
  legend("topright", legend=tolower(names(col_table)), bty="n",
         pch=16, col=col_table, cex=1, title=title)
}
