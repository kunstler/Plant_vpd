process_chave_2009 <- function(filename) {
  ## There are several strategies for reading in an excel file, but this one works
  ## quite well.
  df <- read.xlsx2(filename, sheetIndex = 2, stringsAsFactors = FALSE, check.names = FALSE)
  # remove empty row (no idea why this happen with the xls file ...
  df <- df[df$Family != "", ]
  names(df) <- gsub("_$", "", gsub("__", "_", gsub("[ '('')''/''^'-]", "_", names(df))))
  df$Genus <- sub(" .*", "", df$Binomial)
  df$Species <- sub(".* ", "", df$Binomial)
  names(df) <- c("Number", "Family", "Binomial", "Wood_density",
                 "Region", "Reference_Number", "Genus", "Species")
  for (i in names(df)[c(1, 4, 6)]) {
    df[[i]] <- as.numeric(df[[i]])
  }
  #compute mean WD per species
  require(dplyr)
  df_sp <-  group_by(df, Binomial) %>% mutate(Wood_density = mean(Wood_density, na.rm = TRUE),
                                                  count = n()) %>% slice(1) %>%  ungroup()
  df_sp
}

plot_wood_density_clim <- function(df) {

  df <- df[!is.na(df$MAP), ]
  df <- df[df$MAP>0, ]
  df <- df[!is.na(df$MAT), ]
  df$PDRY <- df$PDRY + 0.006
  df$P_T <- df$PDRY/df$MAT
  mylabel <- function(lab, ...) label(0.05, 0.95, lab, log.y = TRUE, ...)
## TODO FIX LOG PROBLEM
  par(mfrow=c(1,2), mar=c(1,1,1,1), oma=c(4,5,0,0))
  plot_sma(df, var_x = "MAT", var_y = "Wood_density", ylabel = TRUE, xlabel = TRUE, log ="y")
  mtext(expression(paste("Wood density")), 2, line=4)
  mtext(expression(paste("Av. temperature (",degree,"C)")) ,1, line=4)
  mylabel("a")
  plot_sma(df, var_x = "MAP", var_y = "Wood_density", xlabel = TRUE)
  mylabel("b", log.x = TRUE)
  mtext(expression(paste("Precipitation (mm)")), 1, line=4)
}
