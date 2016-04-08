process_zanne_2010<- function(filename){
  require(xlsx)
  ## There are several strategies for reading in an excel file, but
  ## this one works quite well.
  df <-read.xlsx2(filename, sheetIndex=2,
                  stringsAsFactors=FALSE, check.names=FALSE)
  # remove empty row (no idea why this happen with the xls file ...
  df <-  df[df$Family != '', ]
  names(df) <- gsub('_$', '',
                    gsub("__", "_",
                         gsub("[ '('')''/''^'-]", "_",
                              names(df))))
  df$Genus <- sub(" .*", "", df$Binomial)
  df$Species <- sub(".* ", "", df$Binomial)
  for (i in names(df)[3:ncol(df)]) {df[[i]] <-  as.numeric(df[[i]])}
  return(df)
}


process_olson_2014<- function(filename){
  require(xlsx)
  ## There are several strategies for reading in an excel file, but
  ## this one works quite well.
  df <-read.xlsx2(filename, sheetIndex=2,
                  stringsAsFactors=FALSE, check.names=FALSE)

  names(df) <- paste(names(df), df[1, ])
  df <-  df[-1, ]

  for(i in c(6:11, 14:33)) df[, i] <- as.numeric(df[, i])
  return(df)
}


## get climatic niche of the species

get_climate_niche_sp<- function(sp, wc){
require(dismo)
require(rgbif)
##data cleaning TODO IMPROVE JUST A FIRST GUESS
code_issues <- c("cdout", "cdrepf", "cucdmis", "preswcd", "zerocd")
key <- name_suggest(q= sp, rank='species')$key[1]
dat <- tryCatch(occ_search(taxonKey=key, return='data',
                           fields =c('name', 'decimalLatitude',
                                     'decimalLongitude', 'issues')),
                error = function(err) NA)
if(is.data.frame(dat)){
    dat <-  dat %>% occ_issues(-cdout, -cdrepf, -cucdmis, -preswcd, -zerocd)
    dat <- dat[!is.na(dat$decimalLatitude) & dat$decimalLatitude != 0 &
               !is.na(dat$decimalLongitude) & dat$decimalLongitude != 0, ]

   if(nrow(dat)>0){
      #gbifmap(dat)
      df_clim <- GetClimate(dat$decimalLatitude, dat$decimalLongitude, wc)
      if(nrow(df_clim)>0){
        df_wc <- lapply(df_clim, mean, na.rm = TRUE)
        res <- data.frame(Binomial= sp,
                          do.call("cbind", df_wc),
                          nobs = nrow(dat), stringsAsFactors = FALSE)
        names(res) <-  c("Binomial", "MAT", "MAP", "TMAX", "PDRY", "nobs")
        print(paste("species ", sp, " done"))
      }else{
        res  <- data.frame(Binomial= sp,
                           MAT = NA, MAP = NA, TMAX = NA, PDRY = NA,
                           nobs = 0, stringsAsFactors = FALSE)
        print(paste("species ", sp, " error in wc no data with climate"))
     }
   }else{
    res  <- data.frame(Binomial= sp,
                       MAT = NA, MAP = NA, TMAX = NA, PDRY = NA,
                       nobs = 0, stringsAsFactors = FALSE)
    print(paste("species ", sp, " error in gbif no data with coordinates"))
   }
}else{
res  <- data.frame(Binomial= sp,
                   MAT = NA, MAP = NA, TMAX = NA, PDRY = NA,
                   nobs = 0, stringsAsFactors = FALSE)
print(paste("species ", sp, " error in gbif no data"))
}
return(res)
}


get_climate_niche<- function(df, wc){
require(dplyr)
require(parallel)
cl <- makeCluster(20)
list_res <- parLapply(cl, df$Binomial, get_climate_niche_sp, wc = wc)
res <- dplyr::rbind_all(list_res)
stopCluster(cl)
rr <- left_join(df, res, by = "Binomial")
return(rr)
}
##


plot_sma <- function(df, var_x, var_y, log="xy"){
  sm1 <- sma(df[[var_y]] ~ df[[var_x]], log=log)
  par(mar=c(4.6, 4.6, .5, .5))
  plot(NA, type="n", log=log, xlim=range(df[[var_x]], na.rm = TRUE),
       ylim= range(df[[var_y]], na.rm = TRUE),
       xlab="", ylab="", las=1)
  mtext(var_x, line=3, side = 1)
  mtext(var_y, line=3, side = 2)
  points(df[[var_x]], df[[var_y]], col=make_transparent("grey", 0.5), pch=16)
  plot(sm1, add=TRUE, type="l", lwd=1.25, p.lines.transparent=NA)
}

plot_wood_clim <- function(df){
require(smatr)
df <-  df[!is.na(df$MAP), ]
df$PDRY <- df$PDRY + 0.006
df$P_T <- df$PDRY/df$MAT

par(mfrow = c(2,3))
plot_sma(df, var_x = "MAP", var_y = "A_mm_2")
plot_sma(df, var_x = "MAP", var_y = "F_mm_2_mm_2")
plot_sma(df, var_x = "MAP", var_y = "N_mm_2")

plot_sma(df, var_x = "PDRY", var_y = "A_mm_2")
plot_sma(df, var_x = "PDRY", var_y = "F_mm_2_mm_2")
plot_sma(df, var_x = "PDRY", var_y = "N_mm_2")


}
