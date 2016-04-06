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

## get climatic niche of the species

get_climate_niche_sp<- function(sp){
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

    if(nrow(dat)>0){
      #gbifmap(dat)
      dat <- dat[!is.na(dat$decimalLatitude) & dat$decimalLatitude != 0 &
                 !is.na(dat$decimalLongitude) & dat$decimalLongitude != 0, ]
      df_clim <- GetClimate(dat$decimalLatitude, dat$decimalLongitude)

        res <- data.frame(Binomial= sp,
                          do.call("cbind", lapply(df_clim, mean)),
                          nobs = nrow(dat))
      print(paste("species ", sp, " done"))
    }else{
      res  <- data.frame(Binomial= sp,
                         MAT = NA, MAP = NA, TMAX = NA, PDRY = NA,
                         nobs = nrow(dat))
      print(paste("species ", sp, " error wc"))
    }
}else{
    res  <- data.frame(Binomial= sp,
                       MAT = NA, MAP = NA, TMAX = NA, PDRY = NA,
                       nobs = nrow(dat))
    print(paste("species ", sp, " errro gbif"))
}
return(res)
}

get_climate_niche<- function(df){
require(dplyr)
res <- rbind_list(lapply(df$Binomial, get_climate_niche_sp))
rr <- left_join(df, res, by = "Binomial")
return(rr)
}


