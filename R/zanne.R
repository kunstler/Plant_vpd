process_zanne_2010 <- function(filename) {
  ## There are several strategies for reading in an excel file, but this one works
  ## quite well.
  df <- read.xlsx2(filename, sheetIndex = 2,
                   stringsAsFactors = FALSE, check.names = FALSE)
  # remove empty row (no idea why this happen with the xls file ...
  df <- df[df$Family != "", ]
  names(df) <- gsub("_$", "", gsub("__", "_", gsub("[ '('')''/''^'-]", "_", names(df))))
  df$Genus <- sub(" .*", "", df$Binomial)
  df$Species <- sub(".* ", "", df$Binomial)
  for (i in names(df)[3:ncol(df)]) {
    df[[i]] <- as.numeric(df[[i]])
  }
  df
}


process_olson_2014 <- function(filename) {
  ## There are several strategies for reading in an excel file, but this one works
  ## quite well.
  df <- read.xlsx2(filename, sheetIndex = 2,
                   stringsAsFactors = FALSE, check.names = FALSE)

  names(df) <- paste(names(df), df[1, ])
  df <- df[-1, ]

  for (i in c(6:11, 14:33)) df[, i] <- as.numeric(df[, i])
  df
}

xy_extract_wc <- function(dat, wc, sup_lim, sp){
      # TODO BJOERN COMMENT ACCOUNT FOR SAPLING BIAS IN GBIF, for one simple example see http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0055158 use densitylist to get teh sampling density
      df_clim <- GetClimate(dat$decimalLatitude, dat$decimalLongitude, wc)
      if (nrow(df_clim) > 0) {
        df_wc <- lapply(df_clim, mean, na.rm = TRUE)
        ret <- data.frame(Binomial = sp, do.call("cbind", df_wc), nobs = nrow(dat),
                          sup_lim = sup_lim, stringsAsFactors = FALSE)
        names(ret) <- c("Binomial", "MAT", "MAP", "TMAX", "PDRY", "nobs", "sup_lim")
        print(paste("species ", sp, " done", nrow(dat)))
      } else {
        print(paste("species ", sp, " error in wc no data with climate"))
      }
return(ret)
}

clean_gbif_data <- function(dat){
 if(sum(names(dat) %in% "locality") == 0) dat$locality <- NA
 dat <- dat %>% filter(!issues %in% c("cdout", "cdrepf", "cucdmis", "preswcd", "zerocd",
                            "bri", "ccm", "conti", "cdreps", "cuiv", "cum",
                            "gdativ", "preneglat", "preneglon",
                            "txmatnon"),
                       !is.na(decimalLatitude) & decimalLatitude != 0,
                       !is.na(decimalLongitude) & decimalLongitude != 0,
                       !is.na(locality))
 dat <- dplyr::filter(dat, basisOfRecord %in% c("HUMAN_OBSERVATION", "OBSERVATION",
                                         "LIVING_SPECIMEN", "PRESERVED_SPECIMEN"))

# TRY TO REMOVE BOTANICAL GARDENS see https://github.com/ropensci/rgbif/issues/110
 loc_to_remove <- c("Jardín", "Jardin", "Jardí", "Jardim", "Garden",
                    "Garten", "Orto", "Hortu", "Tuin", "Giardino", "Ogród",
                    "Musée", "Museum", "Muséum", "Museo", "Museu",
                    "Arboretum", "Arboreto", "Arboretos",
                    "Herbarium", "Herbario", "Herbier", "Erbario", "Erbario",
                    "Herbário", "Herbários",
                    "Yale Natural",
                    "University", "Université", "Universität", "Universidad",
                    "Università", "Universidade","Universiteit")

 remove_obs <- unique(grep(paste(loc_to_remove,collapse="|"),
                     dat$locality, ignore.case = TRUE))
 if(length(remove_obs) > 0) dat <- dat[-remove_obs, ]
return(dat)
}

## get climatic niche of the species
get_climate_niche_sp <- function(sp, wc) {
  require(dismo)
  require(rgbif)
  ret <- data.frame(Binomial = sp, MAT = NA, MAP = NA, TMAX = NA,
          PDRY = NA, nobs = 0, sup_lim = NA, stringsAsFactors = FALSE)
  sup_lim <- NA
  dat <- NA
  try({dat <- occ_data(scientificName = sp, minimal = FALSE,
                      hasCoordinate = TRUE,
                      limit = 200000, hasGeospatialIssue = FALSE)$data}
      )
  if (is.data.frame(dat)) {
    try(if(nrow(dat)>199999) sup_lim <- "sup")
    try(dat <- clean_gbif_data(dat))
      if (nrow(dat) > 0) {
       try(ret <- xy_extract_wc(dat, wc, sup_lim, sp))
      } else {
      print(paste("species ", sp, " error in gbif no data with coordinates"))
      }
   } else {
    print(paste("species ", sp, " error in gbif no data"))
   }
  ## Sys.sleep(2)
  return(ret)
}

#################
### TODO FOR SPECIES WITH MORE THAN 200000 occurence need to use occ_download
## dir.create(file.path("gbif"), showWarnings = FALSE)
## key <- name_backbone(name="Encelia californica")$speciesKey
## occ_count(taxonKey = key)


## key <- name_backbone(name="Abelia biflora")$speciesKey
## nobs <- occ_count(taxonKey = key, georeferenced = TRUE)
## if (nobs > 200000){
## res <- occ_download( paste("taxonKey =", key),
##                      user = "georges_kunstler",
##                      pwd = "gbif_traits",
##                      email = "georges.kunstler@gmail.com" )
## dat <- occ_download_get(res[[1]], path = "gbif", overwrite = TRUE)  %>% occ_download_import
## }

##     dat <- dat %>% filter(!issue %in% c("cdout", "cdrepf", "cucdmis", "preswcd", "zerocd",
##                                "bri", "ccm", "conti", "cdreps", "cuiv", "cum",
##                                "gdativ", "preneglat", "preneglon",
##                                "txmatnon"),
##                           !is.na(decimalLatitude) & decimalLatitude != 0,
##                           !is.na(decimalLongitude) & decimalLongitude != 0,
##                           !is.na(locality))
##     dat <- dplyr::filter(dat, basisOfRecord %in% c("HUMAN_OBSERVATION", "OBSERVATION",
##                                             "LIVING_SPECIMEN", "PRESERVED_SPECIMEN"))

get_climate_niche <- function(df, wc) {
  list_res <- lapply(df$Binomial, get_climate_niche_sp, wc = wc)
  print("gbif extraction done")
  res <- dplyr::bind_rows(list_res)
  df_out <- left_join(df, res, by = "Binomial")
  return(df_out)
}


plot_sma <- function(df, var_x, var_y, log = "xy", xlabel = FALSE, ylabel = FALSE) {
  sm1 <- sma(df[[var_y]] ~ df[[var_x]], log = log)
  plot(NA, type = "n", log = log, ann = FALSE,
       xlim = range(df[[var_x]], na.rm = TRUE),
       ylim = range(df[[var_y]], na.rm = TRUE),
       xaxt = "n", yaxt = "n", pch = 16, cex = 0.9)
    axis(1, labels=xlabel)
    axis(2, labels=ylabel, las=1)
  points(df[[var_x]], df[[var_y]], col = make_transparent("grey", 0.5), pch = 16)
  plot.sma(sm1, add = TRUE, type = "l", lwd = 1.25, p.lines.transparent = NA)
}

plot_vessel_clim <- function(df) {
  df <- df[!is.na(df$MAP), ]
  df <- df[df$MAP>0, ]
  df <- df[!is.na(df$MAT), ]
  df$PDRY <- df$PDRY + 0.006
  df$P_T <- df$PDRY/df$MAT
  mylabel <- function(lab, ...) label(0.05, 0.95, lab, log.y = TRUE, ...)
## TODO FIX LOG PROBLEM
  par(mfrow=c(3,2), mar=c(1,1,1,1), oma=c(4,5,0,0))
  plot_sma(df, var_x = "MAT", var_y = "A_mm_2", ylabel = TRUE, log ="y")
  mtext(expression(paste("Vessel area, A(",mm^2,")")), 2, line=4)
  mylabel("a")
  plot_sma(df, var_x = "MAP", var_y = "A_mm_2")
  mylabel("b", log.x = TRUE)
  plot_sma(df, var_x = "MAT", var_y = "F_mm_2_mm_2", ylabel = TRUE, log ="y")
  mtext(expression(paste("Lumen fraction")), 2, line=4)
  mylabel("c")
  plot_sma(df, var_x = "MAP", var_y = "F_mm_2_mm_2")
  mylabel("d", log.x = TRUE)
  plot_sma(df, var_x = "MAT", var_y = "N_mm_2", ylabel = TRUE, xlabel = TRUE, log ="y")
  mtext(expression(paste("Av. temperature (",degree,"C)")) ,1, line=4)
  mtext(expression(paste("Vessel number, N(",mm^{-2}, ")")), 2, line=4)
  mylabel("e")
  plot_sma(df, var_x = "MAP", var_y = "N_mm_2", xlabel = TRUE)
  mtext(expression(paste("Precipitation (mm)")), 1, line=4)
  mylabel("f", log.x = TRUE)
}




plot.sma <- function(x, which=c("default","residual","qq"),  use.null=FALSE, add=FALSE, type='o',
	xaxis=NULL, yaxis=NULL, xlab=NULL, ylab=NULL, pch=NULL, col=NULL, lty=NULL,
        from=NULL, to = NULL, log=x$log,
	frame.plot = TRUE, tck=par("tck"),p.lines.transparent=NA, ...){

	# function used to make colours transparent alpha = 0 means fully transparaent
	make.transparent <- function(col, alpha=1) {
  		tmp <- col2rgb(col)/255
	rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
	}

	#preprocessing ------------------------------------------------------------
	obj <- x  # this needed for consistency with plot
	if(obj$gt == "none"){
		ngrps <- 1
	}
	else{
		groups <- levels(obj$data[,3])
		ngrps <- length(groups)
	}

	whichplot <- match.arg(which)

	#---colors--------------------------------
	#user-defined colors
	if(!is.null(col[1])){
		if(length(col)== 1 &&  ngrps > 1)
			col<-rep(col[1],ngrps); #check right vector length
	} else {
	#default colors
		col <- c("blue2",  "goldenrod1", "firebrick2", "chartreuse4", "deepskyblue1", "darkorange1",
		"darkorchid3", "darkgrey", "mediumpurple1", "orangered2", "chocolate", "burlywood3",
		"goldenrod4", "darkolivegreen2", "palevioletred3", "darkseagreen3", "sandybrown", "tan",
		"gold", "violetred4", "darkgreen")
		col <- rep(col, ceiling(ngrps/length(col)))
	}

	#---symbols--------------------------------
	#user-defined symbols
	if(!is.null(pch[1])){
		if(length(pch) == 1 && ngrps > 1)
			pch <- rep(pch[1],ngrps) #check right vector length
	} else #default SYMBOLS
		pch <- rep(1,ngrps)
	#---line type--------------------------------	#user-defined symbol
	if(!is.null(lty[1])){
			if(length(lty) == 1 && ngrps > 1)
			lty<-rep(lty[1],ngrps) #check right vector length
	} else #default SYMBOLS
		lty <- rep("solid", ngrps)
	#-----------------------------------------------------------------------------
	# DATA PLOT
	if(whichplot == "default"){

		 #obtain data--------------------------------
		 Y <- obj$data[,1]
		 X <- obj$data[,2]

	    #log trasnformations--------------------------------
	    log_data <- obj$log	#logging of data based on transformation applied in sma. Allows scaling of axes ot be changed, while maintaing correct transformation of data
		XLog <- YLog <- 0
		if((log_data == "x") | (log_data == "xy")){ XLog=1; X = 10^X}
		if((log_data == "y") | (log_data == "xy")){ YLog=1; Y = 10^Y}

		#axis labels--------------------------------
	    #determine axis labels if not already specified
	    if(is.null(xlab)){
	    	xlab <- names(obj$data)[2]
			if(XLog)
				xlab <- paste(xlab, "[log scale]")
		}
	    if(is.null(ylab)){
	    	ylab <- names(obj$data)[1]
			if(YLog)
				ylab <- paste(ylab, "[log scale]")
		}

		#SETUP AXES--------------------------------
		if(add==FALSE)
		{
			#use nice plotting if appropriate
			if(!is.null(xaxis)  && !is.null(yaxis)){

				#Deteremine axis limits if missing - caluclated on transformed data. 				#add 5% white space around plots, then back transform if necessary
				if (is.null(xaxis$limits)){
					Range_x <-range(obj$data[,2])
					Buffer <- (Range_x[2]-Range_x[1])*0.05
					xaxis$limits <- c(Range_x[1]-Buffer, Range_x[2]+Buffer)
					if(XLog)xaxis$limits <- 10^xaxis$limits
				}
				if (is.null(yaxis$limits)){
					Range_y <-range(obj$data[,1])
					#add 4% white space around plots (like R default in plot)
					Buffer <- (Range_y[2]-Range_y[1])*0.04
					yaxis$limits <- c(Range_y[1]-Buffer, Range_y[2]+Buffer)
					if(YLog) yaxis$limits <- 10^yaxis$limits
				}

				#Make plot
				nicePlot(xaxis,yaxis,log=log,xlab=xlab, ylab=ylab,
					frame.plot = frame.plot, tck=tck,...)
			}
			else
				plot(X,Y, type='n', log=log, xlab=xlab, ylab=ylab,
					frame.plot = frame.plot, tck=tck,...)
		}

		#add datapoints	--------------------------------
		if(type %in% c("o","b", "p")){
			if(obj$gt == "none")
				points(X, Y, col = col[1], pch=pch[1],...)
			else{
				for(i in 1:ngrps){
					iref  <- as.character(obj$data[,3]) == groups[i]
					points(X[iref], Y[iref], col =col[i], pch=pch[i],...)
				}
			}
		}

		#add lines --------------------------------
		if(type %in% c("o","b", "l")){

			#decide which coefficients to use: alternative (default) or null
			if(use.null==FALSE)
				coef <- obj$coef
			else
				coef <- obj$nullcoef

			#determine end points for lines
			if(is.null(from[1])){  #based on fitted values
				for(i in 1:ngrps){
					from[i] <- as.numeric(obj$from[i])
					to[i] <- as.numeric(obj$to[i])
				}
			} else {  #user defined
				if(length(from) == 1){
					from <- rep(from[1], ngrps)
					to <- rep(to[1], ngrps)
				}
			}

			#add lines to plot
			for(i in 1:ngrps){
				#coefficients
				a <- coef[[i]][1,1]
				B <-  coef[[i]][2,1]

			    p <- obj$groupsummary$p[i]

				if(!is.na(p.lines.transparent))
					col.tr <- make.transparent(col[i], max(0, (1 - p/p.lines.transparent)))
				else
					col.tr <-  col[i]

#choose line according to log-trsnaformation used in fitting data, even if different transformation used for axes
				if(log_data=="xy")
	        		curve(10^a*x^B, from[i], to[i], add=TRUE,col = col.tr, lty= lty[i],...)
    	   	 	if(log_data=="x")
        			curve(a+B*log10(x), from[i], to[i], add=TRUE, col = col.tr, lty= lty[i],...)
        		if(log_data=="y")
        			curve(10^(a+x*B), from[i], to[i], add=TRUE, col = col.tr, lty= lty[i],...)
       	    	if(log_data=="")
        			curve(a + x*B, from[i], to[i], add=TRUE,  col = col.tr, lty= lty[i],...)
        	}
        }
	}

	# RESIDUAL PLOT
	if(whichplot == "residual")
	{
		 #obtain data--------------------------------
		Y <- fitted.sma(obj, type = "residuals")
		X <- fitted.sma(obj, type = "fitted")


		#axis labels--------------------------------
	    #determine axis labels if not already specified
	    if(is.null(xlab)) xlab <- paste("Fitted values (",names(obj$data)[2], " v ",
			names(obj$data)[1],")")
	    if(is.null(ylab)) ylab <- paste("Residual values (",names(obj$data)[2], " v ",
			names(obj$data)[1],")")

		#SETUP AXES--------------------------------
		if(!add){ #use default plotting options
			plot(X,Y, type='n', xlab=xlab, ylab=ylab, frame.plot = frame.plot,...)
   		}

		#add datapoints	--------------------------------
		if(type %in% c("o","b", "p")){
			if(obj$gt == "none")
				points(X, Y, col = col[1], pch=pch[1],...)
			else{
				for(i in 1:ngrps){
					iref <- as.character(obj$data[,3]) == groups[i]
					points(X[iref], Y[iref], col =col[i], pch=pch[i],...)
				}
			}
		}
	}

	# QQ PLOT
	if(whichplot == "qq")
	{
		 #obtain data--------------------------------
		Y <- fitted.sma(obj, type = "residuals")

		#axis labels--------------------------------
	    #determine axis labels if not already specified
	    if(is.null(xlab)) xlab <- "Normal quantiles"
	    if(is.null(ylab)) ylab <- paste("Residual values (",names(obj$data)[2], " v ",
			names(obj$data)[1],")")

		#SETUP AXES--------------------------------
		if(add==FALSE){ #use default plotting options
			qqnorm(Y, xlab=xlab, ylab=ylab,...)
			qqline(Y)
   		}
	}

}

