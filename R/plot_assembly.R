
format_list_res_vpd_cluster <- function(name,
                                        seq_vpd = seq(0, 30, by = 5)/10,
                                        folder = "output_cluster"){
  f <- function(vpd, folder, name)  readRDS(file.path(folder,
                                                      paste0(name,vpd,".rds")))

  list_res <- lapply(seq_vpd, f, folder, name)
  return(list_res)
}



plot_trait_gradient <- function(list_assembly_lma, vec_site_prod, var = "lma",
                                varlab = expression(paste("Leaf-mass per area (kg ", m^-2,")"))) {
  list_data <- list_assembly_lma
  grad <- vec_site_prod
  y <- lapply(list_data, function(x) x$community$traits[ ,var])

  plot(NA, type="n", log="y", las=1, xlim=range(grad),
       ylim= range(unlist(y)),
    ylab= varlab,
    xlab="Site productivity")


  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], y[[i]], type='p', col="black", pch=16)
  }
}

plot_trait_vpd <- function(list_assembly_lma, vec_site_vpd, var = "lma",
                           varlab = expression(paste("Leaf-mass per area (kg ", m^-2,")"))) {
  list_data <- list_assembly_lma
  grad <- vec_site_vpd
  y <- lapply(list_data, function(x) x$community$traits[ ,var])

  plot(NA, type="n", log="y", las=1, xlim = range(grad),
       ylim= range(unlist(y)),
    ylab= varlab,
    xlab="VPD")


  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], y[[i]], type='p', col="black", pch=16)
  }
}


plot_trait_gradient_narea_lma<- function(list_assembly_lma,
                                         vec_site_prod,
                                         var = "narea",
               varlab = expression(paste("Nitrogen per area (kg ", m^-2,")"))) {
  list_data <- list_assembly_lma
  grad <- vec_site_prod
  y <- lapply(list_data, function(x) x$community$traits[ ,var])
  z <- lapply(list_data, function(x) x$community$traits[ ,"lma"])


  plot(NA, type="n", log="y", las=1, xlim=range(grad),
       ylim= range(unlist(y)),
    ylab= varlab,
    xlab="Site productivity")

  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], y[[i]], type='p', col="black",
           pch=1, cex = 0.5+10*z[[i]])
  }
}

plot_trait_vpd_narea_lma<- function(list_assembly_lma,
                                         vec_site_prod,
                                         var = "narea",
               varlab = expression(paste("Nitrogen per area (kg ", m^-2,")"))) {
  list_data <- list_assembly_lma
  grad <- vec_site_prod
  y <- lapply(list_data, function(x) x$community$traits[ ,var])
  z <- lapply(list_data, function(x) x$community$traits[ ,"lma"])


  plot(NA, type="n", log="y", las=1, xlim=range(grad),
       ylim= range(unlist(y)),
    ylab= varlab,
    xlab="vpd")

  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], y[[i]], type='p', col="black",
           pch=1, cex = 0.5+10*z[[i]])
  }
}


plot_cor_narea_lma<- function(list_assembly_lma,
                                         vec_site_prod,
                                         var = "narea") {
  list_data <- list_assembly_lma
  y <- lapply(list_data, function(x) x$community$traits[ ,var])
  z <- lapply(list_data, function(x) x$community$traits[ ,"lma"])

  plot(unlist(y), unlist(z), type="n", log="y", las=1,
    xlab = expression(paste("Nitrogen per area (kg ", m^-2,")")),
    ylab = expression(paste("Leaf-mass per area (kg ", m^-2,")")))
colfunc <- colorRampPalette(c("red", "green"))
cols <- colfunc(length(list_data))
  for(i in seq_along(list_data)) {
    points(y[[i]] , z[[i]], type='p', col=cols[i],
           pch=16)
  }
}


plot_trait_gradient_narea_lma2<- function(list_assembly_lma,
                                         vec_site_prod,
                                         var = "narea",
               varlab = expression(paste("Nitrogen per area (kg ", m^-2,")"))) {
  list_data <- list_assembly_lma
  grad <- vec_site_prod
  y <- lapply(list_data, function(x) x$community$traits[ ,var])
  z <- lapply(list_data, function(x) x$community$traits[ ,"lma"])

  par(mfrow = c(2,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,4,1,1) + 0.1)
  plot(NA, type="n", log="y", las=1, xlim=range(grad),
       ylim= range(unlist(y)),  axes = FALSE, ylab = varlab)
  axis(side = 1, labels = FALSE, at = grad)
  axis(side = 2)
  box()

  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], y[[i]], type='p', col="black",
           pch=16)
  }
  plot(NA, type="n", log="y", las=1, xlim=range(grad),
       ylim= range(unlist(z)),ylab= expression(paste("Leaf-mass per area (kg ", m^-2,")")),
    xlab=NA)

  for(i in seq_along(list_data)) {
    points(y[[i]]*0 + grad[i], z[[i]], type='p', col="black",
           pch=16)
  }

title(xlab = "Site productivity",
      outer = TRUE, line = 3)
}


plot_fitness_assembly <- function(list_assembly_lma, vec_site_prod) {
  list_data <- list_assembly_lma
  grad <- vec_site_prod


  cols <- rev(brewer.pal(5, "Blues"))

  plot(NA, log="x", xlim=c(0.01, 1), ylim=c(-2,0.5),
    xlab= expression(paste("Leaf-mass per area (kg ", m^-2,")")),
    ylab="Fitness")
  abline(h=0, col="grey")

  for(i in seq_along(list_data)) {
    data <- list_data[[i]]
    community <- last(data$history)
    ff <- community_fitness_approximate(community)
    lma <- sort(c(seq_log_range(data$community$bounds, 400), data$community$traits))
    w <- ff(lma)
    lma_res <- data$community$traits[,"lma"]
    points(lma,w, type='l', col=cols[i])
    points(lma_res,ff(lma_res), type='p', col=cols[i], pch=19)
  }
  legend("topleft", bty="n", legend=c("low", "med", "high"), col=cols, pch=19)
}

