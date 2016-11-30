### Functions to generate theo tradeoffs simulation


generate_init <- function(nsp = 100, NN = 4*256, Nlandscape = 10){
 df_sp <- GenerateRandSp(nsp)
 list_init <- InitLandscape(df_sp, NN = NN, Nlandscape = Nlandscape,
                            min_ss = 0, max_ss = 1, per_occup = 0.5)
 return(list_init)
}

run_K <- function(list_init, K){
require(TheoTradeOff)
 UpdateIterR(list_init[['mat_sp']], list_init[['mat_suc']],list_init[['c_e']], list_init[['c_l']], list_init[['c_s']],list_init$ss,prob_distur = 0.05, prob_suc = 0.2, K , n = 8000)
}

plot_map <- function(l_res, list_init, val_K){
    image_landscape(l_res,
                     list_init$c_e,
                     list_init$c_l,
                     list_init$c_s)
     mtext(side = 1, text = paste0('Abiotic stress gradient K =', val_K), line = -1.2, cex = 1.3)
     arrows(250, -200, 10000, -200)
}

plot_mean_param <- function(l_res, list_init, sp_param = 'c_e',
                            y_lab = "Mean early successional competitive ability"){
mean_c_l <- mean_landscape_l(l_res, list_init[[sp_param]])
x <- seq_len(length(mean_c_l))
lo <- loess(mean_c_l~x, span = 0.4)
plot(x, mean_c_l,  pch = 16, cex = 0.5,
     xlab = "Abiotic stress gradient",
     ylab = y_lab,
     cex.lab = 1.6)
lines(predict(lo), col='red', lwd=2)
}

generate_data_sp_ranges <- function(l_res){
require(reldist)
l_res$sp[l_res$sp == 0] <- NA
sp_mat <- l_res$sp
sp <- as.numeric(names(table(l_res$sp)))
data_ranges <- data.frame(sp = sp,
                          q_l = rep(NA, length(sp)),
                          q_h = rep(NA, length(sp)),
                          q_m = rep(NA, length(sp)),
                          q_r = rep(NA, length(sp)))

for (i in sp){
sp_obs <- apply(sp_mat == i & !is.na(sp_mat), MARGIN = 2, FUN = sum)
sp_obs[is.na(sp_obs)] <-  0
vec_pos <- seq_len(max(dim(l_res$sp)))
q_l <- wtd.quantile(vec_pos, q = 0.1, weight = sp_obs)
q_h <- wtd.quantile(vec_pos, q = 0.9, weight = sp_obs)
q_m <- wtd.quantile(vec_pos, q = 0.5, weight = sp_obs)
data_ranges[data_ranges$sp == i, "q_l"] <- wtd.quantile(vec_pos, q = 0.001, weight = sp_obs)
data_ranges[data_ranges$sp == i, "q_h"] <- wtd.quantile(vec_pos, q = 0.999, weight = sp_obs)
data_ranges[data_ranges$sp == i, "q_m"] <- wtd.quantile(vec_pos, q = 0.5, weight = sp_obs)
data_ranges[data_ranges$sp == i, "q_r"] <-
    data_ranges[data_ranges$sp == i, "q_h"]  - data_ranges[data_ranges$sp == i, "q_l"]
}

return(data_ranges)
}

plot_range <- function(data_ranges){
plot(data_ranges$q_m, data_ranges$q_r, xlab = "Median of species abiotic stress distribution",
     ylab = "Range in species abiotic stress distribution",
     cex.lab = 1.6)
}
