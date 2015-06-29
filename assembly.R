
source("R/assembly.R")

p <- trait_gradients_base_parameters()

p$disturbance_mean_interval <- 10.0
sys0 <- community(p, bounds_infinite("lma"))
#sys0 <- community(p, bounds(lma= c(-Inf, Inf), stc=c(0, 100)))

obj_m0 <- assembler(sys0, list(birth_move_tol=1))
obj_m <- assembler_run(obj_m0, 20)
obj_m$done

ff <- lapply(obj_m$history, community_fitness_approximate)
cols <- c("black", "blue", "orange")
lma <- seq_log_range(obj_m0$community$bounds, 400)
w <- sapply(ff, function(f) f(lma))
res_lma <- sapply(obj_m$history[-1], function(x) x$traits)
res_w0 <- sapply(seq_along(res_lma), function(i) ff[[i  ]](res_lma[[i]]))
res_w1 <- sapply(seq_along(res_lma), function(i) ff[[i+1]](res_lma[[i]]))
matplot(lma, w, type="l", col=cols, lty=1, ylim=c(-1, max(w)),
        log="x")
abline(h=0, col="grey")
segments(res_lma, res_w0, res_lma, res_w1, col=cols[-1])
points(res_lma, res_w1, col=cols[-1], pch=19)
