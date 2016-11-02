
data_intro_fig <- function(){
##
# Figure to plot species distribution according to a Lotka-Volterra model based a gaussian fundamental niche with symmetric competition or a shifting hierarchy Keddy's like with monotonic increasing fundamental niche and asymmetric competition

## Figures of abstract represtentation of shifting competitive hierarchy vs. fundamental niche differentiation

#shifting competitive hierarchy
Ks_fun<- function(p, a) a*p/(a^4+p)
#gaussian niche
Kb_fun<- function(p, a) dnorm(p, mean = a, sd = 0.7)

# competition matrix
a_seq_s<- c(1, 2/3, 1/3+0.1)
a_seq_b<- rev(c(0.75, 1.5, 2.25))

mat_alpha_s<- matrix(c(1,   1/4, 1/8,
                       1/2, 1,   1/4,
                       1,   1/2, 1),
                     nrow = 3, ncol = 3,
                     byrow = TRUE)
mat_alpha_b<- matrix(c(1,   1/4, 1/4,
                       1/4, 1,   1/4,
                       1/4, 1/4, 1),
                     nrow = 3, ncol = 3,
                     byrow = TRUE)

### OD for LV model and solving

library(deSolve)

LVmod3 <- function(Time, State, Pars) {
with(append(State, Pars), {
  dA<- A*(1 - M[1,1]*A/K[1] - M[1,2]*B/K[1] - M[1,3]*C/K[1])
  dB<- B*(1 - M[2,1]*A/K[2] - M[2,2]*B/K[2] - M[2,3]*C/K[2])
  dC<- C*(1 - M[3,1]*A/K[3] - M[3,2]*B/K[3] - M[3,3]*C/K[3])
  return(list(c(dA, dB, dC)))
 })
}


time <- c(0,100)
N_seq <- 1000
p <- seq(0.0001, 3, length.out = N_seq)

res_s<- matrix(NA, nrow = N_seq, ncol = 3)
res_b<- matrix(NA, nrow = N_seq, ncol = 3)

for (j in 1:N_seq){
# initial conditions: also a named vector
state <- c(A = 0.33, B = 0.33, C = 0.33)
#shift competitive hierarchy
parameters <- list(K = Ks_fun(p[j], a_seq_s),
                   M = sqrt(mat_alpha_s))
out <- ode(y = state, times = time, func = LVmod3, parms = parameters)
res_s[j, ] <- out[2, 2:4]
#gaussian
parameters <- list(K = Kb_fun(p[j], a_seq_b),
                   M = sqrt(mat_alpha_b))
out <- ode(y = state, times = time, func = LVmod3, parms = parameters)
res_b[j, ] <- out[2, 2:4]
}

## plot basic results

Ks_p <- t(sapply(p, Ks_fun, a = a_seq_s))
Kb_p <- cbind(Kb_fun(p, a = a_seq_b[1]),
              Kb_fun(p, a = a_seq_b[2]),
              Kb_fun(p, a = a_seq_b[3]))

df_pred_A<- data.frame(gradient = c(p, p, p, p),
                      Abund = c(Kb_p[, 1], Ks_p[, 1], res_b[, 1], res_s[ , 1]),
                      Model = c(rep('Gaussian niche', length(p)),
                                rep('Shifting competitive hierarchy',
                                    length(p)),
                                rep('Gaussian niche', length(p)),
                                rep('Shifting competitive hierarchy',
                                    length(p))),
                      Pred = c(rep('Fundamental niche', length(p)),
                               rep('Fundamental niche', length(p)),
                               rep('Realised niche', length(p)),
                               rep('Realised niche', length(p))),
                      Species = rep('A', length(p)*4))
df_pred_B<- data.frame(gradient = c(p, p, p, p),
                      Abund = c(Kb_p[, 2], Ks_p[, 2], res_b[, 2], res_s[ , 2]),
                      Model = c(rep('Gaussian niche', length(p)),
                                rep('Shifting competitive hierarchy',
                                    length(p)),
                                rep('Gaussian niche', length(p)),
                                rep('Shifting competitive hierarchy',
                                    length(p))),
                      Pred = c(rep('Fundamental niche', length(p)),
                               rep('Fundamental niche', length(p)),
                               rep('Realised niche', length(p)),
                               rep('Realised niche', length(p))),
                      Species = rep('B', length(p)*4))
df_pred_C<- data.frame(gradient = c(p, p, p, p),
                      Abund = c(Kb_p[, 3], Ks_p[, 3], res_b[, 3], res_s[ , 3]),
                      Model = c(rep('Gaussian niche', length(p)),
                                rep('Shifting competitive hierarchy',
                                    length(p)),
                                rep('Gaussian niche', length(p)),
                                rep('Shifting competitive hierarchy',
                                    length(p))),
                      Pred = c(rep('Fundamental niche', length(p)),
                               rep('Fundamental niche', length(p)),
                               rep('Realised niche', length(p)),
                               rep('Realised niche', length(p))),
                      Species = rep('C', length(p)*4))
df_pred <- rbind(df_pred_A, df_pred_B, df_pred_C)
return(df_pred)
}

plot_intro_fig_0a<- function(df_pred){
library(ggplot2)
df_pred$gradient <- 3.0001 - df_pred$gradient
df_pred$Abund[df_pred$Abund<0] <- 0
df <- df_pred[df_pred$Model != 'Shifting competitive hierarchy', ]

ggplot(df, aes(x = gradient, y = Abund, colour = Species)) +
  geom_line(aes(x = gradient, y = Abund, colour = Species), size = 1.5) +
  facet_grid(Pred ~ .) +  theme_simple() +
  ylab("Abundance")+xlab("Abiotic stess")+ylim(0, 0.75)
ggsave(file.path('figures', 'resLV_0a.pdf'))

}

plot_intro_fig_0b<- function(df_pred){
library(ggplot2)
df_pred$gradient <- 3.0001 - df_pred$gradient
df_pred$Abund[df_pred$Abund<0] <- 0

df <- df_pred[df_pred$Model == 'Shifting competitive hierarchy', ]

ggplot(df, aes(x = gradient, y = Abund, colour = Species)) +
  geom_line(aes(x = gradient, y = Abund, colour = Species), size = 1.5) +
  facet_grid(Pred ~ .) +  theme_simple() +
  ylab("Abundance")+xlab("Abiotic stress")+ylim(0, 0.75)
ggsave(file.path('figures', 'resLV_0b.pdf'))

}


plot_intro_fig_1<- function(df_pred){
library(ggplot2)
df_pred$Abund[df_pred$Abund<0] <- 0

df_pred$Abund[df_pred$Model == 'Shifting competitive hierarchy'] <- NA

ggplot(df_pred, aes(x = gradient, y = Abund, colour = Species)) +
  geom_line(aes(x = gradient, y = Abund, colour = Species), size = 1.5) +
  facet_grid(Pred ~ Model) +  theme_simple() +
  ylab("Abundance")+xlab("Abiotic gradient")+ylim(0, 0.75)
ggsave(file.path('figures', 'resLV_1.pdf'))

}


plot_intro_fig_2<- function(df_pred){
library(ggplot2)
df_pred$Abund[df_pred$Abund<0] <- 0

ggplot(df_pred, aes(x = gradient, y = Abund, colour = Species)) +
  geom_line(aes(x = gradient, y = Abund, colour = Species), size = 1.5) +
  facet_grid(Pred ~ Model) +  theme_simple() +
  ylab("Abundance")+xlab("Abiotic gradient")+ylim(0, 0.75)
ggsave(file.path('figures', 'resLV_2.pdf'), width = 12, height = 7)
}


theme_simple <-  function(){
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    strip.background = element_blank(),
    legend.key = element_blank(),
    legend.position=c(0.045,.4),
    strip.text = element_text(size=14),
   axis.title.x = element_text(size=14),
   axis.title.y = element_text(size=14)
  ) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))
}
