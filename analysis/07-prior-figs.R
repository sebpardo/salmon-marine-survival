library(tidyverse)
library(cowplot)


# S1 and S2 priors
x <- seq(0.01,2,by=0.01)
Z1 <- dlnorm(x, 1, 0.22)
exp(-Z1)
plot(Z1~exp(-Z1))

Z2 ~ lognormal(0.2, 0.3);

S1 <- seq(0.000,0.7, by = 0.001)

Z1 <- dlnorm(-log(S1), 1, 0.22)
Z2 <- dlnorm(-log(S1), 0.2, 0.3)
plot(Z1~S1, type = "l", lwd = 2)
lines(Z2/2~S1, lwd =2, col = "blue")

s1df <- data.frame(Survival = S1,
           Probability = Z1,
           Parameter = "S1")
s2df <- data.frame(Survival = S1,
                   Probability = Z2/2,
                   Parameter = "S2")

df1 <- bind_rows(s1df, s2df) %>% as_tibble() 

# Pr priors

x <- seq(0,1, by = 0.001)

# Prior for Conne
osw <- dlogis(qlogis(x), 2.4, 0.35)
plot(osw/2~x, type = "l", col ="blue")

# Prior for others
tsw <- dlogis(qlogis(x), 0, 0.8)
lines(tsw~x, col = "green")

pr1 <- data.frame(Proportion = x,
                  Probability = osw/2,
                  Population = "1SW rivers")

pr2 <- data.frame(Proportion = x,
                  Probability = tsw,
                  Population = "Other rivers")

df2 <- bind_rows(pr1, pr2) %>% as_tibble() 


prpriors <- ggplot(df2, aes(Proportion, Probability, linetype = Population)) + geom_line(size = 1.5) + 
  theme_cowplot() + xlab("Proportion returning after first winter at sea") 
prpriors
ggsave("~/Dropbox/2020-ASCF-webinar/pr-priors.png", height = 4, width = 6)

s1priors <- ggplot(df1, aes(Survival, Probability, color = Parameter)) + geom_line(size = 1.5) + 
  scale_colour_discrete(labels = expression(italic(S[1]), italic(S[2]))) +
  theme_cowplot()
s1priors
ggsave("~/Dropbox/2020-ASCF-webinar/surv-priors.png", height = 4, width = 6)

plot_grid(s1priors, prpriors,  labels = c("a)", "b)"), nrow = 2, label_size = 16, label_x = c(0.125, 0.125))
ggsave("figures/priors.png", height = 6.5, width = 5)
