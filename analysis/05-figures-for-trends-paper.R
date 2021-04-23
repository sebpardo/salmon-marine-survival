# plots for trends paper

library(tidyverse)
library(cowplot)
library(xtable)

s1post <- readRDS(file = "data/s1-posteriors.rds")
s1quant <- readRDS(file = "data/s1quant.rds")

# Add rows with NA values as to stop geom_line to connect years
# with data gaps
addnas <- tibble(pos = NA,
                 river_name = c("Saint-Jean River", "Trinité River", rep("LaHave River", 2)),
                 year = c(1998, 2007, 2012:2013),
                 median = c(NA),
                 q05 = c(NA),
                 q25 = c(NA),
                 q75 = c(NA),
                 q95 = c(NA),
                 zscore = NA
)

s1quant <- bind_rows(s1quant, addnas)

morat2 <- tibble(year = c(1984, 1984, 1992, 1992, 1992, 2000, 2000),
                 region = c("NB","NS", "NL", "NL", "NL", "QC", "QC"),
                 river_name = c("Nashwaak River", "LaHave River", 
                                "Western Arm Brook", "Campbellton River", "Conne River", 
                                "Trinité River", "Saint-Jean River"))

s1trends <- ggplot(s1quant, aes(year, median)) + 
  geom_line(alpha = 0.5) + geom_point(size = 1, alpha = 0.7) + 
  ylab(expression(paste(S[1]~"estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) +
  geom_vline(data = morat2, aes(xintercept = year), alpha = 0.4, linetype = 2) +
  facet_wrap(~river_name, ncol = 2, scales = "free_y") + theme_bw() +
  xlab("Year") +
  labs(color = "River") +
  theme(legend.position = "none", 
        legend.text =  element_text(size = 12))
s1trends

ggsave("figures/s1-trends-faceted.png", height = 6, width = 5)

s1trends + facet_wrap(~river_name, ncol = 3, scales = "free_y") 
ggsave("figures/s1-trends-presentation.png", height = 4, width = 8)



morat <- tibble(year = c(1984, 1992, 2000),
                region = c("NS+\nNB", "NL", "QC"),
                river_name = c("Nashwaak River", "Western Arm Brook", "Trinité River"))



s1trendsall <- ggplot(s1quant, aes(year, median, color = river_name)) + 
  geom_line(alpha = 0.5) + geom_point(size = 2, alpha = 0.8) +
  ylab(expression(paste(S[1]~"posterior estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95), alpha = 0.5, width = 0.1) +
  theme_cowplot() + xlab("Year") + labs(color = "River") +
  theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 3)) +
  geom_vline(data = morat, aes(xintercept = year, color = river_name), alpha = 0.4, linetype = 2) +
  geom_text(data = morat, aes(year-c(1.5, 1.5, -1.5), 0.23, label = region, color = river_name), show.legend = FALSE)
s1trendsall
ggsave("figures/s1-trends.png", height = 5, width = 7)


s1trendsall <- ggplot(s1quant, aes(year, median, color = river_name)) + 
  geom_line(alpha = 0.5) +
  geom_point(size = 1, alpha = 0.8) +
  ylab(expression(paste("Median"~S[1]~"posterior estimates", sep = ""))) +
#  geom_errorbar(aes(ymin = q05, ymax = q95), alpha = 0.5, width = 0.1) +
  theme_cowplot() + xlab("Year") + labs(color = "River") +
  theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 3)) +
 # geom_vline(data = morat, aes(xintercept = year, color = river_name), alpha = 0.4, linetype = 2) +
#geom_text(data = morat, aes(year-c(1.5, 1.5, -1.5), 0.23, label = region, color = river_name), show.legend = FALSE)
  NULL
s1trendsall
ggsave("figures/s1-trends-all.png", height = 5, width = 7)



s1zscore <- ggplot(s1quant, aes(year, zscore, color = river_name)) + 
  geom_line(alpha = 1) + #geom_point(size = 2) + 
  ylab(expression("Z-scores of median"~paste(S[1]~"estimates", sep = ""))) +
  theme_cowplot() + xlab("Year") + labs(color = "River") +
  #scale_color_manual(values = rep("black", 7)) +
  theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 3))
#  theme(legend.position = "none")
s1zscore
ggsave("figures/s1-zscore-trends.png", height = 5, width = 7)

plot_grid(s1trends, s1zscore, labels = c("a)", "b)"), label_size = 16, label_x = c(0.125, 0.1))

ggsave("figures/s1-trends-dual.png", height = 5, width = 13)


s2quant <- readRDS(file = "data/s2quant.rds") %>%
  bind_rows(addnas)

prquant <- readRDS(file = "data/prquant.rds") %>%
  bind_rows(addnas)

s2trends <- ggplot(s2quant, aes(year, median)) + 
  geom_line(alpha = 0.5) + geom_point(size = 1, alpha = 0.7) + 
  ylab(expression(paste(S[2]~"estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) +
  geom_vline(data = morat2, aes(xintercept = year), alpha = 0.4, linetype = 2) +
  facet_wrap(~river_name, ncol = 2) + theme_bw() +
  xlab("Year") +
  labs(color = "River") +
  theme(legend.position = "none", 
        legend.text =  element_text(size = 12))
s2trends
ggsave("figures/s2-trends-faceted.png", height = 6, width = 5)
s2trends + facet_wrap(~river_name, ncol = 3, scales = "free_y") 
ggsave("figures/s2-trends-presentation.png", height = 4, width = 8)


prtrends <- ggplot(prquant, aes(year, median)) + 
  geom_line(alpha = 0.5) + geom_point(size = 1, alpha = 0.7) + 
  ylab(expression(paste(P[g]~"estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) +
  geom_vline(data = morat2, aes(xintercept = year), alpha = 0.4, linetype = 2) +
  facet_wrap(~river_name, ncol = 2) + theme_bw() +
  xlab("Year") + ylim(c(0,1))+ 
  labs(color = "River") +
  theme(legend.position = "none", 
        legend.text =  element_text(size = 12))
prtrends
ggsave("figures/pr-trends-faceted.png", height = 6, width = 5)
prtrends + facet_wrap(~river_name, ncol = 3, scales = "free_y") 
ggsave("figures/pr-trends-presentation.png", height = 4, width = 8)


# Graph of population level Pr estimates

prmuquant <- readRDS("data/prmuquant.rds")

prmuquant$param <- factor(prmuquant$param, levels = c("Pr_mu", "logitPr_CV", "logitPr_mu", "logitPr_sigma"))
prmuquant$param2 <- factor(prmuquant$param,  levels = c("Pr_mu", "logitPr_CV", "logitPr_mu", "logitPr_sigma"),
                           labels = c( "italic(Pr)[mu]", "logit(italic(Pr))[CV]","logit(italic(Pr))[mu]", "logit(italic(Pr))[sigma]"))
#labels = c( "italic(Pr)[mu]", "logit(italic(Pr)[CV])","logit(italic(Pr)[mu])", "logit(italic(Pr)[sigma])"))


prmuplot <- prmuquant %>%
  ungroup() %>%
  mutate(river_name = fct_rev(river_name)) %>%
 # filter(param == "Pr_mu") %>%
ggplot(aes(median, river_name)) + 
 geom_point(size = 2, alpha = 1) + 
  geom_errorbarh(aes(xmin = q25, xmax = q75), size = 1, height = 0.00) +
  geom_errorbarh(aes(xmin = q05, xmax = q95), size = 0.5, height = 0.00) +
  ylab("River") + theme_bw() +
 xlab("Posterior estimates") +
  # xlab(expression(paste("Posterior"~Pr[mu]~"estimates", sep = ""))) +
  facet_wrap(~param2, ncol = 4, scales = "free_x", shrink = FALSE, labeller = label_parsed) +
  theme(strip.text.x = element_text(size = 12), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
prmuplot
ggsave("figures/pr-mu-posteriors.png", height = 3, width = 10)
  

s1quant %>%
  group_by(river_name) %>%
  summarise(min = min(median), max = max(median))

s1quant %>%
  filter(median <= 0.02)

s1quant %>%
  filter(median >= 0.17)


# Plotting individual population level posteriors of pr based on logitpr mu and logitpr sigma

fithierncp <- readRDS("data/stanfit_hier_ncp.rds") # hierarchical
propgrilse <- readRDS("data/propgrilse.rds")
alldatastan <- readRDS("data/alldatastan.rds")

logitpr_mu <- as.matrix(fithierncp, pars = "logitPr_mu") %>% 
  t %>% as_tibble %>% mutate(river_name = alldatastan$river_name) %>%
  pivot_longer(V1:V7500, names_to = "iteration") %>%
  select(river_name, everything()) %>%
  group_by(river_name) %>%
  rename(`Population` = river_name, logitPr_mu = value)

logitpr_mu

logitpr_mu_quant <- logitpr_mu %>%
  rename(value = logitPr_mu) %>%
  ungroup() %>%
  mutate(Population = fct_relevel(Population, propgrilse$river_name)) %>%
  group_by(Population) %>%
  summarise(Median = median(value), 
            `2.5%` = quantile(value, probs = 0.025), 
            `25%` = quantile(value, probs = 0.25), 
            `75%` = quantile(value, probs = 0.75), 
            `97.5%` = quantile(value, probs = 0.975)
            # n = n()
  )


logitpr_sigma <- as.matrix(fithierncp, pars = "logitPr_sigma") %>% 
  t %>% as_tibble %>% mutate(river_name = alldatastan$river_name) %>%
  pivot_longer(V1:V7500, names_to = "iteration") %>%
  select(river_name, everything()) %>%
  group_by(river_name) %>%
  rename(`Population` = river_name, logitPr_sigma = value)

logitpr_sigma

logitpr_sigma_quant <- logitpr_sigma %>%
  rename(value = logitPr_sigma) %>%
  ungroup() %>%
  mutate(Population = fct_relevel(Population, propgrilse$river_name)) %>%
  group_by(Population) %>%
  summarise(Median = median(value), 
            `2.5%` = quantile(value, probs = 0.025), 
            `25%` = quantile(value, probs = 0.25), 
            `75%` = quantile(value, probs = 0.75), 
            `97.5%` = quantile(value, probs = 0.975)
            # n = n()
  )

# Tables of logitPr mu and sigma

left_join(logitpr_mu_quant, logitpr_sigma_quant, by = "Population")

print.xtable(xtable(logitpr_mu_quant, digits = c(0,0,3,3,3,3,3),
                    caption = "Posterior estimates of mean population-level
                     proportion returning after one winter at sea 
                    ($\\mu_{Pg,r}$, logit-transformed) for the seven populations examined.",
                    label = "tab:prmu",
                    NULL),
             file = "ms/table-logitprmu.tex",
             caption.placement = "top",
             include.rownames = FALSE)

 print.xtable(xtable(logitpr_sigma_quant, digits = c(0,0,3,3,3,3,3),
                    caption = "Posterior estimates of population-level
                    standard deviations of the proportion returning 
                    after one winter at sea 
                    ($\\sigma_{Pg,r}$, logit-transformed) for the seven populations examined.",
                    label = "tab:prsigma"),
             file = "ms/table-logitprsigma.tex",
             caption.placement = "top",
             include.rownames = FALSE)


median_logitpr <- tibble(Population = logitpr_mu_quant$Population, 
                             mu = logitpr_mu_quant$Median,
                             sigma = logitpr_sigma_quant$Median)

logitpr <- logitpr_mu %>% ungroup %>% mutate(logitPr_sigma = logitpr_sigma$logitPr_sigma)

invlogit <- function(x, mu, sig) dlogis(qlogis(x), mu, sig)
# invlogit(0.9, 2.5, 0.4)
x <- seq(0,1, by = 0.001)

sublogitpr <- logitpr %>%
  group_by(Population) %>% 
  sample_n(size = 50) %>%
  rename(mu = logitPr_mu, sigma = logitPr_sigma) %>% 
  group_by(Population, iteration, mu, sigma) %>%
  expand(., x = x) %>%
  ungroup() %>%
  mutate(y = invlogit(x, mu = mu, sig=sigma),
         Population = fct_relevel(Population, propgrilse$river_name),
         iter2 = paste(Population, iteration)) %>% group_by(Population)

p1 <- ggplot(sublogitpr, aes(x, y, color = Population)) + geom_line(aes(group = iter2), alpha = 0.3)
p1

pop_logitpr <- median_logitpr %>%
  group_by(Population, mu, sigma) %>%
  expand(., x = x) %>%
  mutate(y = invlogit(x, mu = mu, sig=sigma))

ggplot(pop_logitpr, aes(x, y, color = Population)) + geom_line(aes(group = Population), size = 2)

p1 + geom_line(inherit.aes = FALSE, data = pop_logitpr, aes(x, y, group = Population), 
               size = 1.3, color = "gray20", alpha = 0.5) + 
  coord_cartesian(ylim = c(0, 5)) + guides(color = guide_legend(override.aes = list(size = 1.3, alpha = 1))) +
  theme_cowplot() +
  xlab("Proportion returning after first winter at sea") + ylab("Density")
ggsave("figures/logitpr-mu-sigma-post.png", width = 7, height = 4)
