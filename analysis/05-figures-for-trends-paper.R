# plots for trends paper

library(tidyverse)
library(cowplot)

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

s1trends <- ggplot(s1quant, aes(year, median)) + 
  geom_line(alpha = 0.5) + geom_point(size = 1, alpha = 0.7) + 
  ylab(expression(paste(S[1]~"estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) +
  facet_wrap(~river_name, ncol = 2, scales = "free_y") + theme_bw() +
  xlab("Year") +
  labs(color = "River") +
  theme(legend.position = "none", 
        legend.text =  element_text(size = 12))
s1trends
ggsave("figures/s1-trends-faceted.png", height = 6, width = 5)

s1trends + facet_wrap(~river_name, ncol = 3, scales = "free_y") 
ggsave("figures/s1-trends-presentation.png", height = 4, width = 8)




morat <- tibble(year = c(1984, 1998, 2000),
                region = c("NS+\nNB", "NL", "QC"),
                river_name = c("Nashwaak River", "Western Arm Brook", "Trinité River"))

s1trends <- ggplot(s1quant, aes(year, median, color = river_name)) + 
  geom_line(alpha = 0.5) + geom_point(size = 2, alpha = 0.8) +
  ylab(expression(paste(S[1]~"posterior estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95), alpha = 0.5, width = 0.1) +
  theme_cowplot() + xlab("Year") + labs(color = "River") +
  theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 3)) +
  geom_vline(data = morat, aes(xintercept = year, color = river_name), alpha = 0.4, linetype = 2) +
  geom_text(data = morat, aes(year-c(1.5, 1.5, -1.5), 0.23, label = region, color = river_name), show.legend = FALSE)
s1trends
ggsave("figures/s1-trends.png", height = 5, width = 7)

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
  ylab(expression(paste(P[r]~"estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) +
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

