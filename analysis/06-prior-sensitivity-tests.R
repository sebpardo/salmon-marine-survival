# Sensitivity analyses of models with differing priors
library(tidyverse)
library(rstan)

alldatastan <- readRDS("data/alldatastan.rds")

# Loading stanfit objects 
fithierncp <- readRDS("data/stanfit_hier_ncp_acv_full.rds")
fithierncp_prweak <- readRDS("data/stanfit_hier_ncp_acv_prweak.rds")
fithierncp_zsdmu <- readRDS("data/stanfit_hier_ncp_acv_zsdmu.rds")


# Indexing positions of years for each river
indexdf <- tibble(pos = as.character(1:length(alldatastan$years)),
                  river_name = rep(alldatastan$river_name, 
                                   times = alldatastan$nyears), 
                  year = alldatastan$years)

# For indexing population level posteriors
riverindexdf <- tibble(pos = as.character(1:length(alldatastan$river_name)),
                       river_name = alldatastan$river_name)

# Extracting posterior distributions
extract_post <- function(stanfit, pars = "S1", index = indexdf) {
  as.matrix(stanfit, pars = pars) %>%
    as_tibble %>%
    gather(parameter, value) %>%
    mutate(pos = str_extract(parameter, "(?<=\\[).+?(?=\\])"),
           param = str_extract(parameter, ".+?(?=\\[)")) %>%
    left_join(index, by = "pos")
}

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

morat2 <- tibble(year = c(1984, 1984, 1992, 1992, 1992, 2000, 2000),
                 region = c("NB","NS", "NL", "NL", "NL", "QC", "QC"),
                 river_name = c("Nashwaak River", "LaHave River", 
                                "Western Arm Brook", "Campbellton River", "Conne River", 
                                "Trinité River", "Saint-Jean River"))


# Calculate quantiles from extracted posteriors
calc_quant <- function(posterior) {
  posterior %>%
    group_by(pos, river_name, year) %>%
    summarise(median = median(value),
              q05 = quantile(value, 0.05),
              q25 = quantile(value, 0.25),
              q75 = quantile(value, 0.75),
              q95 = quantile(value, 0.95)) %>%
    mutate(pos = as.numeric(pos)) %>%
    arrange(pos)
}

s1post <- extract_post(fithierncp)
s1post_prweak <- extract_post(fithierncp_prweak)
s1post_zsdmu <- extract_post(fithierncp_zsdmu)

s1quant <- calc_quant(s1post) %>% bind_rows(addnas) %>%
  mutate(priors = "Strong~priors")
s1quant_prweak <- calc_quant(s1post_prweak) %>% bind_rows(addnas) %>%
  mutate(priors = "Weak~italic(P[g])~priors")
s1quant_zsdmu <- calc_quant(s1post_zsdmu) %>% bind_rows(addnas) %>%
  mutate(priors = "Weak~italic(S[1])~and~italic(S[2])~priors")

stpr <- "Strong~priors"
wkpg <- "Weak~italic(P[g])~priors"
wkz <- "Weak~italic(S[1])~and~italic(S[2])~priors"

s1all <- bind_rows(s1quant, s1quant_prweak, s1quant_zsdmu)

s1priors <- ggplot(s1all, aes(year, median)) + 
  geom_line(alpha = 0.5) + geom_point(size = 1, alpha = 0.7) + 
  ylab(expression(paste(S[1]~"estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) +
  geom_vline(data = morat2, aes(xintercept = year), alpha = 0.4, linetype = 2) +
  facet_grid(river_name ~ priors, scales = "free_y",
             labeller = labeller(river_name = label_wrap_gen(width = 10),
                                 priors = label_parsed)) + 
  theme_bw() +
  xlab("Year") +
  labs(color = "River") +
  theme(legend.position = "none", 
        legend.text =  element_text(size = 12))
s1priors
ggsave("figures/s1-trends-priors.png", height = 8, width = 7)

s2post <- extract_post(fithierncp, pars = "S2")
s2post_prweak <- extract_post(fithierncp_prweak, pars = "S2")
s2post_zsdmu <- extract_post(fithierncp_zsdmu, pars = "S2")

s2quant <- calc_quant(s2post) %>% bind_rows(addnas) %>%
  mutate(priors = stpr)
s2quant_prweak <- calc_quant(s2post_prweak) %>% bind_rows(addnas) %>%
  mutate(priors = wkpg)
s2quant_zsdmu <- calc_quant(s2post_zsdmu) %>% bind_rows(addnas) %>%
  mutate(priors = wkz)

s2all <- bind_rows(s2quant, s2quant_prweak, s2quant_zsdmu)

s2priors <- ggplot(s2all, aes(year, median)) + 
  geom_line(alpha = 0.5) + geom_point(size = 1, alpha = 0.7) + 
  ylab(expression(paste(S[2]~"estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) +
  geom_vline(data = morat2, aes(xintercept = year), alpha = 0.4, linetype = 2) +
  facet_grid(river_name ~ priors, scales = "free_y",
             labeller = labeller(river_name = label_wrap_gen(width = 10),
                                 priors = label_parsed)) + 
  theme_bw() +
  xlab("Year") +
  labs(color = "River") +
  theme(legend.position = "none", 
        legend.text =  element_text(size = 12))
s2priors
ggsave("figures/s2-trends-priors.png", height = 8, width = 7)

prpost <- extract_post(fithierncp, pars = "Pr")
prpost_prweak <- extract_post(fithierncp_prweak, pars = "Pr")
prpost_zsdmu <- extract_post(fithierncp_zsdmu, pars = "Pr")

prquant <- calc_quant(prpost) %>% bind_rows(addnas) %>%
  mutate(priors = stpr)
prquant_prweak <- calc_quant(prpost_prweak) %>% bind_rows(addnas) %>%
  mutate(priors = wkpg)
prquant_zsdmu <- calc_quant(prpost_zsdmu) %>% bind_rows(addnas) %>%
  mutate(priors = wkz)

prall <- bind_rows(prquant, prquant_prweak, prquant_zsdmu)

prpriors <- ggplot(prall, aes(year, median)) + 
  geom_line(alpha = 0.5) + geom_point(size = 1, alpha = 0.7) + 
  ylab(expression(paste(P[g]~"estimates", sep = ""))) +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) +
  geom_vline(data = morat2, aes(xintercept = year), alpha = 0.4, linetype = 2) +
  facet_grid(river_name ~ priors, 
             labeller = labeller(river_name = label_wrap_gen(width = 10),
                                 priors = label_parsed)) + 
  theme_bw() +
  xlab("Year") +
  labs(color = "River") +
  theme(legend.position = "none", 
        legend.text =  element_text(size = 12))
prpriors
ggsave("figures/pr-trends-priors.png", height = 8, width = 7)

#  Plotting S1 and S2 priors
x <- seq(0.01,2,by=0.01)
Z1 <- dlnorm(x, 1, 0.22)
exp(-Z1)
plot(Z1~exp(-Z1))

Z2 ~ lognormal(0.2, 0.3);

S1 <- seq(0.000,0.8, by = 0.001)

Z1 <- dlnorm(-log(S1), 1, 0.22)
Z2 <- dlnorm(-log(S1), 0.2, 0.3)

Z1weak <- dlnorm(-log(S1), 0.5, 0.5)
Z2weak <- dlnorm(-log(S1), 0.1, 0.5)

plot(Z1~S1, type = "l", lwd = 2)
lines(Z2/2~S1, lwd =2, col = "blue")

plot(Z1weak~S1, type = "l", lwd = 2)
lines(Z2weak/1.5~S1, lwd =2, col = "blue")

s1df <- data.frame(Survival = S1,
                   Probability = Z1,
                   Parameter = "S1",
                   Prior = "strong")
s1weakdf <- data.frame(Survival = S1,
                   Probability = Z1weak,
                   Parameter = "S1",
                   Prior = "weak")
s2df <- data.frame(Survival = S1,
                   Probability = Z2/2,
                   Parameter = "S2",
                   Prior = "strong")
s2weakdf <- data.frame(Survival = S1,
                   Probability = Z2weak/2,
                   Parameter = "S2",
                   Prior = "weak")

s1df <- bind_rows(s1df, s1weakdf)
s2df <- bind_rows(s2df, s2weakdf)


df1 <- bind_rows(s1df, s2df) %>% as_tibble() 
ggplot(df1, aes(Survival, Probability, color = Parameter, 
                            linetype = Prior)) + geom_line(size = 1.5) + 
  scale_linetype_manual(values=c(1,3)) +
 scale_colour_discrete(labels = expression(italic(S[1]), italic(S[2]))) +
  theme_cowplot() + theme(legend.key.size=unit(2,"lines"))

ggsave("figures/z-weak-prior.png", height = 3.5, width = 5.5)


# Plotting Pr prior
# Prior for weak pr (all popns)
x <- seq(0,1, by = 0.001)
pr_weak <- dlogis(qnorm(x), 0, 2.8) # probit
# tsw <- dlogis(qlogis(x), 0, 0.8) # logit
lines(pr_weak ~ x, col = "green")

prweak <- data.frame(Proportion = x,
                  Probability = pr_weak,
                  Population = "All",
                  Prior = "weak")


prpriors <- ggplot(prweak, aes(Proportion, Probability)) + 
  geom_line(size = 1.5) + 
  theme_cowplot() + xlab("Proportion returning after first winter at sea") 
prpriors

ggsave("figures/pr-weak-prior.png", height = 3.5, width = 5.5)

# Recaculating correlation of Z1 posteriors with weak Z1 and Z2 priors

estreturns <- tibble(river_name = rep(alldatastan$river_name, times = alldatastan$nyears),
                     year = alldatastan$years,
                     #logsmolts = alldatastan$logsmolts, 
                     log1SW = alldatastan$loggrilse,
                     log2SW = alldatastan$logSW2)

propgrilse <- estreturns %>%
  mutate(propgrilse = exp(log1SW)/(exp(log1SW)+exp(log2SW))) %>%
  group_by(river_name) %>%
  summarise(propgrilse = mean(propgrilse)) %>%
  arrange(desc(propgrilse))

river_names <- rep(alldatastan$river_name, times = alldatastan$nyears)
#unique_rivers <- unique(river_names) %>% sort
unique_rivers <- propgrilse$river_name

z1mat <- as.matrix(fithierncp_zsdmu, pars = "Z1") 
s1mat <- as.matrix(fithierncp_zsdmu, pars = "S1") 

z1tib2 <- z1mat %>% t %>% as_tibble %>% mutate(river_name = river_names, year = alldatastan$years) %>%
  select(river_name, year, everything())
s1tib2 <- s1mat %>% t %>% as_tibble %>% mutate(river_name = river_names, year = alldatastan$years) %>%
  select(river_name, year, everything())


riverpos <- tibble(year = alldatastan$years, river_name = river_names) %>%
  mutate(pos = row_number()) %>% select(pos, everything())

# matchingyears <- intersect(filter(riverpos, river_name == "Saint-Jean River") %>% pull(year), 
# filter(riverpos, river_name == "Western Arm Brook") %>% pull(year))

riverpairs <- combn(unique_rivers, 2, simplify = FALSE)

river_cor <- function(river1, river2) {
  # river1 <- "Saint-Jean River"
  # river2 <- "Western Arm Brook"
  
  matchingyears <- intersect(filter(riverpos, river_name == river1) %>% pull(year), 
                             filter(riverpos, river_name == river2) %>% pull(year))
  
  #browser()
  z1river1 <- filter(z1tib2, river_name == river1 & year %in% matchingyears) %>% 
    dplyr::select(-river_name, -year) 
  z1river2 <- filter(z1tib2, river_name == river2 & year %in% matchingyears) %>% 
    dplyr::select(-river_name, -year) 
  
  # pivot_longer(z1river1, cols = V1:V7500, names_to = "pos", names_prefix = "^V", values_to = "value") %>%
  #   mutate(year = as.character(year), pos = as.numeric(pos)) %>% ungroup() %>%
  #   ggplot(aes(pos, value, color = year)) + geom_line() + facet_wrap(~year)
  
  var1 <- mapply(var, z1river1)
  var2 <- mapply(var, z1river2)
  
  covz1 <- mapply(cov, z1river1, z1river2) 
  
  # estimating error variances and covariances, this is done by estimating variance across 
  # iterations for each year and then averaging across years, hence the transpose
  t_river1 <- z1river1 %>% t %>% as_tibble
  t_river2 <- z1river2 %>% t %>% as_tibble
  
  errorvar1 <- mapply(var, t_river1) %>% mean
  errorvar2 <- mapply(var, t_river2) %>% mean
  errorcov <- mapply(cov, t_river1, t_river2) %>% mean
  
  uncorrected_cor <- (covz1)/(sqrt(var1)*sqrt(var2))
  corrected_cor <- (covz1 - errorcov)/(sqrt(var1 - errorvar1)*sqrt(var2 - errorvar2))
  uncorrected_cor <- (covz1)/(sqrt(var1)*sqrt(var2))
  
  print(paste0("Number of matching years: ", length(matchingyears)))
  print(paste0(river1," and ", river2, ", Done!"))
  
  return(list(cor = unname(corrected_cor), 
              uncor = unname(uncorrected_cor), 
              matchingyears = matchingyears,
              riverpair = c(river1, river2)))
}

vec2mat <- function(M) {
  b <- matrix(0, 7, 7)
  colnames(b) <- unique_rivers
  rownames(b) <- unique_rivers
  b[lower.tri(b, diag = FALSE)] <- M 
  b <- t(b)
  b[lower.tri(b, diag = FALSE)] <- M 
  diag(b) <- 1
  as.matrix(b)
}


# corr1 <- river_cor("Saint-Jean River", "Western Arm Brook")

rivercors <- lapply(riverpairs, function(x) river_cor(x[1], x[2]))
lapply(rivercors, function(x) anyNA(x$cor))

which(is.na(rivercors[[3]]$cor))
which(is.na(rivercors[[4]]$cor))
which(is.na(rivercors[[6]]$cor))

lapply(rivercors, function(x) quantile(x$cor, c(0.025, 0.5, 0.975), na.rm = TRUE))
medcor <- lapply(rivercors, function(x) quantile(x$cor, 0.5, na.rm = TRUE)) %>% unlist %>% vec2mat()
corlow <- lapply(rivercors, function(x) quantile(x$cor, 0.025, na.rm = TRUE)) %>% unlist() %>% vec2mat()
corhi <- lapply(rivercors, function(x) quantile(x$cor, 0.975, na.rm = TRUE)) %>% unlist() %>% vec2mat()
pvals <- lapply(rivercors, function(x) quantile(x$cor, c(0.025, 0.5, 0.975), na.rm = TRUE)) %>%
  lapply(function(x) sign(x[1]*x[3])) %>% unlist %>% unname() %>% {ifelse(. == 1, 0.001, 0.5)} %>%
  vec2mat()

ci.up <- 0.975
ci.low <- 0.025

medcor_un <- lapply(rivercors, function(x) quantile(x$uncor, 0.5, na.rm = TRUE)) %>% unlist %>% vec2mat()
corlow_un <- lapply(rivercors, function(x) quantile(x$uncor, 0.025, na.rm = TRUE)) %>% unlist() %>% vec2mat()
corhi_un <- lapply(rivercors, function(x) quantile(x$uncor, 0.975, na.rm = TRUE)) %>% unlist() %>% vec2mat()
pvals_un <- lapply(rivercors, function(x) quantile(x$uncor, c(0.025, 0.5, 0.975), na.rm = TRUE)) %>%
  lapply(function(x) sign(x[1]*x[3])) %>% unlist %>% unname() %>% {ifelse(. == 1, 0.001, 0.5)} %>%
  vec2mat()


rivp <- lapply(riverpairs, . %>% word(1) %>% substr(1,3) %>% paste0(., collapse = "-")) %>% unlist
rivpmat <- lapply(rivercors, function(x) R.utils::seqToHumanReadable(x$matchingyears, 1)) %>% 
  unlist() %>% vec2mat()
rivpmat[upper.tri(rivpmat, diag = TRUE)] <- NA
rivpmat

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))


png("figures/corrplot-post-z1-corrected-weakz.png", width = 4, height = 3.5, units = "in", res = 300, pointsize = 7)
corrplot(medcor, method = "color", col = col(200),  
         type = "upper", diag = FALSE,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = pvals, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         mar = c(0, 0, 1.8, 0),
         title = expression(paste("Corrected correlation of posterior estimates of ",italic(Z[1]),
                                  " with weak ",italic(S[1])," and ",,italic(S[2]), " priors"))
)
dev.off()

