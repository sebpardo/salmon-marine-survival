library(tidyverse)
library(corrplot)
#source("analysis/99-helper-functions.R")

fithierncp <- readRDS("data/stanfit_hier_ncp.rds") # hierarchical
alldatastan <- readRDS("data/alldatastan.rds")

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

z1mat <- as.matrix(fithierncp, pars = "Z1") 
s1mat <- as.matrix(fithierncp, pars = "S1") 

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


corrplot(medcor, method = "circle", col = col(200), 
         type = "upper", diag = FALSE,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = pvals, sig.level = 0.01, insig = "blank", 
         # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         plotCI = "circle", lowCI.mat = corlow, uppCI.mat = corhi,
         mar = c(0, 0, 1.8, 0),
         title = expression(paste("Corrected correlation of posterior estimates of ",italic(Z[1])))
         )

png("figures/corrplot-post-z1-corrected-all.png", width = 4, height = 3, units = "in", res = 300, pointsize = 7)
corrplot(medcor, method = "color", col = col(200),  
         type = "upper", diag = FALSE,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         # p.mat = p.mat.smolts, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         mar = c(0, 0, 1.8, 0),
         title = expression(paste("Corrected correlation of posterior estimates of ",italic(Z[1])))
         )
dev.off()

png("figures/corrplot-post-z1-corrected.png", width = 4, height = 3, units = "in", res = 300, pointsize = 7)
corrplot(medcor, method = "color", col = col(200),  
         type = "upper", diag = FALSE,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = pvals, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         mar = c(0, 0, 1.8, 0),
         title = expression(paste("Corrected correlation of posterior estimates of ",italic(Z[1])))
         )
  dev.off()

png("figures/corrplot-post-z1-uncorrected.png", width = 4, height = 3, units = "in", res = 300, pointsize = 7)
corrplot(medcor_un, method = "color", col = col(200),  
         type = "upper", diag = FALSE,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = pvals_un, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         mar = c(0, 0, 1.8, 0),
         title = expression(paste("UNCORRECTED correlation of posterior estimates of ",italic(Z[1])))
         )
dev.off()


corrplot(medcor, method = "circle", col = col(200),  
         type = "upper", diag = FALSE, 
         number.cex = 0.8,
         #addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         plotCI = "circle", lowCI.mat = corlow, uppCI.mat = corhi,
         # Combine with significance
         p.mat = pvals, sig.level = 0.01, insig = "label_sig", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         #title = expression(paste("Correlation of standardized Z-scores in ", S[1])),
         mar = c(0, 0, 0, 0))

