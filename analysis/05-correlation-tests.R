# Exploring correlations among trends of S1, smolt, and return abundances

library(tidyverse)
library(corrplot)
source("99-helper-functions.R")

alldatastan <- readRDS("data/alldatastan.rds")
s1quant <- readRDS("data/s1quant.rds")
allreturns <- readRDS("data/allreturns.rds") %>% 
  group_by(river_name) %>%
  mutate(zscoresmolts = scale(log(smolts)), zscoretotal = scale(log(totalreturns))) %>%
  ungroup()
  


river_names <- s1quant %>%
  ungroup() %>%
  select(year, zscore, river_name) %>%
  spread(year, zscore)  %>%
  pull(river_name)

s1mat <- s1quant %>%
  #filter(year <= 2013) %>%
  ungroup() %>%
  select(year, zscore, river_name) %>%
  spread(year, zscore) %>%
  select(-river_name) %>%
  t %>% as.matrix

colnames(s1mat) <- c(word(river_names)[1:6], "WAB")

M <- cor(s1mat, use = "complete.obs", method = c("pearson"))
M

# matrix of the p-value of the correlation
# Pearson's product moment correlation coefficient
p.mat <- cor.mtest(s1mat)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

oldpar <- par() 

png("figures/corr-s1.png", width = 7.5, height = 3, units = "in", res = 300, pointsize = 9.5)
par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
corrplot(M, method = "color", col = col(200),  
         type = "upper", order = "hclust",  diag = FALSE, 
          number.cex = 0.8,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
       #  title = expression(paste("Correlation of standardized Z-scores in ", S[1])),
         mar = c(0, 0, 0, 0)
)

mtext("(a)", line = -2, adj =0.1, cex = 1.6)

corrplot(M, method = "color", col = col(200),  
         type = "upper", order = "hclust",  diag = FALSE, 
         number.cex = 0.8,
         #addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "label_sig", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         #title = expression(paste("Correlation of standardized Z-scores in ", S[1])),
         mar = c(0, 0, 0, 0)
)
  mtext("(b)", line = -2, adj =0.1, cex = 1.6)
dev.off()
    

ggplot(s1quant, aes(year, median, color = river_name)) + 
  geom_line(alpha = 0.5) + geom_point(size = 2) + ylab("S1 estimates") +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1)


ggplot(s1quant, aes(year, median, color = river_name)) + 
  geom_line() + geom_point(size = 2) + ylab("S1 estimates") +
  geom_errorbar(aes(ymin = q05, ymax = q95),width= 0.1) + 
  facet_wrap(~river_name)


### Correlation of smolts and returns

smoltsmat <- allreturns %>%
  select(river_name, year, zscoresmolts) %>%
  spread(year, zscoresmolts) %>%
  ungroup() %>%
  select(-river_name) %>%
  t %>% as.matrix

totalmat <- allreturns %>%
  select(river_name, year, zscoretotal) %>%
  spread(year, zscoretotal) %>%
  ungroup() %>%
  select(-river_name) %>%
  t %>% as.matrix



colnames(smoltsmat) <- river_names
colnames(totalmat) <- river_names


Msmolts <- cor(smoltsmat, use = "complete.obs", method = c("pearson"))
Msmolts
Mtotal <- cor(totalmat, use = "complete.obs", method = c("pearson"))
Mtotal





# matrix of the p-value of the correlation
p.mat.smolts <- cor.mtest(smoltsmat)
p.mat.total <- cor.mtest(totalmat)

tm <- cor.test(smoltsmat[,4], smoltsmat[,5])


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(Msmolts, method = "color", col = col(200),  
         type = "upper", order = "hclust",  diag = FALSE,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         # p.mat = p.mat.smolts, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         mar = c(0, 0, 1.8, 0),
         title = "Correlation of standardized Z-scores of log-transformed smolt abundances"
)

corrplot(Mtotal, method = "color", col = col(200),  
         type = "upper", order = "hclust", diag = FALSE,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat.total, sig.level = 0.01, insig = "blank", # only coloured boxes have significant pairwise p-values.
         # hide correlation coefficient on the principal diagonal
         # diag=FALSE,
         mar = c(0, 0, 1.8, 0),
         title = "Correlation of standardized Z-scores of log-transformed total return abundances"
)

allreturns  %>%
  group_by(river_name) %>%
  mutate(zscoresmolts = scale(smolts), zscoretotal = scale(totalreturns)) %>%
  ggplot(aes(year, smolts, colour = river_name)) + geom_point() + geom_line() +
  facet_wrap(~river_name, scale = "free_y")

allreturns  %>%
  group_by(river_name) %>%
  mutate(zscoresmolts = scale(smolts), zscoretotal = scale(totalreturns)) %>%
  ggplot(aes(year, zscoresmolts, colour = river_name)) + geom_point() + geom_line()

