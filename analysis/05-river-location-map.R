library(tidyverse)
library(ggmap)
library(mapdata)
library(maps)
library(rgdal)
library(sf)
library(rgeos)
library(ggsn)
library(ggrepel)

# from https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family = "Ubuntu Regular", color = "#22211d"),
      #axis.line = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      #axis.ticks = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.minor = element_blank(),
      # plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      panel.background = element_rect(fill = "#f5f5f2", color = NA), 
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      ...
    )
}

# Based on https://semba-blog.netlify.com/10/20/2018/genetic-connectivity-in-western-indian-ocean-region/
theme_map2 <- function(...) {
  theme_bw() + 
    theme(panel.background = element_rect(fill = "lightblue"),
          axis.text = element_text(size = 11, colour = 1),
          panel.grid = element_line(colour = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          ...
    )
}

if (!file.exists("./data/mapdata/ne_50m_admin_1_states_provinces_lakes.dbf")){
  download.file(file.path('http://www.naturalearthdata.com/http/',
                          'www.naturalearthdata.com/download/50m/cultural',
                          'ne_50m_admin_1_states_provinces_lakes.zip'), 
                f <- tempfile())
  unzip(f, exdir = "./data/mapdata")
  rm(f)
}

region <- readOGR("./data/mapdata", 'ne_50m_admin_1_states_provinces_lakes', encoding='UTF-8')


regions <- subset(region, name %in% c("British Columbia", "Alberta", "Saskatchewan", "Manitoba", "Ontario", "Québec", "New Brunswick", "Prince Edward Island", "Nova Scotia", "Newfoundland and Labrador", "Yukon", "Northwest Territories", "Nunavut")) # region is defined in the first part of the code (see above)

regions2 <- subset(region, name %in% c("Québec", "New Brunswick", "Prince Edward Island", "Nova Scotia", "Newfoundland and Labrador", "Maine", "New Hampshire", "Vermont", "Massachusetts", "Connecticut")) # region is defined in the first part of the code (see above)


r2sf <- st_as_sf(regions2)

regions_sf <- st_as_sf(region) %>%
  filter(region_sub %in% c("Atlantic Canada", "New England", "Québec"))

provlabels <- read_csv("raw-data/prov-labels.csv") 
riverlocs <- read_csv("raw-data/river-locs.csv") 
dput(riverlocs)
sealoc <- tibble(sea = c("Gulf of\nSt. Lawrence", "Gulf of\nMaine", "ATLANTIC\nOCEAN"),
                 long = c(-60.7, -68, -56),
                 lat = c(48.6, 43.4, 44.5))

tibble(river = c("Campbellton", "LaHave", "Nashwaak", 
               "Trinité", "Conne", "Saint-Jean", "Western Arm Brook"), 
     long = c(-54.936111, -64.52, -66.586194, -67.302542, -55.743611, -64.466032, -56.768353),
     lat = c(49.286667, 44.37, 45.965694, 49.418167, 47.866667, 48.767402, 51.190807), 
     yadj = c(-0.4, -0.4, -0.4, 0.4, -0.4, -0.4, 0.4), 
     xadj = c(1, 0, 0, 0, 0, 0, 0))


ggplot(regions_sf) + 
  geom_sf(aes(group = name)) + #, color = "gray90", fill = "gray60") + 
 # geom_polygon(fill = "gray60") +
  xlim(c(-71, -52.5)) + ylim(c(43, 52)) +
 #geom_path(color="white") +
  coord_sf() + 
  guides(fill = FALSE) +
 # geom_sf_text(aes(label = postal)) +
  geom_text(data = provlabels, aes(long, lat, label = province), size = 4.5, inherit.aes = FALSE) +
  geom_point(data = riverlocs, aes(long, lat), inherit.aes = FALSE, size = 2.5) +
  ggrepel::geom_label_repel(data = riverlocs, aes(long, lat, label = river), min.segment.length = 0.3,
                            nudge_y = -0.4, nudge_x = 0.5) +
  geom_text(data = sealoc[1:2,], aes(long, lat, fontface = "italic", label = sea), 
            color = "cornflowerblue", size = 3.5, inherit.aes = FALSE) +
    geom_text(data = sealoc[3,], aes(long, lat, family = "serif", fontface = "bold.italic", label = sea), 
              color = "cornflowerblue", size = 4, inherit.aes = FALSE) +
 # geom_text(data = riverlocs, aes(long, lat), inherit.aes = FALSE, size = 3, label="★", family = "HiraKakuPro-W3") +
  #geom_label(data = riverlocs, aes(long + xadj, lat + yadj, label = river)) +
  ggsn::north(location = "topright",  x.min = -71, x.max = -52.5, y.min = 43, y.max = 52,
              scale = 0.07, symbol = 12) +
  xlab("") + ylab("") + theme_map2() +
  ggsn::scalebar(location = "bottomright", x.min = -71, x.max = -52.5, y.min = 43.2, y.max = 52, 
                 dist = 200, dist_unit = "km", transform = TRUE, border.size = 0.5, height = 0.02,
                 model = "WGS84", st.bottom = TRUE, st.dist = 0.02, st.size = 3)

ggsave("figures/rivers-map.png", width = 7, height = 5)  
ggsave("figures/rivers-map2.png", width = 6.4, height = 4.5)  
ggsave("figures/rivers-map2.pdf", width = 6.4, height = 4.5)  

