#### mapping locations of the GSD populations ####
install.packages("usmap")
library(usmap)
install.packages("ggmap")
library(ggmap)
devtools::install_github('oswaldosantos/ggsn') 
library(ggsn) #for scalebar

register_google(key = "AIzaSyCsWRRs09XMpil8vh7HD5zGiXXd46QQ93U", write = T)

coordinates <- read.csv("GSD_RNA-seq_population_coords.csv", header=T)
coordinates$Ecotype <- c("non-dune", "non-dune", "non-dune", "dune", "dune", "dune")

bbox <- make_bbox(lon = coordinates$lon, lat = coordinates$lat)
map <- get_stamenmap(bbox=c(-105.675, 37.65, -105.475,37.85), maptype = c("terrain"), zoom=12)

gsd_sampling_map <- ggmap(map) +
  geom_point(data = coordinates, mapping = aes(x=lon, y=lat, shape=Ecotype), size=4) +
  scale_shape_manual(values = c(24,21)) +
  theme(legend.position = c(0.2, .175),
        legend.key.size = unit(.5, 'cm'),
        legend.title = element_blank(),
        text = element_text(size=14),
        axis.text = element_text(size=10)) +
  labs(x="Longitude", y="Latitude")
  




