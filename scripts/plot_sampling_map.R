#### mapping locations of the GSD populations ####
install.packages("usmap")
library(usmap)
install.packages("ggmap")
library(ggmap)
devtools::install_github('oswaldosantos/ggsn')
library(ggsn) #for scalebar
library(tigris)
library(sf)

register_google(key = "AIzaSyCsWRRs09XMpil8vh7HD5zGiXXd46QQ93U", write = T)

coordinates <- read.csv("GSD_RNA-seq_population_coords.csv", header=T) %>%
  rename(lon="long")
coordinates$Ecotype <- c("non-dune", "non-dune", "non-dune", "dune", "dune", "dune")

bbox <- make_bbox(lon = coordinates$long, lat = coordinates$lat)
map <- get_stamenmap(bbox=c(-105.675, 37.65, -105.475,37.85), maptype = c("terrain"), zoom=13)

scale_bar_coords <- data.frame(long=c(-105.5,-105.6), lat=c(37.67,37.8))

gsd_sampling_map <- ggmap(map) +
  geom_point(data = coordinates, mapping = aes(x=long, y=lat, shape=Ecotype), size=10, alpha=.5) +
  scale_shape_manual(values = c(17,19)) +
  theme(legend.position = c(0.25, .175),
        legend.key.size = unit(.5, 'cm'),
        legend.title = element_blank(),
        text = element_text(size=36),
        axis.text = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="Longitude", y="Latitude") +
  ggsn::scalebar(data = scale_bar_coords, location = "bottomright", transform = T,
                 dist = 1, dist_unit = "km", height = .01, st.dist = .05, box.fill = c("black","white"))

gsd_sampling_map
  
ggsave("figures/gsd_sampling_map.pdf", plot = gsd_sampling_map, device = "pdf", width = 6, height = 8, units = "in", dpi = 600)
ggsave("figures/gsd_sampling_map.png", plot = gsd_sampling_map, device = "png", width = 8, height = 10, units = "in", dpi = 600)


# zoomed-out map of colorado with GSDNP outlined?
x_coords <- c(-105.675, -105.675, -105.475, -105.475, -105.675)
y_coords <- c(37.65, 37.85, 37.85, 37.65, 37.65)

poly1 <- sp::Polygon(cbind(x_coords,y_coords))
poly1 <- sp::Polygons(list(poly1), ID = "A")
poly1 <- sp::SpatialPolygons(list(poly1))

colorado_geo <- tigris::states(class="sf") %>%
  subset(NAME=="Colorado")

par(bg=NA)
plot (x)
outset_map <- tmap::tm_shape(colorado_geo, projection = 5070) +
  tm_polygons() +
  tm_shape(poly1, projection = 5070) +
  tm_polygons(col = "red", alpha = 0.5) +
  tm_layout(frame = F, bg.color = "transparent") 


png("figures/gsdnp_colorado.png", width = 8, height = 10, units = "in", res = 300, bg = "transparent")
outset_map
dev.off()