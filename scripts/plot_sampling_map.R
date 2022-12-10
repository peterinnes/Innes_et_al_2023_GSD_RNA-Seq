#### mapping locations of the GSD populations ####
install.packages("usmap")
library(usmap)
install.packages("ggmap")
library(ggmap)
devtools::install_github('oswaldosantos/ggsn')
library(ggsn) #for scalebar
library(tigris)
library(sf)
library(raster)

coordinates <- read.csv("GSD_RNA-seq_population_coords.csv", header=T) %>%
  rename(lon="long")
coordinates$ecotype <- c("non-dune", "non-dune", "non-dune", "dune", "dune", "dune")

bbox <- extent(c(range(coordinates$long)+c(-1,1)*0.05,range(coordinates$lat)+c(-1,1)*0.05))

#map <- get_stamenmap(bbox=c(-105.6117, 37.62805, -105.48,37.79921),
#                     maptype = "terrain", color="bw", zoom=13)
#scale_bar_coords <- data.frame(long=c(-105.5,-105.6), lat=c(37.65,37.8))

map <- get_stamenmap(bbox=c(-105.625, 37.6725, -105.5,37.775),
                     maptype = "terrain", color="bw", zoom=13)

scale_bar_coords <- data.frame(long=c(-105.5,-105.6125), lat=c(37.7, 37.775))

gsd_sampling_map <- ggmap(map) +
  geom_point(data = coordinates, mapping = aes(x=long, y=lat,
                                               shape=ecotype, fill=ecotype),
             size=5) +
  scale_shape_manual(values = c(24,21)) +
  scale_fill_manual(values=c("gold2", "forestgreen")) +
  theme(#legend.position = "none", #c(0.615, .2),
        #legend.key.size = unit(.125, 'in'),
        #legend.background = element_rect(fill = NA),
        #legend.key = element_rect(fill=NA),
        #legend.title = element_blank(),
        legend.position = "bottom",
        legend.background = element_blank(),
        #legend.box.background = element_rect(color = "black"),
        text = element_text(size=24),
        legend.text = element_text(size=18),
        axis.text = element_text(size=18)) +
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_continuous(expand=c(0,0), breaks = c(-105.6, -105.55, -105.5)) +
  scale_y_continuous(expand=c(0,0), breaks = c(37.65, 37.7, 37.75)) +
  labs(x="Long", y="Lat") +
  ggsn::scalebar(data = scale_bar_coords, location = "bottomleft", transform = T,
                 dist = 2, dist_unit = "km", height = .02, st.dist = .05, box.fill = c("black","white"))

gsd_sampling_map
  
ggsave("figures/gsd_sampling_map.pdf", plot = gsd_sampling_map, device = "pdf",
       height = 4.5, units = "in", dpi = 300)
ggsave("figures/gsd_sampling_map.png", plot = gsd_sampling_map, device = "png",
       height = 4.5, units = "in", dpi = 300, bg = "transparent")


# zoomed-out map of colorado with GSDNP outlined
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

#### using MRLCC land cover data, trying to replicate Huang et al 2020 ####
bg_map <- raster("data2/NLCD_2019_Land_Cover_L48_20210604_4yefJHSIZtp5byhkMSjP.tiff")
bbox <- extent(c(range(coordinates$long)+c(-1,1)*0.05,range(coordinates$lat)+c(-1,1)*0.05))
poly <- as(bbox, 'SpatialPolygons')
proj4string(poly) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
poly <- spTransform(poly,proj4string(bg_map))

bg_box <- crop(bg_map, poly)


bg_df <- rasterToPoints(bg_box)
bg_df <- data.frame(bg_df)
colnames(bg_df) <- c("long","lat","lc")
bg_df$lc <- factor(bg_df$lc)


# Color only two land cover types
lc_col <- rep("gray60",length(unique(bg_df$lc)))
lc_col[c(which(levels(bg_df$lc)==31),which(levels(bg_df$lc)==52))] <- c("burlywood1","yellowgreen") # "Barren Land (Rock/Sand/Clay)" and "Shrub/Scrub" 

p <- ggplot(data = bg_df, aes(long, lat)) +
  geom_raster(aes(fill=lc),alpha=.9) +
  theme_bw() +
  scale_fill_manual(values=lc_col,
                    name="Land Cover Classification\nand\nCluster")
