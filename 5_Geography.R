#####
# Age of Reptiles Analysis
# Voronoi of Fauna
#####

library(jpeg)
library(tidyverse)
library(deldir); library(rgdal)

aor.raw <- readJPEG("AgeOfReptiles.jpg")

aor.raw2 <- bind_rows(
  list(
    (as.data.frame(aor.raw[, , 1]) %>% 
       mutate(y=row_number(), channel = "R") %>% 
       gather(x, value, -y, -channel) %>% 
       mutate(x = as.numeric(gsub("V", "", x)))),
    (as.data.frame(aor.raw[, , 2]) %>% 
       mutate(y=row_number(), channel = "G") %>% 
       gather(x, value, -y, -channel) %>% 
       mutate(x = as.numeric(gsub("V", "", x)))),
    (as.data.frame(aor.raw[, , 3]) %>% 
       mutate(y=row_number(), channel = "B") %>% 
       gather(x, value, -y, -channel) %>% 
       mutate(x = as.numeric(gsub("V", "", x))))
  )
) 

aor <- aor.raw2 %>% 
  spread(channel, value) %>% 
  mutate(shadow = 1 - (R+G+B)/3) %>% 
  mutate(Rp = R, Gp = G, Bp = B) %>% 
  mutate_at(vars(Rp, Gp, Bp), funs( as.character(as.hexmode(round(.*255))))) %>% 
  mutate_at(vars(Rp, Gp, Bp), funs(ifelse(nchar(.) == 1, paste0("0", .), .))) %>% 
  mutate(color = toupper(paste0("#", Rp, Gp, Bp))) %>% 
  select(-Rp, -Gp, -Bp)

###
# Inputs
###

timeline.period <- readRDS("DATA/Timeline.rds")

chart.period <- timeline.period %>% 
  select(label = period, value = x_ef_median) %>% 
  bind_rows(timeline.period %>% 
              select(label = earliest, value = x_ef_earliest) %>% 
              mutate(label = as.character(label))) %>% 
  arrange(value) %>% 
  mutate(label = ifelse("Carboniferous /\n Devonian", "Devonian /\n Carboniferous", label))

####
# Object coords
####

coords_fauna <- read_csv("DATA/Fauna_Coords.csv")
coords_flora <- read_csv("DATA/Flora_Coords.csv")
coords_bg    <- read_csv("DATA/BG_Coords.csv")

coords_object <- coords_fauna %>% 
  mutate(object_type = "fauna") %>% 
  rename(object_name = Fauna) %>% 
  bind_rows(coords_flora %>% 
              select(-note) %>% 
              mutate(object_type = "flora") %>% 
              rename(object_name = Flora)) %>% 
  bind_rows(coords_bg %>% 
              mutate(object_type = "background")) %>% 
  rename(object_id = ID) %>% 
  mutate(y = max(aor.raw2$y)-y)

saveRDS(coords_object, file = "DATA/ObjectCoordinates.RDS")

####
# Voronoi
####

aor.vor <- aor %>%  
  mutate(y = max(y) - y) %>% 
  mutate(year = (362-65)/(max(x))*x + 65) %>% 
  #TIme period, as painted
  mutate(period = case_when(
    year <= timeline.period$ef_earliest[1] ~ timeline.period$period[1],
    year <= timeline.period$ef_earliest[2] ~ timeline.period$period[2], 
    year <= timeline.period$ef_earliest[3] ~ timeline.period$period[3], 
    year <= timeline.period$ef_earliest[4] ~ timeline.period$period[4], 
    TRUE ~ timeline.period$period[5]
  ))

vor_pts <- sp::SpatialPointsDataFrame(cbind(coords_object$x,
                                            coords_object$y),
                                      coords_object, match.ID=TRUE)

SPointsDF_to_voronoi_SPolysDF <- function(sp) {
  
  # tile.list extracts the polygon data from the deldir computation
  vor_desc <- tile.list(deldir(sp@coords[,1], sp@coords[,2]))
  
  lapply(1:(length(vor_desc)), function(i) {
    
    # tile.list gets us the points for the polygons but we
    # still have to close them, hence the need for the rbind
    tmp <- cbind(vor_desc[[i]]$x, vor_desc[[i]]$y)
    tmp <- rbind(tmp, tmp[1,])
    
    # now we can make the Polygon(s)
    Polygons(list(Polygon(tmp)), ID=i)
    
  }) -> vor_polygons
  
  # hopefully the caller passed in good metadata!
  sp_dat <- sp@data
  
  # this way the IDs _should_ match up w/the data & voronoi polys
  rownames(sp_dat) <- sapply(slot(SpatialPolygons(vor_polygons),
                                  'polygons'),
                             slot, 'ID')
  
  sp::SpatialPolygonsDataFrame(SpatialPolygons(vor_polygons),
                               data=sp_dat)
  
}

vor <- SPointsDF_to_voronoi_SPolysDF(vor_pts)
vor_df <- fortify(vor) %>% 
  rename(x = long, y = lat) %>% 
  left_join(data.frame(id = as.character(1:nrow(coords_object)), 
                       object_type = coords_object$object_type,
                       AOR_Fauna = coords_object$object_name,
                       stringsAsFactors = FALSE))

length(unique(vor_df$id))


####
# FIgure out which polygon each point belongs to
####

#Dont actually need the voronoi for this part

coord_map <- coords_object %>% 
  mutate(id = paste0("id", row_number())) %>% 
  gather(coord, value, x, y) %>% 
  select(coord, value, id) %>% 
  spread(id, value)

aor.poly <- readRDS("DATA/AoR_PointsInVor.RDS")

aor.poly2 <- readRDS("DATA/AoR_Voronoi.RDS")

####
# Where on earth are these polygons
####

paleobio.raw <- readRDS("DATA/PaleobioDB_fauna.RDS")

paleobio <- paleobio.raw %>% 
  mutate(cc = ifelse(is.na(cc), "Unknown", cc)) %>% 
  group_by(AOR_Fauna, cc) %>% 
  count() %>% 
  group_by(AOR_Fauna) %>% 
  mutate(Region = ifelse(any(grepl("US|CA|MX|GL", cc)), "North America", cc)) %>% 
  ungroup() %>% 
  mutate(Region = case_when(
    Region == "North America" ~ "North America",
    cc == "Unknown" ~ "Unknown",
    cc %in% c("MA", "TZ", "ZA") ~ "Africa", 
    cc %in% c("AA", "AR", "BO") ~ "South America",
    cc %in% c("CN", "JP", "TJ") ~ "Asia",
    TRUE ~ "Europe"
  )) %>% 
  group_by(AOR_Fauna) %>% 
  filter(n == max(n)) %>% 
  ungroup()

vor_df_pb <- vor_df %>% 
  left_join(paleobio)

ggplot(aes(x=x, y=y),
       data = aor.vor) +
  geom_tile(aes(alpha = shadow), fill = "black") +
  scale_alpha_identity() +
  #geom_point(data = coords_object, mapping = aes(x=x, y=y, color=object_type)) +
  #scale_color_manual(values = c("fauna" = "#ff4040", "flora" = "forestgreen", "background" = "#4040ff")) + 
  geom_map(data=vor_df_pb, map=vor_df,
           aes(x=x, y=y, map_id=id, fill = Region),
           color="#eeeeee", alpha = 0.25, size=1) +
  scale_fill_manual(values = c("North America" = "#4040ff", "Europe" = "forestgreen", 
                               "Africa" = "#ff4040", "Unknown" = "#ff40ff"), na.value = "#cccccc") + 
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, aor.orig.x$value),
                     labels = c(65, aor.orig.x$label),
                     expand = c(0.005, 0.005),
                     name = "Million years ago") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(aor$y))) +
  geom_vline(data = timeline.period, aes(xintercept = x_ef_earliest),
             color="white", alpha = 0.25, size=1.5) +
  labs(title = "Geography of fauna fossils",
       subtitle = "The Age of Reptiles | Rudolph Zallinger",
       caption =  "Blue if any specimens have been found in N.America, otherwise filled by region of discovery\n
       Yale Peabody Museum | paleobiodb.org | @RyanTimpe") +
  coord_fixed() +
  theme_bw() +
  theme( 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 11),
    panel.background=element_blank(),panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank(),
    strip.background = element_rect(fill = "#00436b"),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    plot.title = element_text(color = "#00436b", face = "bold", size = 16),
    plot.subtitle = element_text(color = "#00436b", size = 14),
    plot.caption = element_text(size = 11),
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  )
