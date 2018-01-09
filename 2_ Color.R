#####
# Age of Reptiles Analysis
# Color
#####

library(jpeg)
library(tidyverse)

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

###
# Inputs
###

timeline.period <- readRDS("DATA/Timeline.rds")

chart.period <- timeline.period %>% 
  select(label = period, value = x_ef_median) %>% 
  bind_rows(timeline.period %>% 
              select(label = earliest, value = x_ef_earliest) %>% 
              mutate(label = as.character(label))) %>% 
  arrange(value)

###
# Average COlors
### 
aor.colors.1 <- aor.raw2 %>% 
  spread(channel, value) %>% 
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


aor.color.by.period <- aor.colors.1 %>% 
  group_by(period) %>% 
  summarize_at(vars(R, G, B), mean) %>% 
  gather(channel, value, R, G, B) %>% 
  mutate(value_hex = as.character(as.hexmode(round(value*255)))) %>% 
  rowwise() %>% 
  mutate(value_hex = ifelse(nchar(value_hex) == 1, paste0("0", value_hex), value_hex)) %>% 
  ungroup() %>% 
  select(-value) %>% 
  spread(channel, value_hex) %>% 
  rowwise() %>% 
  mutate(color = toupper(paste0("#", R, G, B))) %>% 
  ungroup()

aor.color.by.ten <- aor.colors.1 %>% 
  mutate(time = ceiling(x/10)) %>% 
  group_by(period, time) %>% 
  summarize_at(vars(R, G, B), mean) %>% 
  gather(channel, value, R, G, B) %>% 
  mutate(value_hex = as.character(as.hexmode(round(value*255)))) %>% 
  rowwise() %>% 
  mutate(value_hex = ifelse(nchar(value_hex) == 1, paste0("0", value_hex), value_hex)) %>% 
  ungroup() %>% 
  select(-value) %>% 
  spread(channel, value_hex) %>% 
  rowwise() %>% 
  mutate(color = toupper(paste0("#", R, G, B))) %>% 
  ungroup()

chart.color.period <- aor.colors.1 %>% 
  left_join(aor.color.by.period %>% select(period, color)) %>% 
  mutate(Version = "(B) Period")

#Color by Time
chart.color.time <- aor.colors.1 %>% 
  mutate(time = ceiling(x/10)) %>% 
  left_join(aor.color.by.ten %>% select(period, time, color)) %>% 
  mutate(Version = "(C) Ten X Pixels")

###
# Faceted Chart
###

chart.color.orig <- aor.colors.1 %>% 
  mutate(Version = "(A) Original") %>% 
  gather(channel, value, R, G, B) %>% 
  mutate(value_hex = as.character(as.hexmode(round(value*255)))) %>% 
  rowwise() %>% 
  mutate(value_hex = ifelse(nchar(value_hex) == 1, paste0("0", value_hex), value_hex)) %>% 
  ungroup() %>% 
  select(-value) %>% 
  spread(channel, value_hex) %>% 
  rowwise() %>% 
  mutate(color = toupper(paste0("#", R, G, B))) %>% 
  ungroup() %>% 
  select(-R, -G, -B)
  
facet_chart <- chart.color.orig %>% 
  bind_rows(chart.color.period %>% select(-R, -G, -B)) %>% 
  bind_rows(chart.color.time %>% select(-R, -G, -B))
  
ggplot(aes(x=x, y=y), 
       data = facet_chart) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     name = "Million years ago") +
  geom_vline(data = timeline.period, aes(xintercept = x_ef_earliest),
             color="white", alpha = 0.5, size=1.5) +
  labs(title = "Average color by selected subsets", 
       subtitle = "The Age of Reptiles | Rudolph Zallinger",
       caption = "Yale Peabody Museum | Ryan Timpe"
  ) +
  coord_fixed() +
  facet_grid(Version~.) +
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
    plot.caption = element_text(size = 11)
  )

####
# Color propensity
####

color_threshold <- 0.34
aor.color.prop <- aor.colors.1 %>% 
  mutate(Rp = R / (R+G+B), 
         Gp = G / (R+G+B),
         Bp = B / (R+G+B)) %>% 
  mutate(color_type = case_when(
    Rp >= color_threshold & Rp > Gp & Rp > Bp ~ "(B) Red",
    Gp >= color_threshold & Gp > Rp & Gp > Bp ~ "(C) Green", 
    Bp >= color_threshold & Bp > Rp & Bp > Gp ~ "(D) Blue",
    TRUE ~ "(A) Grey" #If the 3 channels are similar, the color is a shade of grey
  )) 

aor.color.prop %>% 
  count(color_type)

aor.color.prop2 <- aor.color.prop %>% 
  select(-Rp, -Gp, -Bp) %>% 
  mutate(Ro = R, Go = G, Bo = B) %>% 
  mutate_at(vars(Ro, Go, Bo), funs(as.character(as.hexmode(round(.*255))))) %>% 
  mutate_at(vars(Ro, Go, Bo), funs(ifelse(nchar(.) == 1, paste0("0", .), .))) %>% 
  mutate(color = toupper(paste0("#", Ro, Go, Bo))) %>% 
  select(-Ro, -Go, -Bo)

  
ggplot(aes(x=x, y=y), 
       data = aor.color.prop2) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     name = "Million years ago") +
  geom_vline(data = timeline.period, aes(xintercept = x_ef_earliest),
             color="white", alpha = 0.5, size=1.5) +
  labs(title = "Color propensity in AoR", 
       subtitle = "The Age of Reptiles | Rudolph Zallinger",
       caption = "Yale Peabody Museum | Ryan Timpe"
  ) +
  coord_fixed() +
  facet_grid(color_type~.) +
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
    plot.caption = element_text(size = 11)
  )

aor.color.prop3 <- aor.color.prop2 %>% 
  group_by(x, color_type) %>% 
  mutate( y = row_number()) %>% 
  mutate( y = max(y)-y) %>% 
  ungroup()

ggplot(aes(x=x, y=y), 
       data = aor.color.prop3) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     name = "Million years ago") +
  geom_vline(data = timeline.period, aes(xintercept = x_ef_earliest),
             color="white", alpha = 0.5, size=1.5) +
  labs(title = "Color propensity in AoR - Gravity", 
       subtitle = "The Age of Reptiles | Rudolph Zallinger",
       caption = "Yale Peabody Museum | Ryan Timpe"
  ) +
  coord_fixed() +
  facet_grid(color_type~.) +
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
    plot.caption = element_text(size = 11)
  )



# Color choice by period

aor.color.prop.per <- aor.color.prop3 %>% 
  count(period, color_type) %>% 
  group_by(period) %>% 
  mutate(share = round(n/sum(n), 2)) %>% 
  ungroup()

aor.color.avg.per <- aor.color.prop3 %>% 
  group_by(period, color_type) %>% 
  summarize_at(vars(R, G, B), mean) %>% 
  mutate(brightness = (R+G+B)/3) %>% 
  mutate_at(vars(R, G, B), funs(as.character(as.hexmode(round(.*255))))) %>% 
  mutate_at(vars(R, G, B), funs(ifelse(nchar(.) == 1, paste0("0", .), .))) %>% 
  mutate(color = toupper(paste0("#", R, G, B)))


aor.color.prop4 <- aor.color.prop3 %>% 
  arrange(color_type) %>% 
  group_by(x) %>% 
  mutate( y = row_number()) %>% 
  mutate( y = max(y)-y) %>% 
  ungroup() %>% 
  select(-color) %>% 
  left_join(aor.color.avg.per %>% select(period, color_type, color))

ggplot(aes(x=x, y=y), 
       data = aor.color.prop4) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     name = "Million years ago") +
  geom_vline(data = timeline.period, aes(xintercept = x_ef_earliest),
             color="white", alpha = 0.5, size=1.5) +
  labs(title = "Color propensity in AoR - Displayed by average shade per period", 
       subtitle = "The Age of Reptiles | Rudolph Zallinger",
       caption = "Yale Peabody Museum | Ryan Timpe"
  ) +
  coord_fixed() +
  #facet_grid(color_type~.) +
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
    plot.caption = element_text(size = 11)
  )


