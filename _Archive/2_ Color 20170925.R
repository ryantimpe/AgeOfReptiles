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

aor.color.by.total <- aor.colors.1 %>% 
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

aor.color.by.group <- aor.colors.1 %>% 
  mutate(groupx = ceiling(x/10),
         groupy = ceiling(y/10)) %>% 
  group_by(period, groupx, groupy) %>% 
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

#Color by Period
# chart.color.total <- aor.colors.1 %>% 
#   mutate(id = 1) %>% 
#   left_join(aor.color.by.total %>% select(color) %>% mutate(id=1)) %>% 
#   mutate(Version = "(B) Total") %>% 
#   select(-id)

chart.color.period <- aor.colors.1 %>% 
  left_join(aor.color.by.period %>% select(period, color)) %>% 
  mutate(Version = "(B) Period")

ggplot(aes(x=x, y=y), 
       data = chart.color.period) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     name = "Million years ago") +
  geom_vline(data = timeline.period, aes(xintercept = x_ef_earliest),
             color="white", alpha = 0.5, size=1.5) +
  labs(title = "Average color by Period", 
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



#Color by Time
chart.color.time <- aor.colors.1 %>% 
  mutate(time = ceiling(x/10)) %>% 
  left_join(aor.color.by.ten %>% select(period, time, color)) %>% 
  mutate(Version = "(C) Ten X Pixels")

chart.color.group10 <- aor.colors.1 %>% 
  mutate(groupx = ceiling(x/10),
         groupy = ceiling(y/10)) %>% 
  left_join(aor.color.by.group %>% select(period, groupx, groupy, color)) %>% 
  mutate(Version = "(D) Ten XY Pixels")

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
  #bind_rows(chart.color.total %>% select(-R, -G, -B)) %>% 
  bind_rows(chart.color.period %>% select(-R, -G, -B)) %>% 
  bind_rows(chart.color.time %>% select(-R, -G, -B)) %>% 
  bind_rows(chart.color.group10 %>% select(-R, -G, -B))
  
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
  labs(title = "Average color", 
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
# Relative color
####

rel.color.by.period <- aor.colors.1 %>% 
  group_by(period) %>% 
  summarize_at(vars(R, G, B), mean) %>% 
  ungroup() %>% 
  gather(channel, value, R, G, B) %>% 
  mutate(value_hex = as.character(as.hexmode(round(value*255)))) %>% 
  rowwise() %>% 
  mutate(value_hex = ifelse(nchar(value_hex) == 1, paste0("0", value_hex), value_hex)) %>% 
  ungroup() %>% 
  select(-value) %>% 
  rowwise() %>% 
  mutate(color = case_when(
    channel == "R" ~ toupper(paste0("#", value_hex, "00", "00")),
    channel == "G" ~ toupper(paste0("#", "00", value_hex, "00")),
    channel == "B" ~ toupper(paste0("#", "00", "00", value_hex)),
    TRUE ~ "#000000"
    )
  )%>% 
  ungroup()

chart.rel.period <- aor.colors.1 %>% 
  left_join(rel.color.by.period %>% select(channel, period, color))

ggplot(aes(x=x, y=y), 
       data = chart.rel.period) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     name = "Million years ago") +
  geom_vline(data = timeline.period, aes(xintercept = x_ef_earliest),
             color="white", alpha = 0.5, size=1.5) +
  labs(title = "Average color by Period", 
       subtitle = "The Age of Reptiles | Rudolph Zallinger",
       caption = "Yale Peabody Museum | Ryan Timpe"
  ) +
  coord_fixed() +
  facet_grid(channel~.) +
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

#Original

rel.color.by.orig <- aor.colors.1 %>% 
  mutate(R = as.character(as.hexmode(round(R*255))),
         R = ifelse(nchar(R) == 1, paste0("0", R), R)) %>% 
  mutate(color = toupper(paste0("#", R, "00", "00"))) %>% 
  select(-R, -B, -G) %>% 
  mutate(channel = "  R") %>% 
  bind_rows(
    aor.colors.1 %>% 
      mutate(G = as.character(as.hexmode(round(G*255))),
             G = ifelse(nchar(G) == 1, paste0("0", G), G)) %>% 
      mutate(color = toupper(paste0("#", "00", G, "00"))) %>% 
      select(-R, -B, -G) %>% 
      mutate(channel = " G")
  ) %>% 
  bind_rows(
    aor.colors.1 %>% 
      mutate(B = as.character(as.hexmode(round(B*255))),
             B = ifelse(nchar(B) == 1, paste0("0", B), B)) %>% 
      mutate(color = toupper(paste0("#", "00", "00", B))) %>% 
      select(-R, -B, -G) %>% 
      mutate(channel = "B")
  )


ggplot(aes(x=x, y=y), 
       data = rel.color.by.orig) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     name = "Million years ago") +
  geom_vline(data = timeline.period, aes(xintercept = x_ef_earliest),
             color="white", alpha = 0.5, size=1.5) +
  labs(title = "Average color by Period", 
       subtitle = "The Age of Reptiles | Rudolph Zallinger",
       caption = "Yale Peabody Museum | Ryan Timpe"
  ) +
  coord_fixed() +
  facet_grid(channel~.) +
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
