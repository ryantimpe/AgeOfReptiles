#####
# Age of Reptiles Analysis
# Shadows
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
# Light + Dark
### 
aor.shd.1 <- aor.raw2 %>% 
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

chart.shd <- aor.shd.1 %>% 
  mutate(total = R+G+B) %>% 
  mutate(shadow = 1- (total / max(total)))

chart.shd2 <- chart.shd %>% 
  mutate(shadow = shadow ^(2),
         Version = "  Shadow ^2") %>% 
  mutate(shadow = shadow / max(shadow)) %>% 
  mutate(shadow_reg = floor(shadow*5)) %>% 
  #filter(shadow_reg < max(shadow_reg))
  mutate(shadow_reg = ifelse(shadow_reg >= 3, "3-5", as.character(shadow_reg))) 

ggplot(aes(x=x, y=y), 
       data = chart.shd2) +
  geom_tile(aes(alpha = shadow), fill = "#000000") +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     name = "Million years ago") +
  geom_vline(data = timeline.period[1:4, ], aes(xintercept = x_ef_earliest),
             color="#00436b", alpha = 0.5, size=1.5) +
  labs(title = "Shading by intensity", 
       subtitle = "The Age of Reptiles | Rudolph Zallinger",
       caption = "Yale Peabody Museum | Ryan Timpe"
  ) +
  coord_fixed() +
  facet_grid(shadow_reg~.) +
  theme_bw() +
  theme( 
    legend.position = "none",
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

###
# Reproduce the undercoat
### 

color_group_n <- 20
chart.uc <- aor.shd.1 %>% 
  mutate(Ru = R, Gu = G, Bu = B) %>% 
  mutate_at(vars(Ru, Gu, Bu), funs(
    ifelse(. == 1, color_group_n-1, round(.*color_group_n))
  )) %>% 
  mutate(
    del_RG = Ru - Gu,
    del_RB = Ru - Bu,
    del_GB = Gu - Bu
  ) %>% 
  mutate(
    color_group = paste(del_RG, del_RB, del_GB)
  ) %>% 
  select(-dplyr::starts_with("del_")) %>% 
  mutate(shadow = ((1-R)+(1-G)+(1-B))/3)

test <- chart.uc %>% 
  count(color_group)

chart.uc2 <- chart.uc %>% 
  mutate(shadow_group = case_when(
    shadow < 0.2 ~ "Light",
    shadow < 0.5 ~ "Medium",
    TRUE ~ "Dark"
  )) %>% 
  group_by(color_group, shadow_group, period) %>% 
  mutate(under = shadow > min(shadow)*1.5) %>% 
  ungroup() %>% 
  mutate_at(vars(R, G, B), funs(as.character(as.hexmode(.*255)))) %>% 
  mutate_at(vars(R, G, B), funs(ifelse(nchar(.) == 1, paste0("0", .), .))) %>% 
  mutate(color = toupper(paste0("#", R, G, B)))


chart.period.jt <- chart.period[2:6, ]
timeline.period.jt <- timeline.period[2, ]

ggplot(aes(x=x, y=y), 
       data = chart.uc2 %>% filter(under, period %in% c("Jurassic", "Triassic"))) +
  geom_tile(aes(fill = "black", alpha = shadow^.75)) +
  scale_fill_identity() +
  scale_alpha_identity() +
  scale_x_continuous(#limits = c(0, 1250),
    breaks = c(0, chart.period.jt$value),
    labels = c(65, chart.period.jt$label),
    expand = c(0.005, 0.005),
    name = "Million years ago") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(data = timeline.period.jt, aes(xintercept = x_ef_earliest),
             color="#00436b", alpha = 0.8, size=1.5) +
  labs(title = "Reconstructed underpainting - Jurassic & Triassic", 
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

ggplot(aes(x=x, y=y), 
       data = chart.uc2 %>% filter(under)) +
  geom_tile(aes(fill = "black", 
                alpha = shadow^.75, .25)) +
  scale_fill_identity() +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     expand = c(0.005, 0.005),
                     name = "Million years ago") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(data = timeline.period[1:4, ], aes(xintercept = x_ef_earliest),
             color="#00436b", alpha = 0.8, size=1.5) +
  labs(title = "Reconstructed underpainting", 
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

ggplot(aes(x=x, y=y), 
       data = chart.uc2) +
  geom_tile(aes(fill = ifelse(under, "black", color), 
                alpha = ifelse(under, shadow, .8))) +
  scale_fill_identity() +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, chart.period$value),
                     labels = c(65, chart.period$label),
                     expand = c(0.005, 0.005),
                     name = "Million years ago") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(data = timeline.period[1:4, ], aes(xintercept = x_ef_earliest),
             color="white", alpha = 0.25, size=1.5) +
  labs(title = "Reconstructed undercoat overlaid on the original", 
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

