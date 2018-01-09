#####
# Age of Reptiles Analysis
# Corrected time series
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

aor <- aor.raw2 %>% 
  mutate(value_hex = as.character(as.hexmode(round(value*255)))) %>% 
  mutate(value_hex = ifelse(nchar(value_hex) == 1, paste0("0", value_hex), value_hex)) %>% 
  select(-value) %>% 
  spread(channel, value_hex) %>% 
  mutate(color = toupper(paste0("#", R, G, B)))

saveRDS(aor.raw2, "DATA/AoR_DF.rds")

###
# Inputs
###

timeline.period <- data.frame(
  period = c("Cretaceous", "Jurassic", "Triassic", "Permian", "Carboniferous /\n Devonian"),
  earliest = c(144, 206, 248, 290, 362),
  latest = c(65, 144, 206, 248, 290),
  stringsAsFactors = F
) %>% 
  mutate(median = (earliest - latest)/2 + latest) %>% 
  #Painting isnt to scale... effective earliest
  mutate(ef_earliest = c(138, 236, 264, 325, 362),
         ef_latest = c(65, 138, 236, 264, 325)) %>% 
  mutate(ef_median = (ef_earliest - ef_latest)/2 + ef_latest) %>% 
  #Convert to x... this makes the fixed_coord easier
  mutate(x_earliest = (earliest - 65)*max(aor$x)/(362-65),
         x_latest = (latest - 65)*max(aor$x)/(362-65),
         x_median = (median - 65)*max(aor$x)/(362-65),
         x_ef_earliest = (ef_earliest - 65)*max(aor$x)/(362-65),
         x_ef_latest = (ef_latest - 65)*max(aor$x)/(362-65),
         x_ef_median = (ef_median - 65)*max(aor$x)/(362-65)
         )

saveRDS(timeline.period, "DATA/Timeline.rds")

###
# Original AoR Painting
###
aor.orig <- aor %>%
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

aor.orig.x <- timeline.period %>% 
  select(label = period, value = x_median) %>% 
  bind_rows(timeline.period %>% 
              select(label = earliest, value = x_earliest) %>% 
              mutate(label = as.character(label))) %>% 
  arrange(value)

ggplot(aes(x=x, y=y),
       data = aor.orig) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, aor.orig.x$value),
                     labels = c(65, aor.orig.x$label),
                     expand = c(0.005, 0.005),
                     name = "Million years ago") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(data = timeline.period, aes(xintercept = x_earliest),
             color="white", size=1.5) +
  labs(title = "Original Painting - Assuming 238k years/pixel",
       subtitle = "The Age of Reptiles | Rudolph Zallinger",
       caption = "Yale Peabody Museum | Ryan Timpe") +
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
    plot.caption = element_text(size = 11)
  )

####
# AoR if the timeline were proportional
####

#... every pixel is the same lenght of time

timeline.period.2 <- timeline.period %>% 
  #Years per pix
  mutate(MYApP_orig = (earliest - latest) / (x_ef_earliest - x_ef_latest),
         MYApP_total = (362-65)/(1250),
         MYApP_adjust = MYApP_orig/MYApP_total)

timeline.distort <- timeline.period.2 %>% 
  select(period, MYApP_adjust)

aor.distort <- aor.raw2 %>% 
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
  )) %>% 
  left_join(timeline.distort) %>% 
  left_join(timeline.period %>% 
              select(period, x_earliest, x_latest)) %>% 
  #Distort Xs
  mutate(x_latest = ifelse(x_latest == 0, 1, x_latest)) %>% 
  mutate(x_adj = MYApP_adjust * x) %>% 
  group_by(period) %>% 
  mutate(x2 = x_latest + (x_adj - min(x_adj))*(x_earliest - x_latest)/(max(x_adj) - min(x_adj))) %>% 
  ungroup()

aor.distort2 <- aor.distort %>% 
  select(x=x2, y, R, G, B) %>% 
  mutate(x = round(x)) %>% 
  group_by(x, y) %>% 
  summarize_at(vars(R, G, B), mean, na.rm=T) %>% 
  ungroup() %>% 
  complete(x = 1:1250, y) %>% 
  group_by(y) %>% 
  fill(R, G, B) %>% 
  ungroup() %>% 
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

aor.distort3 <- aor.distort2 %>% 
  mutate(year = (362-65)/(max(x))*x + 65) %>%   
  mutate(period = case_when(
    year <= timeline.period$earliest[1] ~ timeline.period$period[1],
    year <= timeline.period$earliest[2] ~ timeline.period$period[2], 
    year <= timeline.period$earliest[3] ~ timeline.period$period[3], 
    year <= timeline.period$earliest[4] ~ timeline.period$period[4], 
    TRUE ~ timeline.period$period[5]
  ))

aor.distort.x <- timeline.period %>% 
  select(label = period, value = x_median) %>% 
  bind_rows(timeline.period %>% 
              select(label = earliest, value = x_earliest) %>% 
              mutate(label = as.character(label))) %>% 
  arrange(value)

# ggplot(aes(x=x, y=y), 
#        data = aor.distort2) +
#   geom_tile(aes(fill = color)) +
#   scale_fill_identity() +
#   scale_x_continuous(limits = c(0, 1250),
#                      breaks = c(0, aor.distort.x$value),
#                      labels = c(65, aor.distort.x$label),
#                      name = "Million years ago") +
#   geom_vline(data = timeline.period, aes(xintercept = x_earliest),
#              color="white", size=1.5) +
#   labs(title = "Fixed Timeline Painting - 238k years / pixel"#, 
#        #subtitle = "The Age of Reptiles | Rudolph Zallinger",
#        #caption = "Yale Peabody Museum | Ryan Timpe"
#        ) +
#   coord_fixed() +
#   theme_bw() +
#   theme( 
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.title.y = element_blank()
#   )

###
#Merge into one files
###

aor.comp <- aor.orig %>% 
  select(x, y, color, period) %>% 
  mutate(Version = " Original") %>% 
  bind_rows(aor.distort3 %>% 
              select(x, y, color, period) %>% 
              mutate(Version = "Corrected"))

compare.vert <- timeline.period %>% 
  select(vert = x_ef_earliest) %>% 
  mutate(Version = " Original") %>% 
  bind_rows(timeline.period %>% 
              select(vert = x_earliest) %>% 
              mutate(Version = "Corrected")) %>% 
  filter(vert != max(vert))

mya_p_pixel <- timeline.period.2 %>% 
  select(period, MYApP_orig) %>% 
  mutate(KYApP = MYApP_orig*1000) %>% 
  mutate(label = paste0(round(KYApP, 0), "k yr/p")) %>% 
  left_join(timeline.period %>% select(period, x_ef_latest)) %>% 
  mutate(Version = " Original")

ggplot(aes(x=x, y=y), 
       data = aor.comp) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1250),
                     breaks = c(0, aor.distort.x$value),
                     labels = c(65, aor.distort.x$label),
                     name = "Million years ago") +
  geom_vline(data = compare.vert, aes(xintercept = vert),
             color="white", alpha = 0.5, size=1.5) +
  geom_label(data = mya_p_pixel, y=160,
             aes(label = label, x = x_ef_latest+20),
             hjust = 0, vjust=1, alpha = 0.75) +
  labs(title = "Corrected Timeline Comparison - 238k years / pixel", 
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

#Cretaceous + Jurassic plot only

ggplot(aes(x=x, y=y), 
       data = aor.comp %>% filter(period %in% timeline.period$period[1:2])) +
  geom_tile(aes(fill = color)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 719),
                     breaks = c(0, aor.distort.x$value[1:4]),
                     labels = c(65, aor.distort.x$label[1:4]),
                     name = "Million years ago") +
  geom_vline(data = compare.vert[c(1,2,5,6), ], aes(xintercept = vert),
             color="white", alpha = 0.5, size=1.5) +
  geom_label(data = mya_p_pixel %>% filter(period %in% timeline.period$period[1:2]), 
             y=160,
             aes(label = label, x = x_ef_latest+20),
             hjust = 0, vjust=1, alpha = 0.75) +
  labs(title = "Corrected Timeline Comparison - 238k years / pixel \nCretaceous + Jurassic only", 
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

            