#####
# Age of Reptiles
# Data from Paleobio DB
#####

library(tidyverse)

#Animals in AoR
aor_fauna_list <- list()
aor_fauna_list[["carboniferous"]] <- c("meganeuropsis", "eryops", "eogyrinus",
                                       "diplovertebron", "cheirolepis", "eusthenopteron")

aor_fauna_list[["permian"]] <- c("edaphosaurus", "dimetrodon", "limnoscelis", 
                                 "seymouria", "araeoscelis", "varanosaurus", 
                                 "ophiacodon", "sphenacodon")

aor_fauna_list[["triassic"]] <- c("cynognathus", "saltoposuchus", "plateosaurus",
                                  "podokesaurus")

aor_fauna_list[["jurassic"]] <- c("allosaurus", "camptosaurus", "stegosaurus",
                                  "brontosaurus", "rhamphorhynchus",
                                  "archaeopteryx")

aor_fauna_list[["cretaceous"]] <- c("pteranodon", "edmontosaurus", "triceratops", 
                                    "tyrannosaurus", "ankylosaurus", "struthiomimus",
                                    "cimolestes")

saveRDS(aor_fauna_list, file = "DATA/fauna_names.RDS")

#Function to download
# One at a time for now to preserve genus names...
# ... the DB may use updated genus names
create.pbdb <- function(aor_genus){
  
  pbdb.url <- paste0("https://paleobiodb.org/data1.2/occs/list.csv?base_name=", 
                     aor_genus,
                     "&show=class,coll,coords,loc") #Class, collection, coordinates, location
  
  pbdb.df <- read_csv(url(pbdb.url))
  
  pbdb.df <- pbdb.df %>% 
    mutate(AOR_Fauna = aor_genus)

  return(pbdb.df)
}

aor_fauna <- lapply(unlist(aor_fauna_list), create.pbdb)

aor_fauna2 <- lapply(aor_fauna, function(x){
  df <- x %>% 
    mutate_if(is.numeric, as.character)
  return(df)
})

aor_fauna_df <- bind_rows(aor_fauna2) %>% 
  select(AOR_Fauna, name = accepted_name, dplyr::contains("interval"), dplyr::contains("_ma"), 
         lng, lat, cc, state, county,
         phylum:genus)

saveRDS(aor_fauna_df, file = "DATA/PaleobioDB_fauna.RDS")

####
# Some quick analysis
####
aor_fauna_df %>% 
  count(cc) %>% 
  arrange(desc(n))

aor_fauna_df %>% 
  group_by(AOR_Fauna) %>% 
  summarize(in_us = "US" %in% cc) %>% 
  filter(!in_us) %>% 
  pull(AOR_Fauna)
  
