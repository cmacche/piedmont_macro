
# library -----------------------------------------------------------------

library(MASS)
library(MuMIn)
library(tidyverse)

# load data -------------------------------------------
df_wq = read_csv("data_raw/water_quality.csv")
df_invert = read_csv("data_raw/macroinvertebrate2.csv")


# format ------------------------------------------------------------------


df_invert = filter(df_invert, EcoRegion == "P") %>% 
  mutate(order = case_when(Order == "CO" ~ "COLEOPTERA",
                           Order %in% c("DI", "DIM") ~ "DIPTERA",
                           Order == "OD" ~ "ODONATA",
                           Order == "GA" ~ "COLEOPTERA",
                           Order == "ME" ~ "MEGALOPTERA",
                           Order == "HE" ~ "HEMIPTERA",
                           Order == "EP" ~ "EPHEMEROPTERA",
                           Order == "PL" ~ "PLECOPTERA",
                           Order == "TR" ~ "TRICHOPTERA",
                           Order == "OT" ~ "OTHER_TAXA",
                           TRUE ~ as.character(Order)),
         class = case_when(Class == "CR" ~ "CRUSTACEA",
                           Class == "GA" ~ "GASTROPODA",
                           Class == "OL" ~ "OLIGOCHAETA",
                           Class == "PE" ~ "PLECYPODA",
                           TRUE ~ as.character(Class))) %>% 
  filter(!(County %in% c("York, SC",
                         "(Danville), VA",
                         "Tallapoosa, AL",
                         "Laurens, SC",
                         "Newberry, SC")))

df_invert %>% 
  filter(order == "ODONATA") %>% 
  group_by(Latitude, Longitude) %>% 
  summarize(nsp = n_distinct(ScientificName)) %>% 
  pull(nsp) %>% 
  range()

df_table <- df_invert %>% 
  group_by(class, order) %>% 
  summarize(n_na_family = sum(is.na(Family)) + sum(Family == "Unknown", na.rm = T),
            n_na_latin = sum(is.na(ScientificName)) + sum(ScientificName == "Unknown", na.rm = T),
            n = n(),
            p_na_family = round(n_na_family / n, 2),
            p_na_latin = round(n_na_latin / n, 2))

