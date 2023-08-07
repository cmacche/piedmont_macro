# library -----------------------------------------------------------------

library(MASS)
library(MuMIn)
library(tidyverse)
library(vegan)

# load data -------------------------------------------
df_wq = read_csv("data_raw/water_quality.csv")
df_invert = read_csv("data_raw/macroinvertebrate.csv")


# format ------------------------------------------------------------------


df_invert = filter(df_invert, EcoRegion == "P") %>% 
  mutate(Order = case_when(Order == "CO" ~ "COLEOPTERA",
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
         Class = case_when(Class == "CR" ~ "CRUSTACEA",
                           Class == "GA" ~ "GASTROPODA",
                           Class == "OL" ~ "OLIGOCHAETA",
                           Class == "PE" ~ "PLECYPODA",
                           TRUE ~ as.character(Class))) %>% 
  filter(!(County %in% c("York, SC",
                         "(Danville), VA",
                         "Tallapoosa, AL",
                         "Laurens, SC",
                         "Newberry, SC")))


df_invert = subset(df_invert, select = c("Date","Latitude", "Longitude", 
                                         "Genus",
                                         "Family", "Order",
                                         "Class", "Abundance"))
df_table = df_invert %>% 
  group_by(Class, Order) %>% 
  summarize(n_na_family = sum(is.na(Family)) + sum(Family == "Unknown", na.rm = T),
            n_na_genus = sum(is.na(Genus)) + sum(Genus == "Unknown", na.rm = T),
            n = n(),
            p_na_family = round((n_na_family / n) * 100, 2),
            p_na_genus = round((n_na_genus / n) * 100, 2)) %>%
  ungroup() %>% 
  filter(p_na_family <= 10,
         p_na_genus <= 10)

df_wq = df_wq %>% 
  filter(EcoRegion == "P",
         Sp_Cond != 0,
         pH_SU != 0,
         Diss_Oxy != 0,
         Temp_C != 0) %>% 
  dplyr::select(Date,
                Water_Class,
                Latitude,
                Longitude,
                Temp_C,
                Sp_Cond,
                pH_SU,
                Diss_Oxy) %>% 
  mutate(site_id = paste0(round(Latitude, 4),
                          round(Longitude, 4)),
         Date = as.Date(Date, format = "%m/%d/%y"))


# Diversity Class -----------------------------------------------------------

## pick classes with < 10% unknowns in genus
class_name <- df_table %>% 
  filter(Class != "UNKNOWN") %>% 
  drop_na(Class) %>% 
  pull(Class) %>% 
  unique()

df_class0 = df_invert %>% 
  filter(Class %in% class_name) %>% 
  mutate(site_id = paste0(round(Latitude, 4),
                          round(Longitude, 4)),
         Date = as.Date(Date, format = "%m/%d/%y"))

df_class = df_class0 %>% group_by(site_id,
                                  Date,
                                  Latitude, 
                                  Longitude,
                                  Class) %>% 
  summarise(Div = diversity(Abundance,
                            index = "shannon"))

df_invertwq = df_class %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()

# analysis ----------------------------------------------------------------

pre_invertwq = df_invertwq %>%
  filter(Date <= median(df_invertwq$Date))

## pick latest sampling each site
df_prem <- pre_invertwq %>% 
  group_by(site_id) %>% 
  slice(which.max(Date)) %>% 
  ungroup()

df_prem = df_prem %>% 
  dplyr::select(c(Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy))

invertwq_glm = lm(Div~ .,
                  data = df_prem,
                  na.action = "na.fail")

invertwq_dredge = dredge(invertwq_glm)

## post sample
post_invertwq = df_invertwq %>%
  filter(Date > median(df_invertwq$Date))

