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
            p_na_genus = round(n_na_genus / n * 100, 2))

df_table = filter(df_table, p_na_family <= 10)
df_table = filter(df_table, p_na_genus <= 10)

df_wq = filter(df_wq, EcoRegion == "P") %>% 
  filter(!(Sp_Cond == 0)) %>% 
  filter(!(pH_SU == 0)) %>% 
  filter(!(Diss_Oxy == 0)) %>% 
  filter(!(Temp_C == 0))
df_wq = subset(df_wq, select = c( "Date","Water_Class", 
                                  "Latitude","Longitude",
                                  "Temp_C","Sp_Cond", "pH_SU",
                                  "Diss_Oxy"))


df_wq = as.data.frame(unclass(df_wq),                    
                      stringsAsFactors = TRUE)


# Diversity Class -----------------------------------------------------------
table(df_table$Class)
df_class = df_invert %>% filter(Class %in% c("ARACHNIDA", "BIVALVIA",
                                 "CLITELLATA", "ENOPLA", "GASTROPODA",
                                 "HYDROZOA","INSECTA", "MAXILLOPODA",
                                 "OPHIURIOIDEA","POLYCHAETA", "TURBELLARIA"))




df_class2 = df_class %>% group_by(Date,
                                  Latitude, 
                                  Longitude,Class) %>% 
  summarise(Div = diversity(Abundance))

df_invertwq = merge(df_class2, df_wq, by = c("Date", 
                                             "Latitude","Longitude"))



pre_invertwq = df_invertwq %>% filter(Date <= median(df_invertwq$Date))

post_invertwq = df_invertwq %>% filter(Date > median(df_invertwq$Date))


pre_invertwq = subset(pre_invertwq, 
                      select = c("Div", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
invertwq_glm = lm(Div~ ., data = pre_invertwq, na.action = "na.omit")

invertwq_dredge = dredge(invertwq_glm)

