# library -----------------------------------------------------------------
install.packages("data.table")             # Install & load data.table package
library(data.table)
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

#pick latest sample#
df_postm <- post_invertwq %>% 
  group_by(site_id) %>% 
  slice(which.max(Date)) %>% 
  ungroup()

df_postm = df_postm %>% 
  dplyr::select(c(Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy))


#prediction model#

sc_lm_mod = lm(Div ~ Sp_Cond, data = df_prem)

d0 = df_postm %>% dplyr::select(Div, Sp_Cond)

sc_pred = predict(sc_lm_mod, newdata = d0) %>% exp()

chisq.test(df_postm$Div, sc_pred)

cor.test(df_postm$Div, sc_pred, use = "everything")

#all as the difference of AIC is less than 2#
wq_lm_mod = lm(Div ~ pH_SU + Temp_C + Sp_Cond + Diss_Oxy, data = df_prem)

d0 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy)

wq_pred = predict(wq_lm_mod, newdata = d0) %>% exp()

chisq.test(df_postm$Div, wq_pred)

cor.test(df_postm$Div, wq_pred, use = "everything")



table(df_class$Class)


# arachnida ---------------------------------------------------------------

df_ara = filter(df_class, Class == "ARACHNIDA")

df_invertwq = df_ara %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()

pre_invertwq = df_invertwq %>%
  filter(Date <= median(df_invertwq$Date))

## pick latest sampling each site
df_prem <- pre_invertwq %>% 
  group_by(site_id) %>% 
  slice(which.max(Date)) %>% 
  ungroup()

df_prem = df_prem %>% 
  dplyr::select(c(Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy))

invertwq_lm = lm(Div~ .,
                  data = df_prem,
                  na.action = "na.fail")

invertwq_dredge = dredge(invertwq_lm)

#This model does not work#


# loop --------------------------------------------------------------------

for (i in 1:length(class_name)) {
  fit[[i]] <- df_invertwq %>%
    filter(Date <= median(df_invertwq$Date)) %>% 
    filter(Class == class_name[i]) %>% 
    group_by(site_id) %>% 
    slice(which.max(Date)) %>% 
    ungroup() %>% 
    dplyr::select(c( Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy)) %>% 
    lm(Div ~ Sp_Cond)
    
}
fit

lm(fit[[2]]$Div ~ fit[[2]]$Sp_Cond)

list_beta0 <- lapply(species, FUN = function(x) {
  pre <- df_invertwq %>%
    filter(Date <= median(df_invertwq$Date)) %>% 
    filter(Class == x) %>% group_by(site_id) %>% 
    slice(which.max(Date)) %>% 
    ungroup() %>% 
  return(pre)
})


# bivalvia ----------------------------------------------------------------

df_biv = filter(df_class, Class == "BIVALVIA")

df_invertwq = df_biv %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()

pre_invertwq = df_invertwq %>%
  filter(Date <= median(df_invertwq$Date))

## pick latest sampling each site
df_prem <- pre_invertwq %>% 
  group_by(site_id) %>% 
  slice(which.max(Date)) %>% 
  ungroup()

df_prem = df_prem %>% 
  dplyr::select(c(Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy))

invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

dredge(invertwq_lm)

## post sample
post_invertwq = df_invertwq %>%
  filter(Date > median(df_invertwq$Date))

#prediction model#

wq_lm_mod = lm(Div ~ Diss_Oxy + pH_SU + Sp_Cond, data = df_prem)

d0 = post_invertwq %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond)

wq_pred = predict(wq_lm_mod, newdata = d0) %>% exp()

chisq.test(post_invertwq$Div, wq_pred)

cor.test(post_invertwq$Div, wq_pred, use = "everything")






df_biv = filter(df_class, Class == "BIVALVIA")

df_invertwq = df_biv %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()

pre_invertwq = df_invertwq %>%
  filter(Date <= median(df_invertwq$Date))

## pick latest sampling each site
df_prem <- pre_invertwq %>% 
  group_by(site_id) %>% 
  slice(which.max(Date)) %>% 
  ungroup()

df_prem = df_prem %>% 
  dplyr::select(c(Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy))

invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

dredge(invertwq_lm)

## post sample
post_invertwq = df_invertwq %>%
  filter(Date > median(df_invertwq$Date))

#prediction model#

wq_lm_mod = lm(Div ~ Diss_Oxy + pH_SU + Sp_Cond, data = df_prem)

d0 = post_invertwq %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond)

wq_pred = predict(wq_lm_mod, newdata = d0) %>% exp()

chisq.test(post_invertwq$Div, wq_pred)

cor.test(post_invertwq$Div, wq_pred, use = "everything")
# clitellata --------------------------------------------------------------
df_cli = filter(df_class, Class == "CLITELLATA")

df_invertwq = df_cli %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()

pre_invertwq = df_invertwq %>%
  filter(Date <= median(df_invertwq$Date))

## pick latest sampling each site
df_prem <- pre_invertwq %>% 
  group_by(site_id) %>% 
  slice(which.max(Date)) %>% 
  ungroup()

df_prem = df_prem %>% 
  dplyr::select(c(Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy))

invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

dredge(invertwq_lm)

## post sample
post_invertwq = df_invertwq %>%
  filter(Date > median(df_invertwq$Date))

#prediction model#

dis_lm_mod = lm(Div ~ Diss_Oxy, data = df_prem)

d0 = post_invertwq %>% dplyr::select(Div, Diss_Oxy)

dis_pred = predict(dis_lm_mod, newdata = d0) %>% exp()

chisq.test(post_invertwq$Div, dis_pred)

cor.test(post_invertwq$Div, dis_pred, use = "everything")



