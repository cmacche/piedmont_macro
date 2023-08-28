# library -----------------------------------------------------------------
library(data.table)
library(MASS)
library(MuMIn)
library(tidyverse)
library(vegan)
library(ggplot2)
library(data.table) 
library(jtools) 
library(raster)
library(scales)
library(gt)

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
                                  Longitude,Class) %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq4_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])
wq5_lm_mod = lm(get.models(invert_dre, subset = 6)[[1]])


d0 = df_postm %>% dplyr::select(Div, Sp_Cond)
d1 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond)
d3 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond)
d4 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond, Temp_C)
d5 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()
wq4_pred = predict(wq4_lm_mod, newdata = d4) %>% exp()
wq5_pred = predict(wq5_lm_mod, newdata = d5) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")


chisq.test(df_postm$Div, wq4_pred)

cor.test(df_postm$Div, wq4_pred, use = "everything")


chisq.test(df_postm$Div, wq5_pred)

cor.test(df_postm$Div, wq5_pred, use = "everything")




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

#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq4_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])

d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond, Temp_C)
d3 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond, Temp_C)
d4 = df_postm %>% dplyr::select(Div, Diss_Oxy, Temp_C, Sp_Cond)

wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()
wq4_pred = predict(wq4_lm_mod, newdata = d4) %>% exp()


chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")


chisq.test(df_postm$Div, wq4_pred)

cor.test(df_postm$Div, wq4_pred, use = "everything")


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

#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(Div ~ Diss_Oxy + Sp_Cond, data = df_prem)
wq2_lm_mod = lm(Div ~ Diss_Oxy + Temp_C, data = df_prem)

d0 = df_postm %>% dplyr::select(Div, Diss_Oxy)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")






# enopla ------------------------------------------------------------------

df_eno = filter(df_class, Class == "ENOPLA")

df_invertwq = df_eno %>% 
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

#no model is good as has 0 diversity at this level#


# gastropoda --------------------------------------------------------------
df_gas = filter(df_class, Class == "GASTROPODA")

df_invertwq = df_gas %>% 
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


invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])


d0 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond, Temp_C)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond, Temp_C)
d3 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")






# hydro -------------------------------------------------------------------

df_hyd = filter(df_class, Class == "HYDROZOA")

df_invertwq = df_hyd %>% 
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

#zero diversity#


# insecta -----------------------------------------------------------------
df_ins = filter(df_class, Class == "INSECTA")

df_invertwq = df_ins %>% 
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

invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq4_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])

d0 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond, Temp_C)
d2 = df_postm %>% dplyr::select(Div, Sp_Cond)
d3 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, Temp_C)
d4 = df_postm %>% dplyr::select(Div, Temp_C)

wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()
wq4_pred = predict(wq4_lm_mod, newdata = d4) %>% exp()


chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")


chisq.test(df_postm$Div, wq4_pred)

cor.test(df_postm$Div, wq4_pred, use = "everything")



# maxillopoda -------------------------------------------------------------

df_max = filter(df_class, Class == "MAXILLOPODA")

df_invertwq = df_max %>% 
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

invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])


d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Temp_C)

wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")





# ophiurioidea ------------------------------------------------------------

df_oph = filter(df_class, Class == "OPHIURIOIDEA")

df_invertwq = df_oph %>% 
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

#zero diversity#



# polychaeta --------------------------------------------------------------

df_pol = filter(df_class, Class == "POLYCHAETA")

df_invertwq = df_pol %>% 
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

#zero diversity#



# turbellaria -------------------------------------------------------------

df_tur = filter(df_class, Class == "TURBELLARIA")

df_invertwq = df_max %>% 
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

invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])



d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Temp_C)

wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")





# class_grid_results ------------------------------------------------------



df_class_results = read_csv("data_raw/class_results.csv")
df_class_results = arrange(df_class_results, cor_result)


df_class_plot = df_class_results %>% gt() %>% 
  tab_header(title = "class diversity prediction results") %>%
  data_color(columns = chi_sq_p_value,
             rows = chi_sq_p_value < 1,
             method = "numeric",
             palette = c("BuPu"),
             domain = c(0, 1)) %>% 
  data_color(columns = cor_p_value,
             rows = cor_p_value < 1,
             method = "numeric",
             palette = c("BuPu"),
             domain = c(0, 1)) %>% 
  data_color(columns = cor_result,
             rows = chi_sq_p_value < 1,
             method = "numeric",
             palette = c("viridis"),
             domain = c(-1, 1)) %>% 
  data_color(columns = c(Diss_Oxy, pH_SU, Sp_Cond, Temp_C),
             method = "auto",
             palette = c("green", "red"),
             domain = c("n","y"))
            
  )
df_melt = melt(df_class_results)
ggplot(df_class_results, aes(chi_sq_p_value, cor_p_value))

ggplot(df_melt, aes(Diss_Oxy, pH_SU, Sp_Cond, Temp_C)) +                          
  geom_tile(aes(fill = value))

df_class_results1 <- df_class_results %>% remove_rownames() %>% group_by(class) %>%
  mutate(newclass = paste0(class, row_number())) %>% 
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("newclass") %>% 
  ungroup() %>% 
  select(-class)

library(plotly)                                                 
install.packages("plotly")
plot_ly(z = df_melt, type = "heatmap")  

install.packages("plotrix")
library(plotrix)

as.matrix(df_class_results)

scale(df_class_results1)


summary(df_class_results)

heatmap(x, scale = "none")
install.packages("RColorBrewer")
library(RColorBrewer)
df <- scale(mtcars)
df = scale(df_class_results1)
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(df, scale = "none", col =  col, 
        RowSideColors = rep(c("blue", "pink"), each = 16),
        ColSideColors = c(rep("purple", 5), rep("orange", 6)))



mat <- as.matrix(expr[, grep("cell", colnames(expr))])
type <- gsub("s\\d+_", "", colnames(mat))
ha = HeatmapAnnotation(
  df = data.frame(type = type),
  annotation_height = unit(4, "mm")
)

mat <- as.matrix(df_class_results1[, grep("cell", colnames(df_class_results1))])
type <- gsub("s\\d+_", "", colnames(mat))
ha = HeatmapAnnotation(
  df = data.frame(type = type),
  annotation_height = unit(4, "mm")
)
install.packages("circlize")
install.packages("ComplexHeatmap")
library(circlize)
library(ComplexHeatmap)

col1 = colorRamp2(c(0,1), c("green", "red"))
col2 = colorRamp2(c(-1, 0, 1), c("lightblue", "blue", "darkblue"))
col3 = colorRamp2(c(0,0.05,1), c("white", "pink","hotpink"))

ht1 = heatmap(df_class_results1[, group == c("Diss_Oxy","pH_SU")], 
              col = col1, name = "Water Quality")
ht2 = Heatmap(mat[, group == "B"], col = col2, name = "Group_B")
ht1 + ht2

# order data set ----------------------------------------------------------

## pick orderes with < 10% unknowns in genus
order_name <- df_table %>% 
  filter(Order != "UNKNOWN") %>% 
  drop_na(Order) %>% 
  pull(Order) %>% 
  unique()

df_order0 = df_invert %>% 
  filter(Order %in% order_name) %>% 
  mutate(site_id = paste0(round(Latitude, 4),
                          round(Longitude, 4)),
         Date = as.Date(Date, format = "%m/%d/%y"))

df_order = df_order0 %>% group_by(site_id,
                                  Date,
                                  Latitude, 
                                  Longitude,
                                  Order) %>% 
  summarise(Div = diversity(Abundance,
                            index = "shannon"))

df_invertwq = df_order %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()## pick orderes with < 10% unknowns in genus
order_name <- df_table %>% 
  filter(Order != "UNKNOWN") %>% 
  drop_na(Order) %>% 
  pull(Order) %>% 
  unique()

df_order0 = df_invert %>% 
  filter(Order %in% order_name) %>% 
  mutate(site_id = paste0(round(Latitude, 4),
                          round(Longitude, 4)),
         Date = as.Date(Date, format = "%m/%d/%y"))

df_order = df_order0 %>% group_by(site_id,
                                  Date,
                                  Latitude, 
                                  Longitude,
                                  Order) %>% 
  summarise(Div = diversity(Abundance,
                            index = "shannon"))

df_invertwq = df_order %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()

# analysis_order ----------------------------------------------------------------

pre_invertwq = df_invertwq %>%
  filter(Date <= median(df_invertwq$Date))

## pick latest sampling each site
df_prem <- pre_invertwq %>% 
  group_by(site_id) %>% 
  slice(which.max(Date)) %>% 
  ungroup()

df_prem = df_prem %>% 
  dplyr::select(c(Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy))


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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])


d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)
d2 = df_postm %>% dplyr::select(Div, Sp_Cond)
d3 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, pH_SU, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")










# trombidiformes ----------------------------------------------------------

df_trom = filter(df_order, Order == "TROMBIDIFORMES")

df_invertwq = df_trom %>% 
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

#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#


# unionida ----------------------------------------------------------------

df_uni = filter(df_order, Order == "UNIONIDA")

df_invertwq = df_uni %>% 
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

#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])


d0 = df_postm %>% dplyr::select(Div, pH_SU)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy)
d2 = df_postm %>% dplyr::select(Div, Sp_Cond)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")




# veneroida ---------------------------------------------------------------


df_ven = filter(df_order, Order == "VENEROIDA")

df_invertwq = df_ven %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq4_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])
wq5_lm_mod = lm(get.models(invert_dre, subset = 6)[[1]])


d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, pH_SU)
d2 = df_postm %>% dplyr::select(Div, Sp_Cond, pH_SU)
d3 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond, Temp_C)
d4 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond, Temp_C)
d5 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()
wq4_pred = predict(wq4_lm_mod, newdata = d4) %>% exp()
wq5_pred = predict(wq5_lm_mod, newdata = d5) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")


chisq.test(df_postm$Div, wq4_pred)

cor.test(df_postm$Div, wq4_pred, use = "everything")


chisq.test(df_postm$Div, wq5_pred)

cor.test(df_postm$Div, wq5_pred, use = "everything")








# arhynchobdellida --------------------------------------------------------

df_arh = filter(df_order, Order == "ARHYNCHOBDELLIDA")

df_invertwq = df_arh %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#0 diversity#


# haplotaxida -------------------------------------------------------------

df_hap = filter(df_order, Order == "HAPLOTAXIDA")

df_invertwq = df_hap %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])



d0 = df_postm %>% dplyr::select(Div, Sp_Cond)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy)



wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()




chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")



# rhynchobdellida ---------------------------------------------------------

df_rhy = filter(df_order, Order == "RHYNCHOBDELLIDA")

df_invertwq = df_rhy %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])



d0 = df_postm %>% dplyr::select(Div, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)



wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()




chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")





# hoplonemertea -----------------------------------------------------------

df_hop = filter(df_order, Order == "HOPLONEMERTEA")

df_invertwq = df_hop %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#


# architaenioglossa -------------------------------------------------------

df_arch = filter(df_order, Order == "ARCHITAENIOGLOSSA")

df_invertwq = df_arch %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#
# basommatophora ----------------------------------------------------------

df_bas = filter(df_order, Order == "BASOMMATOPHORA")

df_invertwq = df_bas %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq4_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])
wq5_lm_mod = lm(get.models(invert_dre, subset = 6)[[1]])
wq6_lm_mod = lm(get.models(invert_dre, subset = 7)[[1]])
wq7_lm_mod = lm(get.models(invert_dre, subset = 8)[[1]])
wq8_lm_mod = lm(get.models(invert_dre, subset = 9)[[1]])
wq9_lm_mod = lm(get.models(invert_dre, subset = 10)[[1]])


d0 = df_postm %>% dplyr::select(Div, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, Temp_C)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy)
d3 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)
d4 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, Temp_C)
d5 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU)
d6 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Temp_C)
d7 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond)
d8 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond, Temp_C)
d9 = df_postm %>% dplyr::select(Div, pH_SU, Temp_C)



wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()
wq4_pred = predict(wq4_lm_mod, newdata = d4) %>% exp()
wq5_pred = predict(wq5_lm_mod, newdata = d5) %>% exp()
wq6_pred = predict(wq6_lm_mod, newdata = d6) %>% exp()
wq7_pred = predict(wq7_lm_mod, newdata = d7) %>% exp()
wq8_pred = predict(wq8_lm_mod, newdata = d8) %>% exp()
wq9_pred = predict(wq9_lm_mod, newdata = d9) %>% exp()




chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")


chisq.test(df_postm$Div, wq4_pred)

cor.test(df_postm$Div, wq4_pred, use = "everything")


chisq.test(df_postm$Div, wq5_pred)

cor.test(df_postm$Div, wq5_pred, use = "everything")


chisq.test(df_postm$Div, wq6_pred)

cor.test(df_postm$Div, wq6_pred, use = "everything")


chisq.test(df_postm$Div, wq7_pred)

cor.test(df_postm$Div, wq7_pred, use = "everything")


chisq.test(df_postm$Div, wq8_pred)

cor.test(df_postm$Div, wq8_pred, use = "everything")


chisq.test(df_postm$Div, wq9_pred)

cor.test(df_postm$Div, wq9_pred, use = "everything")



# neotaenioglossa ---------------------------------------------------------
df_neo = filter(df_order, Order == "NEOTAENIOGLOSSA")

df_invertwq = df_neo %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])



d0 = df_postm %>% dplyr::select(Div, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, Temp_C)




wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()





chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")






# anthoathecata -----------------------------------------------------------

df_ant = filter(df_order, Order == "ANTHOATHECATA")

df_invertwq = df_ant %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#


# coleoptera --------------------------------------------------------------


df_col = filter(df_order, Order == "COLEOPTERA")

df_invertwq = df_col %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])



d0 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Sp_Cond, Diss_Oxy, Temp_C)




wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()





chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")



# diptera -----------------------------------------------------------------


df_dip = filter(df_order, Order == "DIPTERA")

df_invertwq = df_dip %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])



d0 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond, Temp_C)




wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()





chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")



# ephemeroptera -----------------------------------------------------------


df_eph = filter(df_order, Order == "EPHEMEROPTERA")

df_invertwq = df_eph %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])


d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")



# megaloptera -------------------------------------------------------------


df_meg = filter(df_order, Order == "MEGALOPTERA")

df_invertwq = df_meg %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq4_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])
wq5_lm_mod = lm(get.models(invert_dre, subset = 6)[[1]])
wq6_lm_mod = lm(get.models(invert_dre, subset = 7)[[1]])
wq7_lm_mod = lm(get.models(invert_dre, subset = 8)[[1]])



d0 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond)
d3 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond)
d4 = df_postm %>% dplyr::select(Div, Sp_Cond)
d5 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond)
d6 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, Temp_C)
d7 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()
wq4_pred = predict(wq4_lm_mod, newdata = d4) %>% exp()
wq5_pred = predict(wq5_lm_mod, newdata = d5) %>% exp()
wq6_pred = predict(wq6_lm_mod, newdata = d6) %>% exp()
wq7_pred = predict(wq7_lm_mod, newdata = d7) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")


chisq.test(df_postm$Div, wq4_pred)

cor.test(df_postm$Div, wq4_pred, use = "everything")


chisq.test(df_postm$Div, wq5_pred)

cor.test(df_postm$Div, wq5_pred, use = "everything")


chisq.test(df_postm$Div, wq6_pred)

cor.test(df_postm$Div, wq6_pred, use = "everything")


chisq.test(df_postm$Div, wq7_pred)

cor.test(df_postm$Div, wq7_pred, use = "everything")



# neuroptera --------------------------------------------------------------


df_neu = filter(df_order, Order == "NEUROPTERA")

df_invertwq = df_neu %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#


# odonata -----------------------------------------------------------------


df_odo = filter(df_order, Order == "ODONATA")

df_invertwq = df_odo %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq4_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])



d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Temp_C)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, Temp_C)
d3 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)
d4 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Temp_C)




wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()
wq4_pred = predict(wq4_lm_mod, newdata = d4) %>% exp()





chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")


chisq.test(df_postm$Div, wq4_pred)

cor.test(df_postm$Div, wq4_pred, use = "everything")



# plecoptera --------------------------------------------------------------
df_ple = filter(df_order, Order == "PLECOPTERA")

df_invertwq = df_ple %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])



d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


plot_data = data.frame(predicted_value = wq0_pred,
                       observed_value = d0$Div)

# plot predicted values and actual values
ggplot(plot_data, aes(x = predicted_value, y = observed_value)) +
  geom_point() +
  geom_smooth(method = "lm")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")

plot_data = data.frame(predicted_value = wq1_pred,
                       observed_value = d1$Div)

# plot predicted values and actual values
ggplot(plot_data, aes(x = predicted_value, y = observed_value)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "green")


# trichoptera -------------------------------------------------------------

df_tri = filter(df_order, Order == "TRICHOPTERA")

df_invertwq = df_tri %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])


d0 = df_postm %>% dplyr::select(Div,Diss_Oxy, pH_SU, Sp_Cond, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()


chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")

plot_data = data.frame(predicted_value = wq0_pred,
                       observed_value = d0$Div)

# plot predicted values and actual values
ggplot(plot_data, aes(x = predicted_value, y = observed_value)) +
  geom_point() +
  geom_smooth(method = "lm",
                formula = y ~ x)
# amphipoda ---------------------------------------------------------------
df_amp = filter(df_order, Order == "AMPHIPODA")

df_invertwq = df_amp %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])


d0 = df_postm %>% dplyr::select(Div, Sp_Cond)
d1 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond)
d2 = df_postm %>% dplyr::select(Div, Sp_Cond, Temp_C)
d3 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()


chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")



# isopoda -----------------------------------------------------------------

df_iso = filter(df_order, Order == "ISOPODA")

df_invertwq = df_iso %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#


# ophiurida ---------------------------------------------------------------
df_oph = filter(df_order, Order == "OPHIURIDA")

df_invertwq = df_oph %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#


# canalipalpata -----------------------------------------------------------
df_can = filter(df_order, Order == "CANALIPALPATA")

df_invertwq = df_can %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#


# eunicida ----------------------------------------------------------------
df_eun = filter(df_order, Order == "EUNICIDA")

df_invertwq = df_eun %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#


# proseriata --------------------------------------------------------------
df_pros = filter(df_order, Order == "PROSERIATA")

df_invertwq = df_pros %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#


# gordioidea --------------------------------------------------------------
df_gor = filter(df_order, Order == "GORDIOIDEA")

df_invertwq = df_gor %>% 
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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#zero diversity#




# wetland -----------------------------------------------------------------
df_wq = df_wq %>% 
  filter(EcoRegion == "P",
         Sp_Cond != 0,
         pH_SU != 0,
         Diss_Oxy != 0,
         Temp_C != 0,
         Water_Class %in% c("Freshwater emergent wetland", 
                            "Forested Shrub Wetland",
                            "Freshwater Forested/Shrub Wetland")) %>% 
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
                                  Longitude,Class) %>% 
  summarise(Div = diversity(Abundance,
                            index = "shannon"))

df_invertwq = df_class %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()


# analysis_w ----------------------------------------------------------------

pre_invertwq = df_invertwq %>%
  filter(Date <= median(df_invertwq$Date))

## pick latest sampling each site
df_prem <- pre_invertwq %>% 
  group_by(site_id) %>% 
  slice(which.max(Date)) %>% 
  ungroup()

df_prem = df_prem %>% 
  dplyr::select(c(Div, pH_SU, Sp_Cond, Temp_C, Diss_Oxy))


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


#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])


d0 = df_postm %>% dplyr::select(Div, pH_SU, Temp_C)
d1 = df_postm %>% dplyr::select(Div, pH_SU)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU)



wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()




chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")




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

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])


d0 = df_postm %>% dplyr::select(Div, pH_SU, Temp_C)
d1 = df_postm %>% dplyr::select(Div, pH_SU)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU)



wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()




chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")

#not enough diversity observations#


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

#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])
wq4_lm_mod = lm(get.models(invert_dre, subset = 6)[[1]])
wq5_lm_mod = lm(get.models(invert_dre, subset = 7)[[1]])


d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU)
d3 = df_postm %>% dplyr::select(Div, Sp_Cond)
d4 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Sp_Cond)
d5 = df_postm %>% dplyr::select(Div, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()
wq4_pred = predict(wq4_lm_mod, newdata = d4) %>% exp()
wq5_pred = predict(wq5_lm_mod, newdata = d5) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")


chisq.test(df_postm$Div, wq4_pred)

cor.test(df_postm$Div, wq4_pred, use = "everything")


chisq.test(df_postm$Div, wq5_pred)

cor.test(df_postm$Div, wq5_pred, use = "everything")



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

#dredge model#
invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])


d0 = df_postm %>% dplyr::select(Div, Diss_Oxy)
d1 = df_postm %>% dplyr::select(Div, Temp_C)



wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()




chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")









# enopla ------------------------------------------------------------------

df_eno = filter(df_class, Class == "ENOPLA")

df_invertwq = df_eno %>% 
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

#no model is good as has 0 diversity at this level#


# gastropoda --------------------------------------------------------------
df_gas = filter(df_class, Class == "GASTROPODA")

df_invertwq = df_gas %>% 
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


invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])



d0 = df_postm %>% dplyr::select(Div, Diss_Oxy)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU)
d2 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond)



wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()




chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")







# hydro -------------------------------------------------------------------

df_hyd = filter(df_class, Class == "HYDROZOA")

df_invertwq = df_hyd %>% 
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

#zero diversity#


# insecta -----------------------------------------------------------------
df_ins = filter(df_class, Class == "INSECTA")

df_invertwq = df_ins %>% 
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

invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 3)[[1]])
wq2_lm_mod = lm(get.models(invert_dre, subset = 4)[[1]])
wq3_lm_mod = lm(get.models(invert_dre, subset = 5)[[1]])


d0 = df_postm %>% dplyr::select(Div, Sp_Cond)
d1 = df_postm %>% dplyr::select(Div, pH_SU, Sp_Cond, Temp_C)
d2 = df_postm %>% dplyr::select(Div, Sp_Cond)
d3 = df_postm %>% dplyr::select(Div, Diss_Oxy, Sp_Cond, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()
wq2_pred = predict(wq2_lm_mod, newdata = d2) %>% exp()
wq3_pred = predict(wq3_lm_mod, newdata = d3) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")





# maxillopoda -------------------------------------------------------------

df_max = filter(df_class, Class == "MAXILLOPODA")

df_invertwq = df_max %>% 
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

invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])


d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Temp_C)

wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")





# ophiurioidea ------------------------------------------------------------

df_oph = filter(df_class, Class == "OPHIURIOIDEA")

df_invertwq = df_oph %>% 
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

#zero diversity#



# polychaeta --------------------------------------------------------------

df_pol = filter(df_class, Class == "POLYCHAETA")

df_invertwq = df_pol %>% 
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

#zero diversity#



# turbellaria -------------------------------------------------------------

df_tur = filter(df_class, Class == "TURBELLARIA")

df_invertwq = df_max %>% 
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

invertwq_lm = lm(Div~ .,
                 data = df_prem,
                 na.action = "na.fail")

invert_dre = dredge(invertwq_lm)

#prediction model#

wq0_lm_mod = lm(get.models(invert_dre, subset = 1)[[1]])
wq1_lm_mod = lm(get.models(invert_dre, subset = 2)[[1]])



d0 = df_postm %>% dplyr::select(Div, Diss_Oxy, Temp_C)
d1 = df_postm %>% dplyr::select(Div, Diss_Oxy, pH_SU, Temp_C)

wq0_pred = predict(wq0_lm_mod, newdata = d0) %>% exp()
wq1_pred = predict(wq1_lm_mod, newdata = d1) %>% exp()



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


