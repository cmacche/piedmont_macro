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
df_order1 = df_order0 %>% group_by(site_id,
                                   Date,
                                   Latitude, 
                                   Longitude) %>% 
  summarise(Div = diversity(Abundance,
                            index = "shannon"))

df_invertwq = df_order %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()## pick orderes with < 10% unknowns in genus


df_invertwq1 = df_order1 %>% 
  left_join(df_wq,
            by = c("Date", 
                   "site_id")) %>% 
  drop_na(Latitude.y,
          Longitude.y) %>% 
  ungroup()

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


d0 = df_postm %>% dplyr::select(Diss_Oxy, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Sp_Cond, Temp_C)
d2 = df_postm %>% dplyr::select(Sp_Cond)
d3 = df_postm %>% dplyr::select(Diss_Oxy, Sp_Cond, pH_SU, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) 
wq1_pred = predict(wq1_lm_mod, newdata = d1) 
wq2_pred = predict(wq2_lm_mod, newdata = d2) 
wq3_pred = predict(wq3_lm_mod, newdata = d3) 



chisq.test(df_postm$Div, wq0_pred)

cor.test(df_postm$Div, wq0_pred, use = "everything")


chisq.test(df_postm$Div, wq1_pred)

cor.test(df_postm$Div, wq1_pred, use = "everything")


chisq.test(df_postm$Div, wq2_pred)

cor.test(df_postm$Div, wq2_pred, use = "everything")


chisq.test(df_postm$Div, wq3_pred)

cor.test(df_postm$Div, wq3_pred, use = "everything")










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


d0 = df_postm %>% dplyr::select(Diss_Oxy, pH_SU, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Diss_Oxy, Sp_Cond, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) 
wq1_pred = predict(wq1_lm_mod, newdata = d1) 



cor.test(df_postm$Div, wq0_pred, use = "everything")


cor.test(df_postm$Div, wq1_pred, use = "everything")


plot_data <- data.frame(Predicted_value = wq0_pred,   
                        Observed_value = df_postm$Div) 


ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "hotpink")

plot_data <- data.frame(Predicted_value = wq1_pred,   
                        Observed_value = df_postm$Div) 


ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "hotpink")



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



d0 = df_postm %>% dplyr::select(Diss_Oxy, pH_SU, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Diss_Oxy, pH_SU, Sp_Cond)


wq0_pred = predict(wq0_lm_mod, newdata = d0)
wq1_pred = predict(wq1_lm_mod, newdata = d1)



cor.test(df_postm$Div, wq0_pred, use = "everything")
cor.test(df_postm$Div, wq1_pred, use = "everything")


plot_data <- data.frame(Predicted_value = wq0_pred,   
                        Observed_value = df_postm$Div) 


ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "hotpink")

plot_data <- data.frame(Predicted_value = wq1_pred,   
                        Observed_value = df_postm$Div) 


ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "hotpink")





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


d0 = df_postm %>% dplyr::select(Diss_Oxy, pH_SU, Sp_Cond, Temp_C)


wq0_pred = predict(wq0_lm_mod, newdata = d0) 


cor.test(df_postm$Div, wq0_pred, use = "everything")

plot_data <- data.frame(Predicted_value = wq0_pred,   
                        Observed_value = df_postm$Div) 


ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "hotpink")


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



d0 = df_postm %>% dplyr::select(pH_SU, Sp_Cond, Temp_C)
d1 = df_postm %>% dplyr::select(Diss_Oxy, pH_SU, Sp_Cond, Temp_C)




wq0_pred = predict(wq0_lm_mod, newdata = d0) 
wq1_pred = predict(wq1_lm_mod, newdata = d1) 

cor.test(df_postm$Div, wq0_pred, use = "everything")


cor.test(df_postm$Div, wq1_pred, use = "everything")


plot_data <- data.frame(Predicted_value = wq0_pred,   
                        Observed_value = df_postm$Div) 

# plot predicted values and actual values 
ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "hotpink")

plot_data <- data.frame(Predicted_value = wq1_pred,   
                        Observed_value = df_postm$Div) 

# plot predicted values and actual values 
ggplot(plot_data, aes(x = Predicted_value, y = Observed_value)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "hotpink")


