
#install.packages("MuMIn")

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



df_invert2 = df_invert %>% group_by(Date,
                               Location,
                               County,
                               Latitude, 
                               Longitude,
                               Waterbody) %>% 
  summarise(Tot = sum(Abundance))

df_invert2 = na.omit(df_invert2)
df_invert2 = as.data.frame(unclass(df_invert2), stringsAsFactors = TRUE)

df_wq = filter(df_wq, EcoRegion == "P") %>% 
  filter(!(Sp_Cond == 0)) %>% 
  filter(!(pH_SU == 0)) %>% 
  filter(!(Diss_Oxy == 0)) %>% 
  filter(!(Temp_C == 0))
df_wq = subset(df_wq, select = c( "Date", "County",
                                                "Waterbody","Water_Class", 
                                                "Latitude","Longitude",
                                                "Drainage", 
                                                "Temp_C","Sp_Cond", "pH_SU",
                                                "Diss_Oxy"))


df_wq = as.data.frame(unclass(df_wq),                    
                             stringsAsFactors = TRUE)


df_invertwq = merge(df_invert2, df_wq, by = c("Date", "County",
                                        "Waterbody", 
                                        "Latitude","Longitude"))

df_invertwq$Date <- as.Date(df_invertwq$Date , format = "%m/%d/%y")


pre_invertwq = df_invertwq %>% filter(Date <= median(df_invertwq$Date))

post_invertwq = df_invertwq %>% filter(Date > median(df_invertwq$Date))


pre_invertwq = subset(pre_invertwq, 
                      select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
invertwq_glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail" )

invertwq_dredge = dredge(invertwq_glm)

temp_glm_mod = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot, Temp_C)

temp_pred = predict(temp_glm_mod, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, temp_pred)

cor.test(post_invertwq$Tot, temp_pred, use = "everything")


# habitats ----------------------------------------------------------------


#Lake#

df_lake = filter(df_invertwq, Water_Class == "Lake")

df_lake = na.omit(df_lake)



pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))


pre_lake = subset(pre_lake, 
                      select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail" )

lake_dredge = dredge(lake_glm) # in this case the best model shows
#none of the variables are the best model, going with the second best




ph_glm_mod = glm.nb(Tot ~ pH_SU, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,pH_SU)

ph_pred = predict(ph_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, ph_pred)

cor.test(post_lake$Tot, ph_pred, use = "everything")



#riverine#

df_riv = filter(df_invertwq, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

temp_glm_mod = glm.nb(Tot ~ Temp_C, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C)

temp_pred = predict(temp_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, temp_pred)

cor.test(post_riv$Tot, temp_pred, use = "everything")


#wetland#

df_wet = filter(df_invertwq, Water_Class %in% 
                   c("Freshwater Forested/Shrub Wetland", 
                     "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

spco_glm = glm.nb(Tot ~ Sp_Cond, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, spco_pred)

cor.test(post_invertwq$Tot, spco_pred, use = "everything")





#class#
# Arachnida ---------------------------------------------------------------

#df_arach#
df_arach = filter(df_invert, Class == "ARACHNIDA")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)

pre_arach = df_arach %>% filter(Date <= median(df_arach$Date))

post_arach = df_arach %>% filter(Date > median(df_arach$Date))

pre_arach = subset(pre_arach, 
                      select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
arach_glm = glm.nb(Tot~ ., data = pre_arach, na.action = "na.fail" )

arach_dredge = dredge(arach_glm)



temp_glm_mod = glm.nb(Tot ~ Temp_C, data = pre_arach)

d0 <- post_arach %>% dplyr::select(Tot, Temp_C)

temp_pred = predict(temp_glm_mod, newdata = d0) %>% exp()

chisq.test(post_arach$Tot, temp_pred)

cor.test(post_arach$Tot, temp_pred, use = "everything")


#Lake#

df_lake = filter(df_arach, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



diss_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Diss_Oxy)

diss_pred = predict(diss_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, diss_pred)

cor.test(post_lake$Tot, diss_pred, use = "everything")



#riverine#

df_riv = filter(df_arach, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

temp_glm_mod = glm.nb(Tot ~ Temp_C, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C)

temp_pred = predict(temp_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, temp_pred)

cor.test(post_riv$Tot, temp_pred, use = "everything")


#wetland#

df_wet = filter(df_arach, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

ph_glm_mod = glm.nb(Tot ~ pH_SU, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,pH_SU)

ph_pred = predict(ph_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, ph_pred)

cor.test(post_wet$Tot, ph_pred, use = "everything")

# BIVALVIA ----------------------------------------------------------------

df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                                 Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                         "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                         stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                   select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")
# CLITELLATA --------------------------------------------------------------
df_cli = filter(df_invert, Class == "CLITELLATA")

df_cli = df_cli %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_cli = merge(df_cli, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_cli$Date = as.Date(df_cli$Date , format = "%m/%d/%y")

df_cli = as.data.frame(unclass(df_cli),
                        stringsAsFactors = TRUE)

df_cli$Tot <- as.numeric(df_cli$Tot)


pre_cli = df_cli %>% filter(Date <= median(df_cli$Date))

post_cli = df_cli %>% filter(Date > median(df_cli$Date))

pre_cli = subset(pre_cli, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
cli_glm = glm.nb(Tot~ ., data = pre_cli, na.action = "na.fail" )

cli_dredge = dredge(cli_glm)



wq2_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy, data = pre_cli)

d0 <- post_cli %>% dplyr::select(Tot, Temp_C, Diss_Oxy)

wq2_pred = predict(wq2_glm_mod, newdata = d0) %>% exp()

chisq.test(post_cli$Tot, wq2_pred)

cor.test(post_cli$Tot, wq2_pred, use = "everything")


#Lake#

df_lake = filter(df_cli, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is dissolved oxygen so will try that#



d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, d0_pred)

cor.test(post_lake$Tot, d0_pred, use = "everything")



#riverine#

df_riv = filter(df_cli, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq2_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy)

wq2_pred = predict(wq2_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq2_pred)

cor.test(post_riv$Tot, wq2_pred, use = "everything")


#wetland#

df_wet = filter(df_cli, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)#specific conduct second best#

spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, spco_pred)

cor.test(post_wet$Tot, spco_pred, use = "everything")
# ENOPLA  -----------------------------------------------------------------
df_eno = filter(df_invert, Class == "ENOPLA")

df_eno = df_eno %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_eno = merge(df_eno, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_eno$Date = as.Date(df_eno$Date , format = "%m/%d/%y")

df_eno = as.data.frame(unclass(df_eno),
                        stringsAsFactors = TRUE)

df_eno$Tot <- as.numeric(df_eno$Tot)


pre_eno = df_eno %>% filter(Date <= median(df_eno$Date))

post_eno = df_eno %>% filter(Date > median(df_eno$Date))

pre_eno = subset(pre_eno, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
eno_glm = glm.nb(Tot~ ., data = pre_eno, na.action = "na.fail" )

eno_dredge = dredge(eno_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_eno)

d0 <- post_eno %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_eno$Tot, wq3_pred)

cor.test(post_eno$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_eno, Water_Class == "Lake")



#riverine#

df_riv = filter(df_eno, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, d0_pred)

cor.test(post_riv$Tot, d0_pred, use = "everything")


#wetland#

df_wet = filter(df_eno, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))



# GASTROPODA --------------------------------------------------------------
df_gas = filter(df_invert, Class == "GASTROPODA")

df_gas = df_gas %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_gas = merge(df_gas, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_gas$Date = as.Date(df_gas$Date , format = "%m/%d/%y")

df_gas = as.data.frame(unclass(df_gas),
                        stringsAsFactors = TRUE)

df_gas$Tot <- as.numeric(df_gas$Tot)


pre_gas = df_gas %>% filter(Date <= median(df_gas$Date))

post_gas = df_gas %>% filter(Date > median(df_gas$Date))

pre_gas = subset(pre_gas, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
gas_glm = glm.nb(Tot~ ., data = pre_gas, na.action = "na.fail" )

gas_dredge = dredge(gas_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Sp_Cond + pH_SU, data = pre_gas)

d0 <- post_gas %>% dplyr::select(Tot, Temp_C, Sp_Cond, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_gas$Tot, wq3_pred)

cor.test(post_gas$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_gas, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



#riverine#

df_riv = filter(df_gas, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Sp_Cond + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Sp_Cond, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_gas, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

temp_glm_mod = glm.nb(Tot ~ Temp_C, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot, Temp_C)

temp_pred = predict(temp_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, temp_pred)

cor.test(post_wet$Tot, temp_pred, use = "everything")

# INSECTA -----------------------------------------------------------------
df_inse = filter(df_invert, Class == "INSECTA")

df_inse = df_inse %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_inse = merge(df_inse, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_inse$Date = as.Date(df_inse$Date , format = "%m/%d/%y")

df_inse = as.data.frame(unclass(df_inse),
                        stringsAsFactors = TRUE)

df_inse$Tot <- as.numeric(df_inse$Tot)


pre_inse = df_inse %>% filter(Date <= median(df_inse$Date))

post_inse = df_inse %>% filter(Date > median(df_inse$Date))

pre_inse = subset(pre_inse, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
inse_glm = glm.nb(Tot~ ., data = pre_inse, na.action = "na.fail" )

inse_dredge = dredge(inse_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + Sp_Cond, data = pre_inse)

d0 <- post_inse %>% dplyr::select(Tot, Temp_C, Diss_Oxy, Sp_Cond)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_inse$Tot, wq3_pred)

cor.test(post_inse$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_inse, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



ph_glm_mod = glm.nb(Tot ~ pH_SU, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,pH_SU)

ph_pred = predict(ph_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, ph_pred)

cor.test(post_lake$Tot, ph_pred, use = "everything")



#riverine#

df_riv = filter(df_inse, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + Sp_Cond, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, Sp_Cond)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_inse, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, spco_pred)

cor.test(post_wet$Tot, spco_pred, use = "everything")



# MAXILLOPODA -------------------------------------------------------------
df_maxi = filter(df_invert, Class == "MAXILLOPODA")

df_maxi = df_maxi %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_maxi = merge(df_maxi, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_maxi$Date = as.Date(df_maxi$Date , format = "%m/%d/%y")

df_maxi = as.data.frame(unclass(df_maxi),
                        stringsAsFactors = TRUE)

df_maxi$Tot <- as.numeric(df_maxi$Tot)


pre_maxi = df_maxi %>% filter(Date <= median(df_maxi$Date))

post_maxi = df_maxi %>% filter(Date > median(df_maxi$Date))

pre_maxi = subset(pre_maxi, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
maxi_glm = glm.nb(Tot~ ., data = pre_maxi, na.action = "na.fail" )

maxi_dredge = dredge(maxi_glm)



wq2_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy, data = pre_maxi)

d0 <- post_maxi %>% dplyr::select(Tot, Temp_C, Diss_Oxy)

wq2_pred = predict(wq2_glm_mod, newdata = d0) %>% exp()

chisq.test(post_maxi$Tot, wq2_pred)

cor.test(post_maxi$Tot, wq2_pred, use = "everything")



#Lake#

df_lake = filter(df_maxi, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



temp_glm_mod = glm.nb(Tot ~ Temp_C, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Temp_C)

temp_pred = predict(temp_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, temp_pred)

cor.test(post_lake$Tot, temp_pred, use = "everything")



#riverine#

df_riv = filter(df_maxi, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + Sp_Cond, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, Sp_Cond)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_maxi, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

wq2_glm_mod = glm.nb(Tot ~ Diss_Oxy + pH_SU, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy,pH_SU)

wq2_pred = predict(wq2_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, wq2_pred)

cor.test(post_wet$Tot, wq2_pred, use = "everything")

# TURBELLARIA -------------------------------------------------------------

df_turb = filter(df_invert, Class == "TURBELLARIA")

df_turb = df_turb %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_turb = merge(df_turb, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_turb$Date = as.Date(df_turb$Date , format = "%m/%d/%y")

df_turb = as.data.frame(unclass(df_turb),
                        stringsAsFactors = TRUE)

df_turb$Tot <- as.numeric(df_turb$Tot)


pre_turb = df_turb %>% filter(Date <= median(df_turb$Date))

post_turb = df_turb %>% filter(Date > median(df_turb$Date))

pre_turb = subset(pre_turb, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
turb_glm = glm.nb(Tot~ ., data = pre_turb, na.action = "na.fail" )

turb_dredge = dredge(turb_glm)



wq4_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU + Sp_Cond, data = pre_turb)

d0 <- post_turb %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU, Sp_Cond)

wq4_pred = predict(wq4_glm_mod, newdata = d0) %>% exp()

chisq.test(post_turb$Tot, wq4_pred)

cor.test(post_turb$Tot, wq4_pred, use = "everything")


#Lake#

df_lake = filter(df_turb, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))



#riverine#

df_riv = filter(df_turb, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq4_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU + Sp_Cond, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU, Sp_Cond)

wq4_pred = predict(wq4_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq4_pred)

cor.test(post_riv$Tot, wq4_pred, use = "everything")


#wetland#

df_wet = filter(df_turb, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")





#order#
# AMPHIPODA ---------------------------------------------------------------
df_amph = filter(df_invert, Order == "AMPHIPODA")

df_amph = df_amph %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_amph = merge(df_amph, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_amph$Date = as.Date(df_amph$Date , format = "%m/%d/%y")

df_amph = as.data.frame(unclass(df_amph),
                        stringsAsFactors = TRUE)

df_amph$Tot <- as.numeric(df_amph$Tot)


pre_amph = df_amph %>% filter(Date <= median(df_amph$Date))

post_amph = df_amph %>% filter(Date > median(df_amph$Date))

pre_amph = subset(pre_amph, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
amph_glm = glm.nb(Tot~ ., data = pre_amph, na.action = "na.fail" )

amph_dredge = dredge(amph_glm)



ph_glm_mod = glm.nb(Tot ~ Temp_C + pH_SU, data = pre_amph)

d0 <- post_amph %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_amph$Tot, wq3_pred)

cor.test(post_amph$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_amph, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_amph, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_amph, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# ARCHITAENIOGLOSSA -------------------------------------------------------

df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")
# ARHYNCHOBDELLIDA --------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# BASOMMATOPHORA ----------------------------------------------------------

df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")
# BRANCHIOBDELLIDA --------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# COLEOPTERA --------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# DECAPODA ----------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# DIPTERA -----------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# EPHEMEROPTERA -----------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# HAPLOTAXIDA -------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# HEMIPTERA ---------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                                 Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                         "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                         stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                   select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")
# HOPLONEMERTEA -----------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# ISOPODA -----------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# LEPIDOPTERA -------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# LUMBRICULIDA ------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# MEGALOPTERA -------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# NEOTAENIOGLOSSA ---------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# NEUROPTERA --------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# ODONATA -----------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# OPISTHOPORA -------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")
# OTHER_TAXA --------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# PLECOPTERA --------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# RHYNCHOBDELLIDA ---------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")
# TRICHOPTERA -------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# TRICLADIDA --------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# TROMBIDIFORMES ----------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# UNIONIDA ----------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")

# VENEROIDA ---------------------------------------------------------------
df_biva = filter(df_invert, Class == "BIVALVIA")

df_biva = df_biva %>% group_by(Date, Location, County, Waterbody,Latitude,
                               Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_biva = merge(df_biva, df_wq, by = c("Date", "County",
                                       "Waterbody","Latitude", "Longitude" ))


df_biva$Date = as.Date(df_biva$Date , format = "%m/%d/%y")

df_biva = as.data.frame(unclass(df_biva),
                        stringsAsFactors = TRUE)

df_biva$Tot <- as.numeric(df_biva$Tot)


pre_biva = df_biva %>% filter(Date <= median(df_biva$Date))

post_biva = df_biva %>% filter(Date > median(df_biva$Date))

pre_biva = subset(pre_biva, 
                  select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
biva_glm = glm.nb(Tot~ ., data = pre_biva, na.action = "na.fail" )

biva_dredge = dredge(biva_glm)



wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_biva)

d0 <- post_biva %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_biva$Tot, wq3_pred)

cor.test(post_biva$Tot, wq3_pred, use = "everything")


#Lake#

df_lake = filter(df_biva, Water_Class == "Lake")

df_lake = na.omit(df_lake)


pre_lake = df_lake %>% filter(Date <= median(df_lake$Date))

post_lake = df_lake %>% filter(Date > median(df_lake$Date))

pre_lake= subset(pre_lake, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
lake_glm = glm.nb(Tot~ ., data = pre_lake, na.action = "na.fail") 

lake_dredge = dredge(lake_glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



spco_glm_mod = glm.nb(Tot ~ Sp_Cond, data = pre_lake)

d0 <- post_lake %>% dplyr::select(Tot,Sp_Cond)

spco_pred = predict(spco_glm_mod, newdata = d0) %>% exp()

chisq.test(post_lake$Tot, spco_pred)

cor.test(post_lake$Tot, spco_pred, use = "everything")



#riverine#

df_riv = filter(df_biva, Water_Class == "Riverine")

df_riv = na.omit(df_riv)

pre_riv =  df_riv %>% filter(Date <= median(df_riv$Date))

post_riv = df_riv %>% filter(Date > median(df_riv$Date))


pre_riv= subset(pre_riv, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
riv_glm = glm.nb(Tot~ ., data = pre_riv, na.action = "na.fail") 

riv_dredge = dredge(riv_glm)

wq3_glm_mod = glm.nb(Tot ~ Temp_C + Diss_Oxy + pH_SU, data = pre_riv)

d0 <- post_riv %>% dplyr::select(Tot, Temp_C, Diss_Oxy, pH_SU)

wq3_pred = predict(wq3_glm_mod, newdata = d0) %>% exp()

chisq.test(post_riv$Tot, wq3_pred)

cor.test(post_riv$Tot, wq3_pred, use = "everything")


#wetland#

df_wet = filter(df_biva, Water_Class %in% 
                  c("Freshwater Forested/Shrub Wetland", 
                    "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)


pre_wet =  df_wet %>% filter(Date <= median(df_wet$Date))

post_wet = df_wet %>% filter(Date > median(df_wet$Date))


pre_wet= subset(pre_wet, select = c("Tot", "pH_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wet_glm = glm.nb(Tot~ ., data = pre_wet, na.action = "na.fail") 

wet_dredge = dredge(wet_glm)

d0_glm_mod = glm.nb(Tot ~ Diss_Oxy, data = pre_wet)

d0 <- post_wet %>% dplyr::select(Tot,Diss_Oxy)

d0_pred = predict(d0_glm_mod, newdata = d0) %>% exp()

chisq.test(post_wet$Tot, d0_pred)

cor.test(post_wet$Tot, d0_pred, use = "everything")


#Family#

# Aeshnidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Ameletidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Ancylidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Asellidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Athericidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Baetidae ----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Baetiscidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Belostomatidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")


# Brachycentridae ---------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Branchiobdellidae -------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Caenidae ----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Calamoceratidae ---------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Calopterygidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Cambaridae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Capniidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Ceratopogonidae ---------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Chaoboridae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Chironomidae ------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Chloroperlidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Coenagrionidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Corbiculidae ------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Cordelegastridae  -------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Cordulegastridae --------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Corduliidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Corixidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Corydalidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Crambidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Crangonyctidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Culicidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Dipseudopsidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Dixidae -----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")


# Dryopidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Dugesiidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Dytiscidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Elmidae -----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Empididae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Enchytraeidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Ephemerellidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Ephemeridae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Erpobdellidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Gammaridae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Gerridae ----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Glossiphoniidae ---------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Glossosomatidae ---------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Gomphidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Gyrinidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Haliplidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Hyalellidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Hydrachnidae ------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Hydrobiidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Hydrochidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Hydrometridae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Hydrophilidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Hydropsychidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Hydroptilidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Isonychiidae ------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Lepidostomatidae --------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Leptoceridae ------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Leptohyphidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Leptophlebiidae ---------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Leuctridae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Libellulidae ------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Limnephilidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Limoniidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Lumbriculidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Lymnaeidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Metretopodidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Molannidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Naididae ----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Nemouridae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Neoephemeridae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Nepidae -----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Noteridae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Nymphalidae/Limoniidae --------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Odontoceridae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Palaemonidae ------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Pediciidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Peltoperlidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Perlidae ----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Perlodidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Philopotamidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Phryganeidae ------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Physidae ----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Pisidiidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Planariidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Planorbidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Pleuroceridae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Polycentropodidae -------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Polymitarcyidae ---------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Potamanthidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Psephenidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Psychodidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Psychomyiidae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Pteronarcyidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Ptilodactylidae ---------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Ptilodactylidae ---------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Pyralidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Rhyacophilidae ----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Scirtidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Sialidae ----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Simuliidae --------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Siphlonuridae -----------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Sisyridae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Tabanidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Taeniopterygidae --------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Tanyderidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Tetrastemmatidae --------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Tipulidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Tubificidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Uenoidae ----------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Unionidae ---------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

# Viviparidae -------------------------------------------------------------
df_arach = filter(df_invert, Class == "df_arach")

df_arach = df_arach %>% group_by(Date, Location, County, Waterbody,Latitude,
                                   Longitude) %>% 
  summarise(Tot = sum(Abundance))
df_arach = merge(df_arach, df_wq, by = c("Date", "County",
                                                  "Waterbody","Latitude", "Longitude" ))


df_arach$Date = as.Date(df_arach$Date , format = "%m/%d/%y")

df_arach = as.data.frame(unclass(df_arach),
                          stringsAsFactors = TRUE)

df_arach$Tot <- as.numeric(df_arach$Tot)

summary(df_arach)



#Before#
pre_invertwq = df_arach %>% filter(Date <= '2005-11-24')
#After#
post_invertwq = df_arach %>% filter(Date > '2005-11-24')


pre_invertwq1 = subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
wq.glm = glm.nb(Tot~ ., data = pre_invertwq1, na.action = "na.fail" )
wq = dredge(wq.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")

#df_lake#
df_lake = filter(df_arach, Water_Class == "df_lake")

df_lake = na.omit(df_lake)

summary(df_lake)

pre_invertwq = df_lake %>% filter(Date <= '2006-08-23')

post_invertwq = df_lake %>% filter(Date > '2006-08-23')

pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
l.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

l = dredge(l.glm)#this is showing none of the waterquality would work but the
#next best is ph so will try that#



DO = glm.nb(Tot ~ Diss_Oxy, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,Diss_Oxy)

DO_pred = predict(DO, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, DO_pred)

cor.test(post_invertwq$Tot, DO_pred, use = "everything")


#df_riv#
df_riv = filter(df_arach, Water_Class == "df_riv")

df_riv = na.omit(df_riv)

summary(df_riv)

pre_invertwq =  df_riv %>% filter(Date <= '2006-04-04')

post_invertwq = df_riv %>% filter(Date > '2006-04-04')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
r.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

r = dredge(r.glm)

Temperature = glm.nb(Tot ~ Temp_C, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot:Temp_C)

Temp_pred = predict(Temperature, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, Temp_pred)

cor.test(post_invertwq$Tot, Temp_pred, use = "everything")


#df_wet#

df_wet = filter(df_arach, Water_Class %in% 
                   c("Freshwater Forested/Shrub df_wet", "Freshwater emergent wetland"))

df_wet = na.omit(df_wet)

summary(df_wet)

pre_invertwq =  df_wet %>% filter(Date <= '2005-08-25')

post_invertwq = df_wet %>% filter(Date > '2005-08-25')


pre_invertwq= subset(pre_invertwq, select = c("Tot", "ph_SU", "Sp_Cond","Temp_C","Diss_Oxy"))
w.glm = glm.nb(Tot~ ., data = pre_invertwq, na.action = "na.fail") 

w = dredge(w.glm) #showed none of the water qualtiy would be best but second best
# specfic conductivity#

ph = glm.nb(Tot ~ ph_SU, data = pre_invertwq)

d0 <- post_invertwq %>% dplyr::select(Tot,ph_SU)

ph_pred = predict(ph, newdata = d0) %>% exp()

chisq.test(post_invertwq$Tot, ph_pred)

cor.test(post_invertwq$Tot, ph_pred, use = "everything")

