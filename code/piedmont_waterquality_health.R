# library -----------------------------------------------------------------

## remove all objects
remove(list = ls())

pacman::p_load(tidyverse,
               MASS,
               MuMIn)

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
                         "Newberry, SC"))) %>% 
  dplyr::select(c("Date",
                  "Latitude",
                  "Longitude",
                  "Genus",
                  "Family",
                  "Order",
                  "Class",
                  "Abundance"))

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

df_wq_mu = df_wq %>% 
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
  mutate(Latitude = round(Latitude, 4),
         Longitude = round(Longitude, 4),
         Date = as.Date(Date, format = "%m/%d/%y"),
         period = ifelse(Date > as.Date("2000-01-01"), "post", "pre")) %>% 
  group_by(Latitude, Longitude, period) %>% 
  summarize(Temp_C = mean(Temp_C, na.rm = TRUE),
            Sp_Cond = mean(Sp_Cond, na.rm = TRUE),
            pH_SU =  mean(pH_SU, na.rm = TRUE),
            Diss_Oxy =  mean(Diss_Oxy, na.rm = TRUE)) %>% 
  ungroup()


# order data set ----------------------------------------------------------

## pick orders with < 10% unknowns in genus
order_name <- df_table %>% 
  filter(Order != "UNKNOWN") %>% 
  drop_na(Order) %>% 
  pull(Order) %>% 
  unique()

## Order-level data frame
df_order0 = df_invert %>% 
  filter(Order %in% order_name) %>% 
  mutate(Latitude = round(Latitude, 4),
         Longitude = round(Longitude, 4),
         site_id = factor(paste0(Latitude, Longitude)) %>% 
           as.numeric(),
         Date = as.Date(Date, format = "%m/%d/%y"))

## check if NA present in data
site_with_na <- df_order0 %>% 
  group_by(site_id) %>% 
  summarize(n_na_ab = sum(is.na(Abundance)),
            n_na_genus = sum(is.na(Genus))) %>% 
  filter(n_na_ab > 0 | n_na_genus > 0) %>% 
  pull(site_id)

# ## check what are NAs
df_order0 %>% 
  filter(site_id %in% site_with_na) %>% 
  view()

## remove sites with NA abundance data
## group_by site and Order
## append earliest date of sampling as each site has multiple sampling dates
## `Abundance` column had only 1, 3, 10 entries - is this abundance category?
## Latitude/Longitude were rounded to 4 decimal points for left_join
## if data is double, 1.00001 and 1.00000 are considered different!!!!

df_macro <- df_order0 %>% 
  filter(!(site_id %in% site_with_na)) %>% 
  mutate(period = ifelse(Date > as.Date("2000-01-01"), "post", "pre")) %>% 
  group_by(site_id, # by site
           period,
           Latitude,
           Longitude,
           Class,
           Order) %>% 
  summarize(genus_richness = n_distinct(Genus)) %>% 
  ungroup()

## uncomment if you want to check how many dates in each site
# df_macro %>% 
#   group_by(site_id) %>% 
#   summarize(n_date = n_distinct(Date))

## join invert and water quality data
## na.omit() remove rows that contain NA in any column(s)
df_m <- df_macro %>% 
  left_join(df_wq_mu, by = c("Latitude", "Longitude", "period")) %>% 
  na.omit()

# tasks need to be performed ----------------------------------------------

# 1. split data into pre and post
# 2. develop model by Class or Order
# 3. Now models are no longer Normal distribution. Select appropriate distribution

## develop predictions models for EPT

taxa <- c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA")

list_m1 <- lapply(taxa, function(x) {
  
  df_pre <- df_m %>% 
    filter(Order == x,
           period == "pre")
  
  df_glm <- df_pre %>% 
    dplyr::select(c(genus_richness, pH_SU, Sp_Cond, Temp_C, Diss_Oxy))
  
  m <- glm(genus_richness ~ .,
           data = df_glm,
           family = "poisson",
           na.action = "na.fail")
  
  summary(m)
  
  ms <- dredge(m, rank = "AIC")
  
  ## if you want to get only the best
  m1 <- get.models(ms, subset = 1)[[1]]
  
  ## if you want to the full set of models with delta < 2
  # m1 <- get.models(ms, subset = delta < 2)
  
  return(m1)
})

names(list_m1) <- taxa

## predict post genus richness

list_pred <- lapply(seq_len(length(taxa)), function(i) {
  
  df_post <- df_m %>% filter(period == "post",
                             Order == taxa[i])
  
  log_y_hat <- predict(list_m1[[i]],
                       newdata = df_post)
  
  df_out <- df_post %>% 
    mutate(y_hat = exp(log_y_hat))
 
   
  return(df_out)
})
 
# #prediction model#

names(list_pred) <- taxa


# correlation test --------------------------------------------------------

list_rho <- lapply(list_pred,
                   function(subdf) {
                     with(subdf, cor.test(genus_richness,
                                          y_hat,
                                          method = "spearman",
                                          exact = FALSE))
                   })


# prediction vs. observation ----------------------------------------------

list_plot <- lapply(list_pred, function(subdf) {
  
  taxon_name <- stringr::str_to_sentence(unique(subdf$Order))
  
  subdf %>% 
    ggplot(aes(x = y_hat,
               y = genus_richness)) + 
    labs(y = "Predicted",
         x = "Observed") +
    labs(title = paste(taxon_name, # change taxon name based on input
                       "Predicted vs. Observed")) +
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, color = "hotpink")
  
})

## uncomment the following to see each plot
list_plot[[1]]
list_plot[[2]]
list_plot[[3]]

## to save plot PDFs to "output/" sub-directory
lapply(seq_len(length(list_plot)), function(i) {
  filename <- paste0("output/fig_", str_to_lower(taxa[i]), ".pdf")
  ggsave(plot = list_plot[[i]],
         filename = filename,
         width = 5,
         height = 5)
})


# water quality plot ------------------------------------------------------

## example for ephemeloptera
## pick names of selected variables
## [-1] to remove "(Intercept)"

label <- c("Temp_C" = "Temperature (Celsius)",
           "Sp_Cond" = "Specific conductivity",
           "pH_SU" = "pH",
           "Diss_Oxy" = "Dissolved oxygen")

list_wq_plot <- lapply(taxa, function(subx) {
  
  x_name <- names(list_m1[[subx]]$coefficients)[-1]
  
  df_pre <- df_m %>% 
    filter(period == "pre",
           Order == subx) %>% 
    pivot_longer(cols = all_of(x_name),
                 names_to = "x_name",
                 values_to = "value")
  
  g_wq <- df_pre %>% 
    ggplot(aes(x = value,
               y = genus_richness)) +
    facet_wrap(facets = ~x_name,
               scales = "free",
               strip.position = "bottom",
               labeller = as_labeller(label),
               nrow = 2,
               ncol = 2) +
    geom_point() +
    labs(y = "Genus richness") +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle(str_to_sentence(subx))
  
  return(g_wq)
})

names(list_wq_plot) <- taxa
list_wq_plot$EPHEMEROPTERA
list_wq_plot$PLECOPTERA
list_wq_plot$TRICHOPTERA


lapply(seq_len(length(list_wq_plot)), function(i) {
  filename <- paste0("output/fig_", str_to_lower(taxa[i]), "_wq", ".pdf")
  ggsave(plot = list_wq_plot[[i]],
         filename = filename,
         width = 7,
         height = 7)
})

