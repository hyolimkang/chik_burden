library(rjags)
library(ggpubr)
library(scico)

case_comparison <- read_excel("C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/Manuscript/Appendix/case_comparison.xlsx", 
                              sheet = "Sheet2")
# Extract columns
symptomatic <- case_comparison$symp
confirmed <- case_comparison$conf
N <- nrow(case_comparison)

# Model construct
model_code <- "
model {
  for (i in 1:N) {
    confirmed[i] ~ dbinom(confirmed_rate[i], symptomatic[i])
    confirmed_rate[i] ~ dunif(0, 1) # Non-informative prior for each confirmed rate
  }
}
"

# data list
data_list <- list(
  confirmed = confirmed,
  symptomatic = symptomatic,
  N = N
)

# Initial values for the parameters
inits <- function() {
  list(confirmed_rate = confirmed / symptomatic)
}

# Initialize the JAGS model
jags_model <- jags.model(textConnection(model_code), data = data_list, inits = inits, n.chains = 3, n.adapt = 1000)

# Burn-in period
update(jags_model, n.iter = 1000)

# Run the MCMC
mcmc_samples <- coda.samples(model = jags_model, variable.names = "confirmed_rate", n.iter = 5000)

# Summary of MCMC results
summary(mcmc_samples)
gelman.diag(mcmc_samples)

combined_samples <- as.mcmc(do.call(rbind, mcmc_samples))

####
library(scales)
library(readxl)
case_comparison <- read_excel("case_comparison.xlsx")

## remove travel associated countries
case_comparison <- case_comparison %>% 
  filter(remark != "travel")

case_comparison <- case_comparison %>% 
                          mutate(diff = ifelse(symp_med > annual_average, 
                                               symp_med / annual_average, 
                                               annual_average / symp_med))
  

long_case_comparison <- case_comparison %>%
  pivot_longer(cols = c(symp_med, annual_average), names_to = "case_type", values_to = "case_value")

africa_long <- long_case_comparison %>% filter(continent == "Africa")
america_long <- long_case_comparison %>% filter(continent == "Americas")
asia_long <- long_case_comparison %>% filter(continent == "Asia")
europe_long <- long_case_comparison %>% filter(continent == "Europe")
oceania_long <- long_case_comparison %>% filter(continent == "Oceania")

africa_error <- long_case_comparison %>% filter(continent == "Africa")
america_error <- long_case_comparison %>% filter(continent == "Americas")
asia_error <- long_case_comparison %>% filter(continent == "Asia")
europe_error <- long_case_comparison %>% filter(continent == "Europe")
oceania_error <- long_case_comparison %>% filter(continent == "Oceania")

africa_long_filtered <- africa_long %>%
  filter(case_type == "symp_med") %>%
  arrange(case_value)
africa_long$name <- factor(africa_long$name, levels = unique(africa_long_filtered$name))

america_long_filtered <- america_long %>%
  filter(case_type == "symp_med") %>%
  arrange(case_value)
america_long$name <- factor(america_long$name, levels = unique(america_long_filtered$name))

asia_long_filtered <- asia_long %>%
  filter(case_type == "symp_med") %>%
  arrange(case_value)
asia_long$name <- factor(asia_long$name, levels = unique(asia_long_filtered$name))

europe_long_filtered <- europe_long %>%
  filter(case_type == "symp_med") %>%
  arrange(case_value)
europe_long$name <- factor(europe_long$name, levels = unique(europe_long_filtered$name))

oceania_long_filtered <- oceania_long %>%
  filter(case_type == "symp_med") %>%
  arrange(case_value)
oceania_long$name <- factor(oceania_long$name, levels = unique(oceania_long_filtered$name))

################################################################################
afr <- ggplot(data = africa_long, aes(x = name, y = case_value, shape = case_type, color = case_type)) +
  geom_point(size = 3) +
  geom_errorbar(data = africa_error, aes(x = name, ymin = symp_lo, ymax = symp_hi, color = "#0099B499"), width = 0.2, inherit.aes = FALSE) +
  scale_y_log10(labels = comma, breaks = c(100, 1000, 10000, 100000, 1000000)) +
  scale_shape_manual(name = "Case Type", 
                     values = c("annual_average" = 17, "symp_med" = 16),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2011-22 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)")) +
  scale_color_manual(name = "Case Type", 
                     values = c("annual_average" = "#ED000099", "symp_med" = "#0099B499"),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2011-22 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Africa", x = "Country", y = "Symptomatic infections")

ggsave(filename = paste0("Results_figs/afr_valid",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), afr, height=10, width=15, dpi=900)


amr <- ggplot(data = america_long, aes(x = name, y = case_value, shape = case_type, color = case_type)) +
  geom_point(size = 3) +
  geom_errorbar(data = america_error, aes(x = name, ymin = symp_lo, ymax = symp_hi, color = "#0099B499"), width = 0.2, inherit.aes = FALSE) +
  scale_y_log10(labels = comma, breaks = c(100, 1000, 10000, 100000, 1000000)) +
  scale_shape_manual(name = "Case Type", 
                     values = c("annual_average" = 17, "symp_med" = 16),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2013-23 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)")) +
  scale_color_manual(name = "Case Type", 
                     values = c("annual_average" = "#ED000099", "symp_med" = "#0099B499"),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2013-23 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Americas", x = "Country", y = "Symptomatic infections")

ggsave(filename = paste0("Results_figs/amr_valid",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), amr, height=10, width=15, dpi=900)

eur <- ggplot(data = europe_long, aes(x = name, y = case_value, shape = case_type, color = case_type)) +
  geom_point(size = 3) +
  geom_errorbar(data = europe_error, aes(x = name, ymin = symp_lo, ymax = symp_hi, color = "#0099B499"), width = 0.2, inherit.aes = FALSE) +
  scale_y_log10(labels = comma, breaks = c(100, 1000, 10000, 100000, 1000000)) +
  scale_shape_manual(name = "Case Type", 
                     values = c("annual_average" = 17, "symp_med" = 16),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2011-22 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)")) +
  scale_color_manual(name = "Case Type", 
                     values = c("annual_average" = "#ED000099", "symp_med" = "#0099B499"),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2011-22 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Europe", x = "Country", y = "Symptomatic infections")

ggsave(filename = paste0("Results_figs/eur_valid",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), eur, height=10, width=15, dpi=900)


oce <- ggplot(data = oceania_long, aes(x = name, y = case_value, shape = case_type, color = case_type)) +
  geom_point(size = 3) +
  geom_errorbar(data = oceania_error, aes(x = name, ymin = symp_lo, ymax = symp_hi, color = "#0099B499"), width = 0.2, inherit.aes = FALSE) +
  scale_y_log10(labels = comma, breaks = c(100, 1000, 10000, 100000, 1000000)) +
  scale_shape_manual(name = "Case Type", 
                     values = c("annual_average" = 17, "symp_med" = 16),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2011-22 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)")) +
  scale_color_manual(name = "Case Type", 
                     values = c("annual_average" = "#ED000099", "symp_med" = "#0099B499"),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2011-22 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Oceania", x = "Country", y = "Symptomatic infections")

ggsave(filename = paste0("Results_figs/oce_valid",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), oce, height=10, width=15, dpi=900)


asi <- ggplot(data = asia_long, aes(x = name, y = case_value, shape = case_type, color = case_type)) +
  geom_point(size = 3) +
  geom_errorbar(data = asia_error, aes(x = name, ymin = symp_lo, ymax = symp_hi, color = "#0099B499"), width = 0.2, inherit.aes = FALSE) +
  scale_y_log10(labels = comma, breaks = c(100, 1000, 10000, 100000, 1000000)) +
  scale_shape_manual(name = "Case Type", 
                     values = c("annual_average" = 17, "symp_med" = 16),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2011-22 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)")) +
  scale_color_manual(name = "Case Type", 
                     values = c("annual_average" = "#ED000099", "symp_med" = "#0099B499"),
                     labels = c("annual_average" = "Globally reported symptomatic cases (2011-22 average)", 
                                "symp_med" = "Model predicted symptomatic (95%UI median)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Asia", x = "Country", y = "Symptomatic infections")

ggsave(filename = paste0("Results_figs/asi_valid",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), asi, height=10, width=15, dpi=900)

################################################################################
merged_case_comp <- case_comparison %>%
  left_join(admin_boundaries_sf, by = "name")

merged_case_comp <- st_as_sf(merged_case_comp)

comb_burden_sf <- combined_burden %>%
  left_join(admin_boundaries_sf, by = "name")

comb_burden_sf <- st_as_sf(comb_burden_sf)

merged_case_comp <- merged_case_comp %>%
  mutate(annual_average_interval = case_when(
    is.na(annual_average) ~ "NA",
    annual_average > 0 & annual_average <= 100 ~ "0 to 100",
    annual_average > 100 & annual_average <= 1000 ~ "100 to 1,000",
    annual_average > 1000 & annual_average <= 10000 ~ "1,000 to 10,000",
    annual_average > 10000 & annual_average <= 50000 ~ "10,000 to 50,000",
    annual_average > 50000 & annual_average <= 100000 ~ "50,000 to 100,000",
    annual_average > 100000 ~ "> 100,000"
  ))

merged_case_comp <- merged_case_comp %>%
  mutate(symp_interval = case_when(
    is.na(symp_med) ~ "NA",
    symp_med > 0 & symp_med <= 100 ~ "0 to 100",
    symp_med > 100 & symp_med <= 1000 ~ "100 to 1,000",
    symp_med > 1000 & symp_med <= 10000 ~ "1,000 to 10,000",
    symp_med > 10000 & symp_med <= 50000 ~ "10,000 to 50,000",
    symp_med > 50000 & symp_med <= 100000 ~ "50,000 to 100,000",
    symp_med > 100000 ~ "> 100,000"
  ))

colors <- c(
  "NA" = "darkgray",
  "0 to 100" = "#D4EBF2",
  "100 to 1,000" = "#A8C6DB",
  "1,000 to 10,000" = "#7DA1C5",
  "10,000 to 50,000" = "#4D79B3",
  "50,000 to 100,000" = "#3A4B93",
  "> 100,000" = "#2A2A72"
)

world <- map_data("world")

reported_map <- base_map(world) +
   geom_sf(data = merged_case_comp,
            aes(fill = annual_average_interval), alpha = 1)+
  scale_fill_manual(values = colors, na.value = "darkgray",
                    breaks = c("0 to 100", "100 to 1,000", "1,000 to 10,000", 
                               "10,000 to 50,000", "50,000 to 100,000", "> 100,000"),
                    labels = c("0 to 100", "100 to 1,000", "1,000 to 10,000", 
                               "10,000 to 50,000", "50,000 to 100,000", "> 100,000"))+
  coord_sf(xlim = c(-180, 180), ylim = c(-55, 50))+
  labs(fill = "Average globally reported cases",
       caption = "(*American continent: average of 2013-2023, Other countries: average of 2011-2022. Grey areas are countries without globally reported cases at the same period. 
                                         Coloured countries are countries that report locally acquired chikungunya cases(n=102). Four countries (Taiwan, United States of America, Nepal, and China) with no established evidence of local transmission are excluded.)")+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 60, l = 20),  # Adjust bottom margin
    plot.caption = element_text(hjust = 0.3, vjust = 1, size = 8)  # Center align the caption
  )

ggsave(filename = paste0("Results_figs/reported_map",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), reported_map, height=6, width=12, dpi=900)


predicted_map <- base_map(world) +
   geom_sf(data = merged_case_comp,
            aes(fill = symp_interval), alpha = 1)+
  scale_fill_manual(values = colors, na.value = "darkgray",
                    breaks = c("0 to 100", "100 to 1,000", "1,000 to 10,000", 
                               "10,000 to 50,000", "50,000 to 100,000", "> 100,000"),
                    labels = c("0 to 100", "100 to 1,000", "1,000 to 10,000", 
                               "10,000 to 50,000", "50,000 to 100,000", "> 100,000"))+
  coord_sf(xlim = c(-180, 180), ylim = c(-55, 50))+
  labs(fill = "Model predicted symptomatic cases (annual, 95% median)")+
  theme(
    plot.margin = margin(t = 20, r = 20, b = 60, l = 20),  # Adjust bottom margin
    plot.caption = element_text(hjust = 0.3, vjust = 1, size = 8)  # Center align the caption
  )
  
ggsave(filename = paste0("Results_figs/predicted_map",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), predicted_map, height=6, width=12, dpi=900)

combined_plot <- ggarrange(reported_map, predicted_map, ncol = 1, nrow = 2, labels = c("A", "B"))

ggsave(filename = paste0("Results_figs/combined_case_comp",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), combined_plot, height=6, width=12, dpi=900)

### specific country map
brazil_sf <- comb_burden_sf %>% filter(name == "Brazil")

world <- ne_countries(scale = 50, type = "countries", returnclass = "sf")

# Filter for Brazil
brazil <- subset(world, name == "Brazil")

brazil_rast <- foi_sf[foi_sf$name == "Brazil", ]
brazil_sf <- st_as_sf(brazil_rast)
brazil_coords <- st_coordinates(brazil_sf)
brazil_df <- as.data.frame(brazil_sf)

combined_burden_inf <- comb_burden_sf 

# Combine coordinates with data frame
brazil_df <- cbind(brazil_df, brazil_coords)

quantiles <- quantile(df1$layer, probs = c(0, 1), na.rm = TRUE)
palette <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))

# Create the plot
brazil_map <- ggplot(data = brazil)+
  geom_sf(fill = "transparent", color = "grey70")+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title   = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks   = element_blank(),
        panel.grid.major = element_blank())+
  geom_tile(data = brazil_df %>% filter(foi_mid == 0), # df1 for FOI
            aes(x = X, y = Y), fill = "#C8C8C8", alpha = 1) +
  geom_tile(data = brazil_df %>% filter(foi_mid != 0), # df1 for non-zero FOI values
            aes(x = X, y = Y, fill = foi_mid), alpha = 1) +
  scale_fill_scico(palette = "roma")+
  labs(fill = "FOI (95% mid)")

ggsave(filename = paste0("Results_figs/brazil_map",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), brazil_map, height=6, width=12, dpi=900)


### Italy
italy_sf <- comb_burden_sf %>% filter(name == "Italy")

world <- ne_countries(scale = 50, type = "countries", returnclass = "sf")

# Filter for Italy
italy <- subset(world, name == "Italy")
italy_rast <- foi_sf[foi_sf$name == "Italy", ]
italy_sf <- st_as_sf(italy_rast)
italy_coords <- st_coordinates(italy_sf)
italy_df <- as.data.frame(italy_sf)
italy_df <- cbind(italy_df, italy_coords)

# Create the plot
italy_map <- ggplot(data = italy)+
  geom_sf(fill = "transparent", color = "grey70")+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title   = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks   = element_blank(),
        panel.grid.major = element_blank())+
  geom_tile(data = italy_df %>% filter(foi_mid == 0), # df1 for FOI
            aes(x = X, y = Y), fill = "#C8C8C8", alpha = 1) +
  geom_tile(data = italy_df %>% filter(foi_mid != 0), # df1 for non-zero FOI values
            aes(x = X, y = Y, fill = foi_mid), alpha = 1) +
  scale_fill_scico(palette = "roma")+
  labs(fill = "FOI (95% mid)")

ggsave(filename = paste0("Results_figs/italy_map",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), italy_map, height=6, width=12, dpi=900)

## count pixels
counts <- italy_df %>%
  filter(mask == 1) %>%
  count(foi_mid)
tot_pop <- italy_df %>%
            filter(mask == 1) %>%
              summarise(tot = sum(tot))
  
  
combined_med_inf  <- comb_burden_sf %>% select(geometry, med_inf)

italy_sf <- st_join(italy_sf, combined_med_inf)

### france
france <- subset(world, name == "France")
france_rast <- foi_sf[foi_sf$name == "France", ]
france_sf <- st_as_sf(france_rast)
france_coords <- st_coordinates(france_sf)
france_df <- as.data.frame(france_sf)
france_df <- cbind(france_df, france_coords)

# Create the plot
france_map <- ggplot(data = france)+
  geom_sf(fill = "transparent", color = "grey70")+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title   = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks   = element_blank(),
        panel.grid.major = element_blank())+
  geom_tile(data = france_df %>% filter(foi_mid == 0), # df1 for FOI
            aes(x = X, y = Y), fill = "#C8C8C8", alpha = 1) +
  geom_tile(data = france_df %>% filter(foi_mid != 0), # df1 for non-zero FOI values
            aes(x = X, y = Y, fill = foi_mid), alpha = 1) +
  scale_fill_scico(palette = "roma")+
  labs(fill = "FOI (95% mid)")

ggsave(filename = paste0("Results_figs/france_map",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), france_map, height=6, width=12, dpi=900)

## count pixels
counts <- france_df %>%
  filter(mask == 1) 

tot_pop <- france_df %>%
  summarise(tot = sum(tot))

tot_pop_pixel <- france_df %>%
  filter(mask == 1) %>%
  summarise(tot = sum(tot))


combined_med_inf  <- comb_burden_sf %>% select(geometry, med_inf)

france_sf <- st_join(france_sf, combined_med_inf)

### occ map
chik_occ_sf <- st_as_sf(chik_occ, coords = c("Longitude", "Latitude"), crs = 4326)
chik_occ_sf$PA <- 1

occ_map <- base_map(world) +
  geom_sf(data = chik_occ_sf,
          aes(fill = PA), alpha = 0.7, color = "#A8C6DB")

ggsave(filename = paste0("Results_figs/occ_map",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), occ_map, height=6, width=12, dpi=900)

## foi map
foi_comb_sf <- st_as_sf(foi_comb_all_cdc)
foi_comb_rast <- rasterize(foi_comb_sf, tsuit, field = 'foi_mid')
world_coords <- st_coordinates(foi_comb_sf)
world_df <- as.data.frame(foi_comb_sf)
world_df <- cbind(world_df, world_coords)

# Create the plot
base_map(world)+
  geom_sf(fill = "transparent", color = "grey70")+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title   = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks   = element_blank(),
        panel.grid.major = element_blank())+
  geom_tile(data = world_df %>% filter(foi_mid == 0), # df1 for FOI
            aes(x = X, y = Y), fill = "#C8C8C8", alpha = 1) +
  geom_tile(data = world_df %>% filter(foi_mid != 0), # df1 for non-zero FOI values
            aes(x = X, y = Y, fill = foi_mid), alpha = 1) +
  scale_fill_scico(palette = "roma")+
  labs(fill = "FOI (95% mid)")


