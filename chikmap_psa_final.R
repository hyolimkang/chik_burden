load("ensemble_list_hyper_raw_0415.RData")
load("ensemble_list_hyper_raw_0625.RData")
load("foi_mat.RData")
load("foi_comb_all_cdc.RData")
load("foi_comb_all_cdc_0623.RData")
load("foi_comb_all_df.RData")
load("combined_burden_0427.RData")
load("combined_burden_0623.RData")
load("combined_burden_recent5.RData")
load("allfoi.RData")
load("foi_mat.RData")
load("burden_1.RData")
load("burden_2.RData")
load("burden_3.RData")
load("burden_4.RData")
source("CHIK_mapping/Functions/BurdenFunctions_v2.R")
source("library.R")
# 1. preprare random FOIs per pixel (100 FOIs)
mapRasCon <- function(pred_vec) { 
  output <- tsuit
  values(output) = as.numeric(pred_vec)
  return(output)
}

ensemble_raster <- list()
for(i in 1:length(ensemble_list_hyper_raw)){
  mapresult          <- mapRasCon(ensemble_list_hyper_raw[[i]])
  ensemble_raster[[i]] <- mapresult 
}

foi_dfs <- list()
for (i in seq_along(ensemble_raster)) {
  # Convert the current raster to points, retaining spatial information
  points_spatial <- rasterToPoints(ensemble_raster[[i]], spatial = TRUE)
  
  # Add the resulting SpatialPointsDataFrame to the list
  foi_dfs[[i]] <- points_spatial
}
foi_dat <- lapply(foi_dfs, function(x) data.frame(x@data))
foi_mat <- do.call(cbind, foi_dat)
foi_mat$long <- coordinates(foi_dfs[[1]])[, 1]
foi_mat$lat <- coordinates(foi_dfs[[1]])[, 2]
new_names <- paste0("foi", 1:100)
colnames(foi_mat)[1:100] <- new_names
colnames(foi_mat)[101:102] <- c("x","y")

save(foi_mat, file = "foi_mat.RData")
save(foi_mat, file = "foi_mat_0623.RData") # without CHIK Risk
### FOI LHS
n_samp <- nrow(foi_mat)
n_foi_cols <- ncol(foi_mat) - 2  # Assuming the first two columns are long and lat
lhs_fois <- matrix(NA, nrow = n_samp, ncol = 1)

sampled_fois <- apply(lhs_fois, 1, function(row) {
  sample(foi_mat[row, -c(101,102)], 1)  # Sample one FOI value for each pixel
})


### input param LHS
library(lhs)
set.seed(123)
runs = 1000
A <- randomLHS (n = runs, 
                k = 16) 

symp <- c(0.741, 0.364, 0.997) # mid, 95%UI
fatal <- c(0.00024) # mid
hosp <- c(0.04, 0.03, 0.05) # mid, 95%UI 
lt  <- c(0.51, 0.45, 0.58)
le_lost <- 52
dw_hosp <- c(0.81, 0.6, 0.92)
dw_nonhosp <- c(0.71, 0.49, 0.9)
dw_chronic <- c(0.349, 0.117, 0.581)
nonhosp_dur <- c(6/365, 2/365, 21/365)
hosp_dur <- c(7/365, 3/365, 15/365) 
chronic_dur <- c(1.717, 1, 4.416)

# Calculate mu (Log Mean)
(mu_nonhospdur <- log(6/365))

# Calculate sigma (Log SD)
(sigma_nonhospdur <- (log(21/365) - log(2/365)) / (2 * 1.96))

beta <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
med_symp_overall <- 1-0.476
ui_symp_overall  <- c(1-0.359, 1-0.594)
var_symp_overall <- ((ui_symp_overall[2] - ui_symp_overall[1])/3.92)^2
beta(med_symp_overall, var_symp_overall)

med_symp_asia <- 1-0.41
ui_symp_asia  <- c(1-0.306, 1-0.516)
var_symp_asia <- ((ui_symp_asia[2] - ui_symp_asia[1])/3.92)^2
beta(med_symp_asia, var_symp_asia)

med_symp_africa <- 1-0.483
ui_symp_africa  <- c(1-0.364, 1-0.603)
var_symp_africa <- ((ui_symp_africa[2] - ui_symp_africa[1])/3.92)^2
beta(med_symp_africa, var_symp_africa)

med_symp_america <- 1-0.471
ui_symp_america  <- c(1-0.355, 1-0.588)
var_symp_america <- ((ui_symp_america[2] - ui_symp_america[1])/3.92)^2
beta(med_symp_america, var_symp_america)

med_fatal_nonhosp <- 0.0010947685
ui_fatal_nonhosp  <- c(0.0010171846, 0.0011782034)
var_fatal_nonhosp <- ((ui_fatal_nonhosp[2] - ui_fatal_nonhosp[1])/3.92)^2
beta(med_fatal_nonhosp, var_fatal_nonhosp)

med_dwhosp <- 0.81
ui_dwhosp  <- c(0.6, 0.92)
var_dwhosp <- ((ui_dwhosp[2] - ui_dwhosp[1])/3.92)^2
beta(med_dwhosp, var_dwhosp)

# Given data
fatalities <- c(39, 74, 12, 2)
total_cases <- c(410, 1821, 450, 171)

# Calculate successes and failures
successes <- sum(fatalities)
failures <- sum(total_cases - fatalities)

# Calculate alpha and beta parameters for the Beta distribution
a<- successes + 1
b <- failures + 1

# Propagate uncertainty using the Beta distribution
# For example, generating a sample of fatality rates
set.seed(123)  # For reproducibility
fatality_samples <- rbeta(10000, a, b)

# Summary statistics
(mean_fatality_rate <- mean(fatality_samples))
(ci_95 <- quantile(fatality_samples, c(0.025, 0.975)))

med_fatal_hosp <- 0.0010947685
ui_fatal_hosp  <- c(0.0010171846 , 0.0011782034)
var_fatal_hosp <- ((ui_fatal_hosp[2] - ui_fatal_hosp[1])/3.92)^2
beta(med_fatal_hosp, var_fatal_hosp)

med_hosp <- 30076/1672214.4992
ui_hosp  <- c(30076/1101748.3173, 30076/2370925.1626)
var_hosp <- ((ui_hosp[2] - ui_hosp[1])/3.92)^2
beta(med_hosp, var_hosp)

##
med_lt <- 0.51
ui_lt  <- c(0.45, 0.58)
var_lt <- ((ui_lt[2] - ui_lt[1])/3.92)^2
beta(med_lt, var_lt)

med_le <- 52
ui_le  <- c(45, 60)
var_le <- ((ui_le[2] - ui_le[1])/3.92)^2
beta(med_lt, var_le)

med_dwc <- 0.233
ui_dwc  <- c(0.117, 0.581)
var_dwc <- ((ui_dwc[2] - ui_dwc[1])/3.92)^2
beta(med_dwc, var_dwc)

med_cdur <- 1.676667
ui_cdur  <- c(0.25, 2)
var_cdur <- ((ui_cdur[2] - ui_cdur[1])/3.92)^2
((var_cdur)^(1/2))
(log(med_cdur))
(sigma = (log(2) - log(0.25)) / (2 * 1.96))

med_dwhosp <- 0.133
ui_dwhosp <- c(0.088, 0.19)
var_dwhosp <- ((ui_dwhosp[2] - ui_dwhosp[1])/3.92)^2
beta(med_dwhosp, var_dwhosp)

med_dwnonhosp <- 0.092
ui_dwnonhosp <- c(0.031, 0.133)
var_dwnonhosp <- ((ui_dwnonhosp[2] - ui_dwnonhosp[1])/3.92)^2
beta(med_dwnonhosp, var_dwnonhosp) 

med_ppv <- 0.683
ui_ppv <- c(0.654, 0.711)
var_ppv <- ((ui_ppv[2] - ui_ppv[1])/3.92)^2
beta(med_ppv, var_ppv) 

med_le <- 82-62
ui_le <- c(82-44, 82-74)
log(med_le)
(log(ui_le[1]) - log(ui_le[2]))/4

### 
lhs_sample <- matrix (nrow = nrow(A), ncol = ncol(A))

lhs_sample [,1]   <- qbeta (p = A[,1], shape1 = 49.14034, shape2 = 34.14837, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,2]   <- qbeta (p = A[,2], shape1 = 34.21298, shape2 = 31.963, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,3]   <- qbeta (p = A[,3], shape1 = 36.77819, shape2 = 32.7458, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,4]   <- qbeta (p = A[,4], shape1 = 35.84287, shape2 = 32.55955, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,5]   <- qbeta (p = A[,5], shape1 = 593.7136, shape2 = 238748.2, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,6]   <- qbeta (p = A[,6], shape1 = 709.5568, shape2 = 647424.5, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,7]   <- qbeta (p = A[,7], shape1 = 58.96698, shape2 = 1415.207, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,8]   <- qbeta (p = A[,8], shape1 = 115.3736, shape2 = 110.8491, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,9]   <- qlnorm (p = A[,9], meanlog = 2.995732, sdlog =  0.3895362, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,10]   <- qbeta (p = A[,10], shape1 = 2.738963, shape2 = 9.016243, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,11]  <- qgamma(p = A[,11], shape = ((1.676667/0.4464286)^2), rate = (1.676667/(0.4464286)^2), lower.tail = TRUE,
                            log.p = FALSE)
lhs_sample [,12]  <- qbeta (p = A[,12], shape1 = 22.51835, shape2 = 146.7926, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,13]  <- qlnorm (p = A[,13], meanlog =  -3.953987, sdlog =  0.4105709, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,14]  <- qbeta (p = A[,14], shape1 = 11.25898, shape2 = 111.1212, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,15]  <- qlnorm (p = A[,15], meanlog = -4.108138, sdlog =  0.5998406, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,16]  <- qbeta (p = A[,16], shape1 = 698.7129, shape2 = 324.2928, ncp = 0, lower.tail = TRUE, log.p = FALSE)


#check 
quantile(lhs_sample[,1], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,2], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,3], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,4], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,5], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,6], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,7], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,8], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,9], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,10], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,11], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,12], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,13], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,14], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,15], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,16], c(0.025, 0.5, 0.975))

cols <- c("symp_asia", "symp_africa", "symp_america", "symp_overall", "fatal_hosp", "fatal_nonhosp", "hosp", 
          "lt", "le_lost", "dw_chronic", "dur_chronic",
          "dw_hosp", "dur_hosp", "dw_nonhosp", "dur_nonhosp", "ppv")

colnames (lhs_sample) <- cols

## lhs for country-specific medical care
infection_by_iso3 <- read_excel("infection_by_iso3.xlsx")

B <- randomLHS (n = 1000, 
                k = nrow(infection_by_iso3)) 

lhs_sample_medcare <- matrix (nrow = nrow(B), ncol = ncol(B))

calculate_beta_params <- function(mean, ci_low, ci_high) {
  var <- ((ci_high - ci_low) / (2 * 1.96))^2
  alpha <- ((1 - mean) / var - 1 / mean) * mean^2
  beta <- alpha * (1 / mean - 1)
  return(list(alpha = alpha, beta = beta))
}

generate_lhs_medcare <- function(infection_df, lhs_sample) {
  num_samples <- 1000
  medcare_samples <- matrix(NA, nrow = num_samples, ncol = nrow(infection_df))
  
  for (i in 1:nrow(infection_df)) {
    medcare_mean <- infection_df$medcare_mid[i]
    ci_low <- infection_df$medcare_lo[i]
    ci_high <- infection_df$medcare_hi[i]
    
    beta_params <- calculate_beta_params(medcare_mean, ci_low, ci_high)
    
    medcare_samples[,i] <- qbeta(B[, i], shape1 = beta_params$alpha, shape2 = beta_params$beta)
  }
  
  return(medcare_samples)
}

med_care <- generate_lhs_medcare(infection_by_iso3, lhs_sample_medcare)


## merge all FOI data with population 
chik_binary <- raster("CHIK_binmap_2024_04_24.tif") # Risk map
chik_occ <- raster("CHIK_binmap_by_country2024_04_24.tif") # occurrence map (occ data to country level)

foi_comb_f <- merge(foi_mat, f_pop, by = c("x", "y")) # total female population combined with foi data
foi_comb_m <- merge(foi_mat, m_pop, by = c("x", "y")) # total female population combined with foi data

mask_values <- extract(chik_binary, foi_comb_f[, c("x", "y")])
mask_val_occ <- extract(chik_occ, foi_comb_f[, c("x", "y")])

foi_comb_f <- foi_comb_f %>%
  mutate(mask_val_bin = mask_values, mask_val_occ = mask_val_occ)
foi_comb_m <- foi_comb_m %>%
  mutate(mask_val_bin = mask_values, mask_val_occ = mask_val_occ)

foi_comb_f_mask <- foi_comb_f[!is.na(foi_comb_f$"f 0"), ]
foi_comb_m_mask <- foi_comb_m[!is.na(foi_comb_m$"m 0"), ]

## 2. merge with admin boundaries
admin_boundaries <- st_read("World_shape/world-administrative-boundaries.shp")
admin_boundaries_sf <- st_as_sf(admin_boundaries, wkt = "geometry", crs = 4326)

foi_comb_f_mask <- st_as_sf(foi_comb_f_mask, coords = c("x", "y"), crs = 4326)
foi_comb_m_mask <- st_as_sf(foi_comb_m_mask, coords = c("x", "y"), crs = 4326)

foi_comb_f_mask <- st_join(foi_comb_f_mask, admin_boundaries_sf, join = st_nearest_feature)
foi_comb_m_mask <- st_join(foi_comb_m_mask, admin_boundaries_sf, join = st_nearest_feature)

foi_comb_f_df <- as.data.frame(foi_comb_f_mask)
foi_comb_m_df <- as.data.frame(foi_comb_m_mask)

save("foi_comb_f_df", "foi_comb_m_df", file = "foi_comb_df.RData")


## combine pops for both sexes
all_pop <- f_pop[,3:21] + m_pop[,3:21]
all_pop <- cbind(f_pop[,1:2], all_pop)
colnames(all_pop)[3:20] <- c(1:18)
foi_comb_all <- merge(foi_mat, all_pop, by = c("x", "y")) # total female population combined with foi data
mask_values <- extract(chik_binary, foi_comb_all[, c("x", "y")])
mask_val_occ <- extract(chik_occ, foi_comb_all[, c("x", "y")])
foi_comb_all <- foi_comb_all %>%
  mutate(mask_val_bin = mask_values, mask_val_occ = mask_val_occ)
foi_comb_all <- foi_comb_all[!is.na(foi_comb_all$"1"), ]

chik_occ_shp <- st_read("chik_data_ahyoung/chik_occ_country.shp")
chik_occ_shp_sf <- st_as_sf(chik_occ_shp, wkt = "geometry", crs = 4326)
invalid_geometries <- st_is_valid(chik_occ_shp, reason = TRUE)
chik_occ_shp <- st_make_valid(chik_occ_shp)
chik_occ_shp_sf <- st_as_sf(chik_occ_shp, wkt = "geometry", crs = 4326)

foi_comb_all_sf <- st_as_sf(foi_comb_all, coords = c("x", "y"), crs = 4326)

foi_comb_all_sf <- st_join(foi_comb_all_sf, chik_occ_shp_sf, join = st_within)

foi_comb_all_df <- as.data.frame(foi_comb_all_sf)

foi_comb_bin1 <- foi_comb_all_df %>% filter(chik_bn == 1)

unique(foi_comb_bin1$geounit)
unique(foi_comb_all_df$geounit)

save(foi_comb_all_df, file = "foi_comb_all_df.RData")
save(foi_comb_all_df, file = "foi_comb_all_df_0623.RData") # without CHIK Risk

chik_risk_cdc <- read.csv(file = "cdc_chik_iso3.csv")  # 1st model (with 106 countries)

allfoi = foi_comb_df
foi_cols <- allfoi[,27:126]
allfoi <- allfoi %>% 
  mutate(
    foi_mid = ifelse(chik_bn == 0 | chik_binmap == 0, 0, apply(foi_cols, 1, median, na.rm = TRUE)),
    foi_lo  = ifelse(chik_bn == 0 | chik_binmap == 0, 0, apply(foi_cols, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))),
    foi_hi  = ifelse(chik_bn == 0 | chik_binmap == 0, 0, apply(foi_cols, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))),
    foi_sd  = ifelse(chik_bn == 0 | chik_binmap == 0, 0, apply(foi_cols, 1, sd, na.rm = TRUE))
  )

allfoi_sf <- st_as_sf(allfoi, coords = c("x", "y"), crs = crs(tsuit))
foi_sf <- st_as_sf(allfoi_sf)

usa <- ne_states(country = "United States of America", returnclass = "sf")
alaska <- usa[usa$name == "Alaska",]
within_alaska <- st_within(foi_sf, alaska, sparse = F)
foi_sf$foi_mid[within_alaska] <- 0
foi_sf$foi_lo[within_alaska] <- 0
foi_sf$foi_hi[within_alaska] <- 0
foi_sf$foi_sd[within_alaska] <- 0


save(allfoi, file = "allfoi.RData")
save(foi_sf, file = "foi_sf.RData")

# make map
foi_rast <- rasterize(foi_sf, tsuit, field = 'foi_mid')
foi_cmask <- rasterize(foi_sf, tsuit, field = 'chik_bn')

# size issue
tsuit_terra <- rast(tsuit)
foi_rast <- terra::rasterize(foi_sf, tsuit_terra, field = 'foi_mid', filename = "foi_raster.tif", overwrite = TRUE)
foi_cmask <- terra::rasterize(foi_sf, tsuit_terra, field = 'chik_bn', filename = "foi_cmask.tif", overwrite = TRUE)

# get rid of Alaska from c_mask
alaska <- st_transform(alaska, crs = crs(foi_cmask))
alaska_mask <- rasterize(alaska, foi_cmask, field = 1, background = NA) # make Alaska = 0 in the foi_cmask
foi_cmask[!is.na(alaska_mask)] <- 0

# convert the raster to SpatRaster
foi_layer1 <- rast("foi_raster.tif")
foi_layer2 <- rast("foi_cmask.tif")
# rasterize SD
sd_rast <- rasterize(foi_sf, tsuit, field = 'foi_sd')
sd_layer1 <- as(sd_rast, "SpatRaster")

# draw map
make_pixel_foi_map(foi_layer1, foi_layer2)
make_pixel_foi_map_sprast(foi_rast, foi_cmask, foi_layer1, foi_layer2)
make_pixel_sd_map_sprast(sd_rast, foi_cmask, sd_layer1, foi_layer2)
make_pixel_sd_map(sd_rast)

# foi dat with only cdc countries
# check data
summary(foi_sf$c_mask)
summary(foi_sf$mask) # 47021 NAs
foi_mask_na <- foi_sf[is.na(foi_sf$mask),]

final_mask <- (foi_sf$chik_bn == 1) & (foi_sf$chik_binmap == 1) 
foi_sf_mask <- foi_sf[final_mask, ]
foi_sf_mask <- na.omit(foi_sf_mask)
foi_comb_all <- as.data.frame(foi_sf_mask)
                       

foi_mid <- foi_sf_mask %>%
  group_by(country) %>%
  summarise(
    Lower = quantile(foi_mid, 0),
    Upper = quantile(foi_mid, 1),
    Mid   = quantile(foi_mid, 0.5),
    country = first(country),
    continent = first(continent)
  ) %>%
  ungroup()

ggplot(foi_mid, aes(y = reorder(country, Upper), x = (Lower + Upper) / 2)) +
  geom_tile(aes(height = 0.8, width = Upper - Lower, fill = Upper), color = "black") +
  labs(x = "Force of Infection", y = "Country") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_wrap(~continent, scales = "free_y", ncol = 3)


allfoi_mask1 <- allfoi_mask1[!allfoi_mask1$name %in% c("United States of America" , "China"),]

allfoi_asia <- foi_comb_all[foi_comb_all$continent == "Asia", ]
allfoi_africa <- foi_comb_all[foi_comb_all$continent == "Africa", ]
allfoi_northamerica <- foi_comb_all[foi_comb_all$continent == "North America", ]
allfoi_southamerica <- foi_comb_all[foi_comb_all$continent == "South America", ]
allfoi_europe <- foi_comb_all[foi_comb_all$continent == "Europe", ]
allfoi_oceania <- foi_comb_all[foi_comb_all$continent == "Oceania", ]

quantile(allfoi_mask1$foi_mid, c(0.5, 0, 1))
quantile(allfoi_asia$foi_mid, c(0.5, 0, 1))
quantile(allfoi_africa$foi_mid, c(0.5, 0, 1))
quantile(allfoi_northamerica$foi_mid, c(0.5, 0, 1))
quantile(allfoi_southamerica$foi_mid, c(0.5, 0, 1))
quantile(allfoi_europe$foi_mid, c(0.5, 0, 1))
quantile(allfoi_oceania$foi_mid, c(0.5, 0, 1))


save(foi_comb_all, file = "foi_comb_all_0707.RData")

#### burden function process (1: mid)
start_time <- Sys.time()

burden_1 <- age_strat_burden_psa(df = foi_comb_all, 26,50)

burden_2 <- age_strat_burden_psa(df = foi_comb_all, 51,75)

burden_3 <- age_strat_burden_psa(df = foi_comb_all, 76,100)

burden_4 <- age_strat_burden_psa(df = foi_comb_all, 101,125)

end_time <- Sys.time()

(duration <- end_time - start_time)

all_results <- c(burden_1, burden_2, burden_3, burden_4)

numbs <- length(all_results)

tot_inf<-do.call(cbind, lapply(all_results, function(df) df$updated_df$total_infection))
colnames(tot_inf) <- paste("tot_inf", seq_len(numbs))
common_cols <- all_results[[1]]$updated_df[, 1:25]
combined_burden <- cbind(common_cols, tot_inf)
tot_inf_cols <- combined_burden[,26:125]
combined_burden$med_inf <- apply(tot_inf_cols, 1, median, na.rm = T)
combined_burden$lo_inf <- apply(tot_inf_cols, 1, function(x) quantile(x, probs = 0.025))
combined_burden$hi_inf <- apply(tot_inf_cols, 1, function(x) quantile(x, probs = 0.975))
combined_burden$sd <- apply(tot_inf_cols, 1, function(x) sd(x, na.rm = TRUE))
age_cols <- all_results[[1]]$updated_df[,7:24]
tot_pop <- rowSums(age_cols[, 1:18])
combined_burden <- cbind(combined_burden, tot_pop)

save(combined_burden, file = "combined_burden_0427.RData")
save(combined_burden, file = "combined_burden_0707.RData")

save(burden_1, file = "burden_1_0709.RData")
save(burden_2, file = "burden_2_0709.RData")
save(burden_3, file = "burden_3_0709.RData")
save(burden_4, file = "burden_4_0709.RData")

## check calculation
burden_samp <- burden_1[[1]]
allinf <- burden_samp$updated_df
incid <- burden_samp$incidence_per_band
infec <- burden_samp$infection_per_band

age_g1_pop <- 22.3551016
incid_g1 <- 0.024750348
(age_g1_pop * incid_g1)

# burden by country (95% mid)
infection_by_iso3 <- combined_burden %>% group_by(country) %>%
  summarise(tot_infec_med  = sum(med_inf),
            tot_infec_lo   = sum(lo_inf),
            tot_infec_hi   = sum(hi_inf),
            iso3           = first(iso3),
            country        = first(country),
            continent      = first(continent),
            tot_pop        = sum(tot_pop)
  ) %>% 
  as.data.frame()
infec_per_100k_mid <- (infection_by_iso3$tot_infec_med/ infection_by_iso3$tot_pop)*100000
infec_per_100k_lo  <- (infection_by_iso3$tot_infec_lo/ infection_by_iso3$tot_pop)*100000
infec_per_100k_hi  <- (infection_by_iso3$tot_infec_hi/ infection_by_iso3$tot_pop)*100000
infection_by_iso3  <- cbind(infection_by_iso3, infec_per_100k_mid, infec_per_100k_lo, infec_per_100k_hi)

infection_by_iso3 <- infection_by_iso3[!infection_by_iso3$country %in% c("United States of America" , "China"),]

infection_asia <- infection_by_iso3[infection_by_iso3$continent == "Asia", ]
infection_africa <- infection_by_iso3[infection_by_iso3$continent == "Africa", ]
infection_america <- infection_by_iso3[infection_by_iso3$continent == "Americas", ]

# filter out US and China
# infection global total/ continent
combined_burden <- combined_burden[!combined_burden$name %in% c("United States of America" , "China"),]

infection_global <- combined_burden %>%
  summarise(
    tot_infec_med  = sum(med_inf, na.rm = TRUE),
    tot_infec_lo   = sum(lo_inf, na.rm = TRUE),
    tot_infec_hi   = sum(hi_inf, na.rm = TRUE),
    country        = first(name),
    continent      = first(continent),
    tot_pop        = sum(tot_pop, na.rm = TRUE),
    annual_incidence_mid = (tot_infec_med/tot_pop)*100000,
    annual_incidence_lo = (tot_infec_lo/tot_pop)*100000,
    annual_incidence_hi = (tot_infec_hi/tot_pop)*100000,
  ) %>% 
  as.data.frame()

infection_by_cont <- combined_burden %>% group_by(continent) %>%
  summarise(
    tot_infec_med  = sum(med_inf, na.rm = TRUE),
    tot_infec_lo   = sum(lo_inf, na.rm = TRUE),
    tot_infec_hi   = sum(hi_inf, na.rm = TRUE),
    country        = first(name),
    continent      = first(continent),
    tot_pop        = sum(tot_pop, na.rm = TRUE),
    annual_incidence_mid = (tot_infec_med/tot_pop)*100000,
    annual_incidence_lo = (tot_infec_lo/tot_pop)*100000,
    annual_incidence_hi = (tot_infec_hi/tot_pop)*100000,
  ) %>% 
  as.data.frame()

write.csv(infection_global, file = "infection_global.csv")
write.csv(infection_by_cont, file = "infection_by_cont.csv")
write.xlsx(infection_by_iso3, "infection_by_iso3.xlsx")

# total infection vs. annual incidence rank methods
top_total_infections <- infection_by_iso3 %>%
  arrange((tot_infec_med)) %>%
  mutate(Rank_Total = row_number()) 

top_per_100k <- infection_by_iso3 %>%
  arrange((infec_per_100k_mid)) %>%
  mutate(Rank_Per_100k = row_number())

data <- full_join(top_total_infections, top_per_100k, by = "country", suffix = c("_total", "_per_100k"))

top_countries <- data %>% group_by(continent_total) %>%
                   top_n(10, Rank_Total)

pdf("infec_incidence_rank.pdf", height = 12, width = 20)

p <- ggplot(top_countries, aes(x = Rank_Total, y = Rank_Per_100k, label = country)) +
  geom_tile(aes(fill = continent_total), color = "white") +  # Adds the tiles with white borders
  geom_text(data = top_countries, aes(label = country), size = 3.5, color = "black") +  # Adds the country labels
  scale_fill_brewer(palette = "Set2") +  # Color scale for continents
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme(axis.text.x = element_text(hjust = 1),  # Rotate x axis labels if needed
        axis.title = element_blank(),  # Remove axis titles if desired
        axis.text = element_text(size = 12),  # Adjust text size
        axis.ticks = element_blank())+
  xlab("Total infection rank")+
  ylab("Annual incidence per 100k rank")+
  labs(fill = "Continent")+
  theme_bw()

ggsave(filename = paste0("Results_figs/Infection_incidence_rank",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), p, height=10, width=15, dpi=900)

dev.off()

infection_by_iso3$country <- factor(infection_by_iso3$country, levels = infection_by_iso3$country[order(infection_by_iso3$infec_per_100k_mid, decreasing = F)])

pdf("incidence_all.pdf", width = 24, height = 13)

incidence_all <- ggplot(infection_by_iso3, aes(x = infec_per_100k_mid, y = country, fill = infec_per_100k_mid)) +
  geom_col() +  # Adjusted point size for better visibility
  geom_errorbar(aes(xmin = infec_per_100k_lo, xmax = infec_per_100k_hi), width = 0.2) +
  #scale_x_continuous(breaks = trans_breaks("log10", function(x) 10^x),
                #labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_distiller(name = "Annual incidence per 100k population", direction = 1,
                       limits = c(min(infection_by_iso3$infec_per_100k_lo), max(infection_by_iso3$infec_per_100k_hi)),
                       labels = comma) +
  theme_classic()+
  theme(axis.text.y = element_text(size = 6))+
  labs(x = "Annual incidence per 100k population",
       y = "Country")
  #facet_wrap(~continent, scales = "free_y") 
ggsave("C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/incidence_all.pdf", plot = incidence_all, width = 24, height = 13)

dev.off()

top_countries$country <- factor(top_countries$country, levels = top_countries$country[order(top_countries$infec_per_100k_mid_per_100k, decreasing = F)])

incidence_selected <- ggplot(top_countries, aes(x = infec_per_100k_mid_per_100k, y = country, fill = infec_per_100k_mid_per_100k)) +
  geom_col() +  # Adjusted point size for better visibility
  geom_errorbar(aes(xmin = infec_per_100k_lo_per_100k, xmax = infec_per_100k_hi_per_100k), width = 0.2) +
  #scale_x_continuous(breaks = trans_breaks("log10", function(x) 10^x),
  #labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_distiller(name = "Annual incidence per 100k population", direction = 1,
                       limits = c(min(top_countries$infec_per_100k_lo_per_100k), max(top_countries$infec_per_100k_hi_per_100k)),
                       labels = comma) +
  theme_classic()+
  theme(axis.text.y = element_text(size = 6))+
  labs(x = "Annual incidence per 100k population",
       y = "Country")+
  facet_wrap(~continent_total, scales = "free_y") 


ggsave(filename = paste0("Results_figs/incidence_all",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), incidence_all, height=8, width=13, dpi=900)

ggsave(filename = paste0("Results_figs/incidence_selected",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), incidence_selected, height=8, width=13, dpi=900)

# cumulative infection graph
cum_df <- infection_by_iso3 %>% 
  ungroup()%>%
  arrange(desc(tot_infec_med))%>%
  mutate(cum_sum  = cumsum(tot_infec_med),
         cum_perc = cum_sum / sum(tot_infec_med) * 100)

cum_df$country <- factor(cum_df$country, levels = cum_df$country)

color <- pal_lancet("lanonc")(5)

cum_graph <- ggplot(cum_df, aes(x = reorder(country, -tot_infec_med), y = tot_infec_med, fill = continent)) +
  geom_bar(stat = "identity", alpha = 0.6) +
  geom_errorbar(aes(ymin = tot_infec_lo, ymax = tot_infec_hi), alpha = 0.6)+
  geom_line(aes(y = cum_perc * max(tot_infec_hi) / 100, group = 1), color = "darkgrey", size = 0.6) +
  geom_point(aes(y = cum_perc * max(tot_infec_hi) / 100), color = "darkgrey", size = 1, alpha = 0.6) +
  scale_y_continuous(
    name = "Total Infections",
    trans = "sqrt", 
    breaks = c(0, 1000000, 4000000, seq(4000000, max(12000000, na.rm = TRUE), by = 4000000)),
    labels = scales::comma,
    sec.axis = sec_axis(~ . / max(cum_df$tot_infec_hi), name = "Cumulative Percentage", labels = scales::percent_format(),
                        breaks = c(0.10, 0.25, 0.5, 0.75, 1))
  ) +
  scale_fill_manual(values = color)+
  #scale_y_log10(scales::comma)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  xlab("Country")+
  labs(fill = "Continent")

ggsave(filename = paste0("Results_figs/cum_graph",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), cum_graph, height=12, width=20, dpi=900)


### sub-burdens 
lhs_sample <- as.data.frame(lhs_sample)
med_care   <- as.data.frame(med_care)
sub_burden_psa <- age_subinf_lineage_psa(infection_by_iso3, lhs_sample, med_care)

fatal_psa <- sub_burden_psa$fatal
symp_psa  <- sub_burden_psa$symp
symp_opd_psa <- sub_burden_psa$symp_opd
ipd_tp_psa <- sub_burden_psa$ipd_tp
opd_tp_psa <- sub_burden_psa$opd_tp
asymp_psa <- sub_burden_psa$asymp
lt_psa    <- sub_burden_psa$chronic
hosp_psa  <- sub_burden_psa$hosp
nonhosp_psa <- sub_burden_psa$nonhosp
yll_psa   <- sub_burden_psa$yll
yldch_psa <- sub_burden_psa$yld_ch
yldhosp_psa <- sub_burden_psa$yld_hosp
yldnonhosp_psa <- sub_burden_psa$yld_nonhosp
daly_psa  <- sub_burden_psa$daly

fatal_psa$tot_pop <- infection_by_iso3$tot_pop
symp_psa$tot_pop  <- infection_by_iso3$tot_pop
symp_opd_psa$tot_pop <- infection_by_iso3$tot_pop
ipd_tp_psa$tot_pop <- infection_by_iso3$tot_pop
opd_tp_psa$tot_pop <- infection_by_iso3$tot_pop
asymp_psa$tot_pop <- infection_by_iso3$tot_pop
lt_psa$tot_pop    <- infection_by_iso3$tot_pop
hosp_psa$tot_pop  <- infection_by_iso3$tot_pop
nonhosp_psa$tot_pop <- infection_by_iso3$tot_pop
yll_psa$tot_pop <- infection_by_iso3$tot_pop
yldch_psa$tot_pop <- infection_by_iso3$tot_pop
yldhosp_psa$tot_pop <- infection_by_iso3$tot_pop
yldnonhosp_psa$tot_pop <- infection_by_iso3$tot_pop
daly_psa$tot_pop <- infection_by_iso3$tot_pop

fatal_glob_count <- sub_inf_global_count(fatal_psa)
symp_glob_count <- sub_inf_global_count(symp_psa)
symp_opd_glob_count <- sub_inf_global_count(symp_opd_psa)
ipd_tp_glob_count <- sub_inf_global_count(ipd_tp_psa)
opd_tp_glob_count <- sub_inf_global_count(opd_tp_psa)
asymp_glob_count <- sub_inf_global_count(asymp_psa)
lt_glob_count <- sub_inf_global_count(lt_psa)
hosp_glob_count <- sub_inf_global_count(hosp_psa)
nonhosp_glob_count <- sub_inf_global_count(nonhosp_psa)
yll_glob_count <- sub_inf_global_count(yll_psa)
yldch_glob_count <- sub_inf_global_count(yldch_psa)
yldhosp_glob_count <- sub_inf_global_count(yldhosp_psa)
yldnonhosp_glob_count <- sub_inf_global_count(yldnonhosp_psa)
daly_glob_count <- sub_inf_global_count(daly_psa)

fatal_global <- fatal_glob_count$global
fatal_cont <- fatal_glob_count$cont
symp_global <- symp_glob_count$global
symp_cont <- symp_glob_count$cont
symp_nat <- symp_glob_count$country
asymp_global <- asymp_glob_count$global
asymp_cont <- asymp_glob_count$cont
lt_global <- lt_glob_count$global
lt_cont <- lt_glob_count$cont
hosp_global <- hosp_glob_count$global
hosp_cont <- hosp_glob_count$cont
nonhosp_global <- nonhosp_glob_count$global
nonhosp_cont <- nonhosp_glob_count$cont
daly_global <- daly_glob_count$global
daly_cont <- daly_glob_count$cont
daly_country <- daly_glob_count$country
yll_cont <-yll_glob_count$cont
yldch_cont <- yldch_glob_count$cont
yldhosp_cont <- yldhosp_glob_count$cont
yldnonhosp_cont <- yldnonhosp_glob_count$cont
yll_cont$rate_mid <- (yll_cont$tot_med / yll_cont$tot_pop)*1000000
yll_cont$rate_lo <- (yll_cont$tot_lo / yll_cont$tot_pop)*1000000
yll_cont$rate_hi <- (yll_cont$tot_hi / yll_cont$tot_pop)*1000000
yldch_cont$rate_mid <- (yldch_cont$tot_med / yldch_cont$tot_pop)*1000000
yldch_cont$rate_lo <- (yldch_cont$tot_lo / yldch_cont$tot_pop)*1000000
yldch_cont$rate_hi <- (yldch_cont$tot_hi / yldch_cont$tot_pop)*1000000
yldnonhosp_cont$rate_mid <- (yldnonhosp_cont$tot_med / yldnonhosp_cont$tot_pop)*1000000
yldnonhosp_cont$rate_lo <- (yldnonhosp_cont$tot_lo / yldnonhosp_cont$tot_pop)*1000000
yldnonhosp_cont$rate_hi <- (yldnonhosp_cont$tot_hi / yldnonhosp_cont$tot_pop)*1000000
yldhosp_cont$rate_mid <- (yldhosp_cont$tot_med / yldhosp_cont$tot_pop)*1000000
yldhosp_cont$rate_lo <- (yldhosp_cont$tot_lo / yldhosp_cont$tot_pop)*1000000
yldhosp_cont$rate_hi <- (yldhosp_cont$tot_hi / yldhosp_cont$tot_pop)*1000000
daly_cont$rate_mid <- (daly_cont$tot_med / daly_cont$tot_pop)*1000000
daly_cont$rate_lo <- (daly_cont$tot_lo / daly_cont$tot_pop)*1000000
daly_cont$rate_hi <- (daly_cont$tot_hi / daly_cont$tot_pop)*1000000
yll_cont$type <- "YLL"
yldch_cont$type <- "YLD_chronic"
yldnonhosp_cont$type <- "YLD_nonhosp"
yldhosp_cont$type <- "YLD_hosp"
daly_cont$type <- "DALY"

glob_sub_burden_df <- cbind(fatal_global, symp_global, asymp_global, lt_global,
                            hosp_global, nonhosp_global, daly_global)
cont_sub_burden_df <- cbind(fatal_cont, symp_cont, asymp_cont, lt_cont,
                            hosp_cont, nonhosp_cont, daly_cont)

yll_yld <- rbind(yll_cont, yldch_cont, yldnonhosp_cont, yldhosp_cont, daly_cont)

yll_yld_all <- yll_yld %>% summarise()

write.csv(glob_sub_burden_df, file = "glob_sub_burden_df_s4.csv")
write.csv(cont_sub_burden_df, file = "cont_sub_burden_df_s4.csv")
write.csv(yll_yld, file = "yll_yld.csv")

daly_cont$type <- "DALY"
yll_cont$type <- "YLL"
yldch_cont$type <- "YLD_chronic"
yldhosp_cont$type <- "YLD_hosp"
yldnonohosp_cont$type <- "YLD_nonhosp"

stack_daly <- rbind(yll_cont, yldch_cont, yldhosp_cont, yldnonhosp_cont)
stack_daly <- stack_daly %>% 
                group_by(continent) %>% 
                   mutate(tot_continent = sum(tot_med),
                          percent       = (tot_med / tot_continent)*100) %>%
                            ungroup()

# stack bar graph by continent (DALY)
pdf("daly_continent.pdf", width = 10, height = 12)

color <- c("#00468B99", "#ED000099", "#42B54099", "#FDAF9199")
daly_cont <- ggplot(stack_daly, aes(x = continent, y = tot_med, fill = type))+
  geom_bar(position = "dodge", stat = "identity")+
  geom_errorbar(position = position_dodge(width = 0.9), aes(ymin = tot_lo, ymax = tot_hi), width = 0.2, alpha = 0.6) +
  scale_fill_manual(values = color,
                    name = "DALY composition",  # Change the legend title
                    labels = c("YLD chronic", "YLD hospitalisation", "YLD non-hospitalisation", "YLL"))+
  xlab("Continent")+
  ylab("DALY")+
  scale_y_continuous(label = comma)+
  theme_bw()

daly_cont_ribbon <- ggplot(stack_daly, aes(x = continent, y = tot_med, group = type, color = type, fill = type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = tot_lo, ymax = tot_hi), alpha = 0.2) +
  scale_color_manual(values = color,
                     name = "DALY composition",  # Change the legend title
                     labels = c("YLD chronic", "YLD hospitalisation", "YLD non-hospitalisation", "YLL")) +
  scale_fill_manual(values = color,
                    name = "DALY composition",  # Change the legend title
                    labels = c("YLD chronic", "YLD hospitalisation", "YLD non-hospitalisation", "YLL")) +
  xlab("Continent") +
  ylab("DALY") +
  scale_y_continuous(labels = comma)+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x))+
  theme_bw()




#-------------------------------------------------------------------------------
## age stratified burden

agelist1 <- lapply(burden_1, function(x) x$infection_per_band)
agelist2 <- lapply(burden_2, function(x) x$infection_per_band)
agelist3 <- lapply(burden_3, function(x) x$infection_per_band)
agelist4 <- lapply(burden_4, function(x) x$infection_per_band)

agelist1 <- lapply(agelist1, as.data.frame)
agelist2 <- lapply(agelist2, as.data.frame)
agelist3 <- lapply(agelist3, as.data.frame)
agelist4 <- lapply(agelist4, as.data.frame)


col_names <- c(1:18)
agelist1 <- lapply(agelist1, function(df){
  colnames(df) <- col_names
  return(df)
})
agelist2 <- lapply(agelist2, function(df){
  colnames(df) <- col_names
  return(df)
})
agelist3 <- lapply(agelist3, function(df){
  colnames(df) <- col_names
  return(df)
})
agelist4 <- lapply(agelist4, function(df){
  colnames(df) <- col_names
  return(df)
})

comb_df_age1 <- list()
comb_df_age2 <- list()
comb_df_age3 <- list()
comb_df_age4 <- list()

for(i in col_names){
  comb_df_age1[[i]] <- do.call(cbind, lapply(agelist1, function(df) df[, i, drop = F]))
}
for(i in col_names){
  comb_df_age2[[i]] <- do.call(cbind, lapply(agelist2, function(df) df[, i, drop = F]))
}
for(i in col_names){
  comb_df_age3[[i]] <- do.call(cbind, lapply(agelist3, function(df) df[, i, drop = F]))
}
for(i in col_names){
  comb_df_age4[[i]] <- do.call(cbind, lapply(agelist4, function(df) df[, i, drop = F]))
}

num_dfs <- length(comb_df_age1)
all_agelist <- vector("list", length = num_dfs)

for(i in 1:num_dfs){
  all_agelist[[i]] <- do.call(cbind, list(comb_df_age1[[i]],
                                          comb_df_age2[[i]],
                                          comb_df_age3[[i]],
                                          comb_df_age4[[i]]))
}

all_agelist <- lapply(all_agelist, function(df) {
  mid_val <- apply(df[,1:100], 1, function(x) quantile(x, 0.5, na.rm = T))
  lo_val <- apply(df[,1:100], 1, function(x) quantile(x, 0.025, na.rm = T))
  hi_val <- apply(df[,1:100], 1, function(x) quantile(x, 0.975, na.rm = T))
  df$lo_val <- lo_val
  df$hi_val <- hi_val
  df$mid_val <- mid_val
  
  df <- df[!(df$name %in% c("United States of America", "China")), ]
  
  return(df)
  
})

## annual incidence
tot_pop <- foi_comb_f_df[,101:118] + foi_comb_m_df[,101:118]
tot_pop$tot_pop <- rowSums(tot_pop[,1:18], na.rm = T)
common_cols <- foi_comb_f_df[,119:129]
tot_pop <- cbind(common_cols, tot_pop)
tot_pop <- tot_pop[tot_pop$geometry %in% foi_comb_all_cdc$geometry, ]

all_agelist <-lapply(all_agelist, function(df){
  cbind(tot_pop, df)
  
})

extract_vals <- lapply(1: length(all_agelist), function(i) {
  df <- all_agelist[[i]]
  common_cols <- df[,1:10]
  age_group <- df[,12:29]
  tot_pop <- df[,30]
  result_df <- data.frame(age_group = i,
                          mid       = df$mid,
                          lo        = df$lo,
                          hi        = df$hi
  )
  result_df <- cbind(common_cols, age_group, tot_pop, result_df)
  return(result_df)
})


extract_vals <- lapply(1: length(all_agelist), function(i) {
  df <- all_agelist[[i]]
  data.frame(age_group = i,
             mid       = df$mid,
             lo        = df$lo,
             hi        = df$hi)
})

all_agelist_global <- lapply(1: length(extract_vals), function(i){
  df <- extract_vals[[i]]
  data.frame(age_group = i,
             mid = sum(df$mid),
             lo  = sum(df$lo),
             hi  = sum(df$hi)
  )
})

all_agelist_cont <- lapply(seq_along(extract_vals), function(i){
  df <- extract_vals[[i]]
  
  pop_cols <- c("f 0", "f 1", paste0("f ", seq(5, 80, by = 5)))
  
  sum_df <- df %>% group_by(continent, age_group) %>%
      summarise(age_group = i,
             mid = sum(mid),
             lo  = sum(lo),
             hi  = sum(hi)
  )
  
  pop_sum_df <- df %>%
    group_by(continent) %>%
    summarise(across(all_of(pop_cols), sum, na.rm = TRUE),
              .groups = 'drop')
  
  final_df <- left_join(sum_df, pop_sum_df, by = "continent")
  
  return(final_df)
})

combined_agedf <- do.call(rbind, all_agelist_global)
combined_agedf_cont <- do.call(rbind, all_agelist_cont)
unique_continent_data <- combined_agedf_cont %>%
  dplyr::group_by(continent) %>%
  dplyr::slice(1) %>%
  dplyr::select(c("continent", "f 0", "f 1", paste0("f ", seq(5, 80, by = 5))))

custom_order <- c("f 0", "f 1", paste0("f ", seq(5, 80, by = 5)))
long_pop <- unique_continent_data %>%
  tidyr::pivot_longer(
    cols = starts_with("f "),  # This selects all columns that start with 'f '
    names_to = "pop_group",
    values_to = "pop_size"
  ) %>%
  mutate(
    pop_group = factor(pop_group, levels = custom_order)  # Set the factor levels to the custom order
  ) %>%
  arrange(pop_group, continent) 

combined_agedf_cont <- cbind(combined_agedf_cont[,1:5], long_pop[,3])

combined_agedf$group <- factor(combined_agedf$group, levels = c("<1", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
                                                                 "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
                                                                 "70-74", "75-79", "80+"))
age_labs <- c("<1", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
              "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
              "70-74", "75-79", "80+")

pdf("age_strat_infection.pdf", height = 12, width = 18)
age_strat_infec <- ggplot(combined_agedf, aes(x = group))+
  geom_ribbon(aes(ymin = lo, ymax = hi, group = 1), fill = "lightblue")+
  geom_line(aes(y = mid, group = 1), color = "black")+
  scale_y_continuous(labels = comma) +
  xlab("Age group") +
  ylab("Total infections")+
  theme_bw()
dev.off()

color <- pal_lancet("lanonc")(5)
pdf("agestrat_cont.pdf", height = 10, width = 18)
agestrat_infect_cont <- ggplot(combined_agedf_cont, aes(x = age_group, y = mid, fill = continent))+
  geom_area(alpha = 0.5)+
  #geom_ribbon(aes(ymin = lo, ymax = hi, fill = continent), alpha = 0.3)+
  #geom_crossbar(aes(ymin = lo, ymax = hi), width = 0.2)+
  scale_y_continuous(labels = comma) +
  xlab("Age group") +
  ylab("Total infections")+
  scale_x_continuous(breaks = 1:18, labels = age_labs)+
  scale_fill_manual(values = color)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


agestrat_infect_cont <- ggplot(combined_agedf_cont, aes(x = age_group, y = mid, group = continent, color = continent, fill = continent)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  scale_color_manual(values = color,
                     name = "Continent") +
  scale_fill_manual(values = color,
                    name = "Continent") +
  xlab("Age group") +
  ylab("Total infection") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = 1:18, labels = age_labs)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_agedf_cont$incidence_mid <- (combined_agedf_cont$mid/ combined_agedf_cont$pop_size)*100000
combined_agedf_cont$incidence_lo <- (combined_agedf_cont$lo/ combined_agedf_cont$pop_size)*100000
combined_agedf_cont$incidence_hi <- (combined_agedf_cont$hi/ combined_agedf_cont$pop_size)*100000

agestrat_inc_cont <- ggplot(combined_agedf_cont, aes(x = age_group, y = incidence_mid, group = continent, color = continent, fill = continent)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = incidence_lo, ymax = incidence_hi), alpha = 0.2) +
  scale_color_manual(values = color,
                     name = "Continent") +
  scale_fill_manual(values = color,
                    name = "Continent") +
  xlab("Age group") +
  ylab("Annual incidence per 100 000 people") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = 1:18, labels = age_labs)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf("combined_age_strat_graph.pdf", height = 5, width = 20)
ggarrange(agestrat_infect_cont, agestrat_inc_cont, daly_cont_ribbon,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()


all_agelist_global <- lapply(1: length(extract_vals), function(i){
  df <- extract_vals[[i]]
  data.frame(age_group = i,
             mid = sum(df$mid),
             lo  = sum(df$lo),
             hi  = sum(df$hi)
  )
})

all_agelist_cont <- lapply(1: length(extract_vals), function(i){
  df <- extract_vals[[i]]
  df <- df %>% group_by(continent) %>%
    summarise(age_group = i,
              mid = sum(mid),
              lo  = sum(lo),
              hi  = sum(hi),
              tot_pop = sum(tot_pop)
    )
  df$incidence_mid <- (df$mid/ df$tot_pop)*100000
  df$incidence_lo <- (df$lo/ df$tot_pop)*100000
  df$incidence_hi <- (df$hi/ df$tot_pop)*100000
  return(df)
})

combined_agedf <- do.call(rbind, all_agelist_global)
combined_agedf_cont <- do.call(rbind, all_agelist_cont)

pop_size <- as.data.frame(colSums(extract_vals[[1]][,11:28]))

combined_agedf <- cbind(combined_agedf, pop_size)
names(combined_agedf) <- c("age_group","mid", "lo","hi", "tot_pop")
combined_agedf$group <- c("<1", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
                          "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
                          "70-74", "75-79", "80+")
combined_agedf$group <- factor(combined_agedf$group, levels = c("<1", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
                                                                "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
                                                                "70-74", "75-79", "80+"))
combined_agedf$incidence_mid <- (combined_agedf$mid/ combined_agedf$tot_pop)*100000
combined_agedf$incidence_lo <- (combined_agedf$lo/ combined_agedf$tot_pop)*100000
combined_agedf$incidence_hi <- (combined_agedf$hi/ combined_agedf$tot_pop)*100000



pdf("age_strat_incidence.pdf", height = 12, width = 18)
agestrat_incidence <- ggplot(combined_agedf, aes(x = group))+
  geom_ribbon(aes(ymin = incidence_lo, ymax = incidence_hi, group = 1), fill = "lightblue")+
  geom_line(aes(y = incidence_mid, group = 1), color = "black")+
  scale_y_continuous(labels = comma) +
  xlab("Age group") +
  ylab("Annual incidence per 100 000 people")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

agestrat_inc_cont <- ggplot(combined_agedf_cont, aes(x = age_group, y = incidence_mid, fill = continent))+
  geom_area(alpha = 0.5)+
  #geom_ribbon(aes(ymin = lo, ymax = hi, fill = continent), alpha = 0.3)+
  #geom_crossbar(aes(ymin = lo, ymax = hi), width = 0.2)+
  scale_y_continuous(labels = comma) +
  xlab("Age group") +
  ylab("Annual incidence per 100 000 people")+
  scale_x_continuous(breaks = 1:18, labels = age_labs)+
  scale_fill_manual(values = color)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined_agedf_cont$incidence_mid <- (combined_agedf_cont$mid/ combined_agedf_cont$pop_size)*100000
combined_agedf_cont$incidence_lo <- (combined_agedf_cont$lo/ combined_agedf_cont$pop_size)*100000
combined_agedf_cont$incidence_hi <- (combined_agedf_cont$hi/ combined_agedf_cont$pop_size)*100000

agestrat_inc_cont <- ggplot(combined_agedf_cont, aes(x = age_group, y = incidence_mid, group = continent, color = continent, fill = continent)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = incidence_lo, ymax = incidence_hi), alpha = 0.2) +
  scale_color_manual(values = color,
                     name = "Continent") +
  scale_fill_manual(values = color,
                    name = "Continent") +
  xlab("Age group") +
  ylab("Annual incidence per 100 000 people") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = 1:18, labels = age_labs)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf("combined_age_strat_graph.pdf", height = 5, width = 20)
ggarrange(agestrat_infect_cont, agestrat_inc_cont, daly_cont_ribbon,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

dev.off()
save(all_agelist, file = "all_agelist.RData")

## age stratified DALYs
all_agelist_country <- lapply(seq_along(extract_vals), function(i, df){
  df <- extract_vals[[i]]
  sum <- df %>% group_by(iso3, continent, name) %>% 
        summarise(
                  mid = sum(mid),
                  lo  = sum(lo),
                  hi  = sum(hi),
                  .groups = 'drop'
  )
  sum$age_group <- i
  return(sum)
})

combined_agedf_country <- do.call(rbind, all_agelist_country)
names(combined_agedf_country) <- c("iso3", "continent", "country", "tot_infec_med", "lo", "hi", "age_group")

age_sub_burden <- age_subinf_psa(combined_agedf_country, lhs_sample)

age_sub_burden_all <- lapply(seq_along(age_sub_burden), function(i, df){
  df <- age_sub_burden[[i]]
  age_group <- rep(1:18, each = 106)
  df$age_group <- age_group
  return(df)
})
names <- c("fatal", "surv", "symp", "asymp", "lt","ac", "hosp",
           "yll", "yld_ch", "yld_hosp", "yld_nonhosp", "daly")
names(age_sub_burden_all) <- names

age_yll_psa   <- age_sub_burden_all$yll
age_yldch_psa <- age_sub_burden_all$yld_ch
age_yldhosp_psa <- age_sub_burden_all$yld_hosp
age_yldnonhosp_psa <- age_sub_burden_all$yld_nonhosp

age_yll_glob_count <- age_inf_global_count(age_yll_psa)
age_yldch_glob_count <- age_inf_global_count(age_yldch_psa)
age_yldhosp_glob_count <- age_inf_global_count(age_yldhosp_psa)
age_yldnonhosp_glob_count <- age_inf_global_count(age_yldnonhosp_psa)

age_yll_glob <- age_yll_glob_count$global
age_yll_glob$type <- "YLL"
age_yldch_glob <- age_yldch_glob_count$global
age_yldch_glob$type <- "YLD_chronic"
age_yldhosp_glob <- age_yldhosp_glob_count$global
age_yldhosp_glob$type <- "YLD_hosp"
age_yldnonhosp_glob <- age_yldnonhosp_glob_count$global
age_yldnonhosp_glob$type <- "YLD_nonhosp"

age_stack_daly <- rbind(age_yll_glob, age_yldch_glob, age_yldhosp_glob, age_yldnonhosp_glob)

# stack bar graph by continent (DALY)
pdf("daly_continent.pdf", width = 10, height = 12)

age_stack_daly$group <- c(rep(c("<1", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
                          "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
                          "70-74", "75-79", "80+"), times = 4))
age_stack_daly$group <- factor(age_stack_daly$group, levels = c("<1", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
                                                                    "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
                                                                    "70-74", "75-79", "80+"))


age_strat_daly <- ggplot(age_stack_daly, aes(x = group, y = tot_med, fill = type))+
  geom_bar(stat = "identity", alpha = 0.8, linewidth = .5)+
  scale_fill_viridis(discrete = T)+
  xlab("Age group")+
  ylab("Total DALY count")+
  scale_y_continuous(labels = comma)+
  theme_bw()

pdf("age_strat_graphs.pdf", height = 12, width = 26)
age_strat <- ggarrange(age_strat_infec, age_strat_incidence, age_strat_daly,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()

ggsave(filename = paste0("Results_figs/age_strat_graphs",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), age_strat, height=12, width=28, dpi=900)

  