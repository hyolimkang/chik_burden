# load data
chik_occ_shp <- st_read("chik_data_ahyoung/chik_occ_country.shp")
chik_binary <- raster("chik_data_ahyoung/CHIK_binmap_2024_04_24.tif") # Risk map

load("worldpop/f_pop.RData")
load("worldpop/m_pop.RData")

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

samp_rast <- ensemble_raster[[1]]

extracted_values <- exact_extract(samp_rast, chik_occ_shp, fun = NULL, include_xy = TRUE)
combined_values <- lapply(seq_along(extracted_values), function(i) {
  df <- extracted_values[[i]]
  if (nrow(df) > 0) {
    df$country <- chik_occ_shp$geounit[i]
    df$continent <- chik_occ_shp$contnnt[i]
    df$region <- chik_occ_shp$region[i]
    df$iso3 <- chik_occ_shp$iso_a3[i]
    df$chik_bn <- chik_occ_shp$chik_bn[i]
  }
  return(df)
})

# Combine the list of data frames into one data frame
combined_values_df <- do.call(rbind, combined_values)

# population data cleaning
all_pop <- f_pop[,3:21] + m_pop[,3:21]
all_pop <- cbind(f_pop[,1:2], all_pop)
colnames(all_pop)[3:20] <- c(1:18)
vars_to_rast <-  c(as.character(1:18), "tot")
raster_stack <- stack()
for (var in vars_to_rast) {
  temp_raster <- rasterFromXYZ(all_pop[, c("x", "y", var)])
  names(temp_raster) <- var
  raster_stack <- stack(raster_stack, temp_raster)
}

crs_string <- "+proj=longlat +datum=WGS84"
# Set CRS for raster stack
crs(raster_stack) <- crs_string
pop_extract_list <- exact_extract(raster_stack, chik_occ_shp, include_xy = TRUE)

pop_rbind_list <- do.call(rbind, pop_extract_list)

## for 99 layers
extract_fun <- function(r) {
  extracted_values <- exact_extract(r, chik_occ_shp, fun = NULL, include_xy = TRUE)
  combined_values <- lapply(seq_along(extracted_values), function(i) {
    df <- extracted_values[[i]]
    if (nrow(df) > 0) {
      df <- dplyr::select(df, value, x, y)  # Select only the value column
    }
    return(df)
  })
  combined_values_df <- do.call(rbind, combined_values)
  return(combined_values_df)
}

remaining_values_list <- lapply(ensemble_raster, extract_fun)

# Prepare the list of data frames for 
for (i in seq_along(remaining_values_list)) {
  col_name <- paste0("foi", i)
  remaining_values_list[[i]] <- rename(remaining_values_list[[i]], !!col_name := value)
}

# Extracting x, y coordinates from the first element of the list
xy_coords <- remaining_values_list[[1]][, c("x", "y")]

# Combining values by cbind
values_combined <- do.call(cbind, lapply(remaining_values_list, function(df) dplyr::select(df, -x, -y)))
combined_df <- cbind(xy_coords, combined_values_df[,5:9], pop_rbind_list[,1:19], values_combined)
colnames(combined_df)[8:25] <- as.character(1:18)

# combine with chik_bin map
extract_binmap <- exact_extract(chik_binary, chik_occ_shp, fun = NULL, include_xy = TRUE)
binmap_vals <- do.call(rbind, extract_binmap)

chik_binmap <- binmap_vals[,1]
foi_comb_df <- cbind(chik_binmap, combined_df)

foi_comb_df <- foi_comb_df[!is.na(foi_comb_df$"1"), ]
foi_comb_df <- foi_comb_df[!is.na(foi_comb_df$"foi1"), ]
foi_comb_bin1 <- foi_comb_df %>% filter(chik_bn == 1)

unique(foi_comb_bin1$country)

save(combined_df, file = "combined_foi_df.RData")
