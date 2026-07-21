# ==============================================================================
# MASTER PROJECT ARCHIVE: NCR HABITAT & CONNECTIVITY DATA GENERATION
# Includes: Raster Extraction, Structural Metrics, Grainscape LCP, 
# Euclidean Water Barriers, and Shiny Dashboard
# ==============================================================================

# ------------------------------------------------------------------------------
# PART 1: BULLETPROOF DEPENDENCY INSTALLATION & SETUP
# ------------------------------------------------------------------------------
cat("Starting dependency check and installation...\n")

options(download.file.method = "wininet")
options(repos = c(CRAN = "https://cloud.r-project.org/"))
options(install.packages.check.source = "no")

required_packages <- c("terra", "sf", "dplyr", "data.table", "grainscape", "igraph", "shiny", "shinydashboard", "leaflet")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing missing package: %s\n", pkg))
    install.packages(pkg, type = "binary")
  }
}

library(terra)
library(sf)
library(dplyr)
library(data.table)
library(grainscape)
library(igraph)
library(shiny)
library(shinydashboard)
library(leaflet)

# Optimize terra memory settings to prevent crashes
terraOptions(memfrac = 0.8, progress = 10) 
cat("All dependencies successfully installed and loaded.\n")


# ------------------------------------------------------------------------------
# PART 2: BASE DATA EXTRACTION (WATER, CONCRETE, & GREEN SPACES)
# ------------------------------------------------------------------------------
cat("\n--- PART 2: BASE DATA EXTRACTION ---\n")
evi_path <- "D:/Data/EVI_MultiClass.tif"
if (!file.exists(evi_path)) stop("Error: Cannot find GeoTIFF at D:/Data/EVI_MultiClass.tif.")
evi_raster <- rast(evi_path)
pixel_area_sqm <- res(evi_raster)[1] * res(evi_raster)[2]

cat("Extracting Water and Concrete boundaries...\n")
water_mask <- ifel(evi_raster == 1, 1, NA)
paved_mask <- ifel(evi_raster == 2, 1, NA)

cat("Converting Water and Concrete to Web Vectors...\n")
water_poly <- as.polygons(water_mask, dissolve = TRUE)
water_sf <- st_as_sf(water_poly) %>% st_set_crs(32651) %>% st_transform(4326) %>% st_simplify(dTolerance = 0.00005)
st_write(water_sf, "D:/Data/NCR_Water.geojson", driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)

paved_poly <- as.polygons(paved_mask, dissolve = TRUE)
paved_sf <- st_as_sf(paved_poly) %>% st_set_crs(32651) %>% st_transform(4326) %>% st_simplify(dTolerance = 0.00005)
st_write(paved_sf, "D:/Data/NCR_Concrete.geojson", driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)

cat("Isolating green spaces (Classes 3 & 4) and applying Rook's Connectivity...\n")
green_mask <- ifel(evi_raster %in% c(3, 4), 1, NA)
gc() 

green_patches <- patches(green_mask, directions = 4, filename = "D:/Data/EVI_Patches_Rook.tif", overwrite = TRUE)


# ------------------------------------------------------------------------------
# PART 3: STRUCTURAL METRICS & HABITAT CHARACTERISTICS
# ------------------------------------------------------------------------------
cat("\n--- PART 3: CALCULATING HABITAT CHARACTERISTICS ---\n")
cat("Identifying Core Area Pixels...\n")
focal_matrix <- matrix(1, nrow = 3, ncol = 3)
core_calc <- focal(green_mask, w = focal_matrix, fun = sum, na.rm = FALSE)
core_mask <- ifel(core_calc == 9, 1, NA)
core_patches <- mask(green_patches, core_mask)

# Compile statistics
all_patch_stats <- freq(green_patches)
setDT(all_patch_stats) 
setnames(all_patch_stats, c("layer", "patch_id", "total_pixel_count"))

core_stats <- freq(core_patches)
setDT(core_stats)
setnames(core_stats, c("layer", "patch_id", "core_pixel_count"))

master_table <- merge(all_patch_stats, core_stats[, .(patch_id, core_pixel_count)], by = "patch_id", all.x = TRUE)
master_table[is.na(core_pixel_count), core_pixel_count := 0]

master_table[, Patch_Area_sqm := total_pixel_count * pixel_area_sqm]
master_table[, Core_Area_sqm := core_pixel_count * pixel_area_sqm]
master_table[, Edge_Area_sqm := Patch_Area_sqm - Core_Area_sqm]
master_table[, has_core := ifelse(core_pixel_count > 0, TRUE, FALSE)]
master_table[, is_islet := ifelse(core_pixel_count == 0, TRUE, FALSE)]

veg_stats <- zonal(evi_raster, green_patches, fun = "modal", na.rm = TRUE)
setDT(veg_stats); setnames(veg_stats, c("patch_id", "veg_class"))
master_table <- merge(master_table, veg_stats, by = "patch_id", all.x = TRUE)
master_table[, Veg_Type := ifelse(veg_class == 4, "Dense Canopy", "Sparse Vegetation")]

cat("Converting ALL patches to vectors (This takes a while)...\n")
all_polygons <- as.polygons(green_patches, dissolve = TRUE)
all_sf <- st_as_sf(all_polygons)
names(all_sf)[1] <- "patch_id"
final_sf <- merge(all_sf, master_table, by = "patch_id")

# Save Web-Optimized GeoJSON and RDS
green_sf_web <- final_sf %>% filter(total_pixel_count >= 10) %>% st_set_crs(32651) %>% st_transform(4326) %>% st_simplify(dTolerance = 0.00005)
st_write(green_sf_web, "D:/Data/NCR_GreenSpaces_Data.geojson", driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)

cat("Saving mathematically corrected RDS file for Leaflet...\n")
patches_wgs84 <- final_sf %>% st_set_crs(32651) %>% st_transform(4326)
saveRDS(patches_wgs84, "D:/Data/All_MetroManila_Patches.rds")


# ------------------------------------------------------------------------------
# PART 4: GRAINSCAPE LEAST-COST PATH CONNECTIVITY
# ------------------------------------------------------------------------------
cat("\n--- PART 4: GRAINSCAPE LCP CONNECTIVITY ---\n")
green_data <- st_read("D:/Data/NCR_GreenSpaces_Data.geojson", quiet = TRUE)
patch_raster <- rast("D:/Data/EVI_Patches_Rook.tif")

rcl_matrix <- matrix(c(1, 10, 2, 5, 3, 2, 4, 1), ncol = 2, byrow = TRUE)
cost_raster <- classify(evi_raster, rcl_matrix, others = NA, filename = "D:/Data/Cost_Surface.tif", overwrite = TRUE)
gc()

grainscape_scenarios <- list(
  ideal = green_data %>% filter(Patch_Area_sqm >= 8000),
  core = green_data %>% filter(Patch_Area_sqm >= 8000, Core_Area_sqm > 0),
  master = green_data 
)

for (name in names(grainscape_scenarios)) {
  cat(sprintf("\nProcessing Grainscape scenario: %s...\n", name))
  ids <- grainscape_scenarios[[name]]$patch_id
  binary_mask <- classify(patch_raster, matrix(c(ids, rep(1, length(ids))), ncol=2), others=0)
  
  mpg <- MPG(cost = cost_raster, patch = binary_mask)
  corridors <- st_as_sf(mpg@lcpLink)
  st_crs(corridors) <- 32651
  
  st_write(corridors %>% filter(lcpPerimWeight <= 100) %>% st_transform(4326), sprintf("D:/Data/Corridors_%s_100.geojson", name), driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  st_write(corridors %>% filter(lcpPerimWeight <= 500) %>% st_transform(4326), sprintf("D:/Data/Corridors_%s_500.geojson", name), driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  st_write(corridors %>% filter(lcpPerimWeight <= 1000) %>% st_transform(4326), sprintf("D:/Data/Corridors_%s_1000.geojson", name), driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
}


# ------------------------------------------------------------------------------
# PART 5: TIERED HYBRID EUCLIDEAN CONNECTIVITY (MEMORY-SAFE WATER INTEGRATION)
# ------------------------------------------------------------------------------
cat("\n--- PART 5: TIERED WATER-INTEGRATED CONNECTIVITY ---\n")
green_data_proj <- green_data %>% st_transform(32651) 
water_sf_proj <- st_read("D:/Data/NCR_Water.geojson", quiet = TRUE) %>% st_transform(32651)

tiered_scenario <- list(
  tier_filtered = green_data_proj %>% 
    filter(Patch_Area_sqm >= 1000) %>%
    filter(!(Patch_Area_sqm < 1500 & Core_Area_sqm == 0))
)

generate_corridors_with_water <- function(patches, max_dist, water_layer) {
  cat("    Finding neighbors within", max_dist, "meters...\n")
  neighbors <- st_is_within_distance(patches, patches, dist = max_dist)
  total_links <- sum(sapply(seq_along(neighbors), function(i) sum(neighbors[[i]] > i)))
  
  if (total_links == 0) return(st_sf(id1 = integer(), id2 = integer(), base_cost = numeric(), cost = numeric(), geometry = st_sfc(), crs = 32651))
  
  cat("    Drawing", total_links, "spatial links...\n")
  geom <- st_geometry(patches)
  ids <- patches$patch_id
  id1_vec <- integer(total_links)
  id2_vec <- integer(total_links)
  geom_list <- vector("list", total_links)
  
  pb1 <- txtProgressBar(min = 0, max = length(neighbors), style = 3, width = 50)
  counter <- 1
  for (i in seq_along(neighbors)) {
    targets <- neighbors[[i]][neighbors[[i]] > i] 
    for (j in targets) {
      id1_vec[counter] <- ids[i]
      id2_vec[counter] <- ids[j]
      geom_list[[counter]] <- st_nearest_points(geom[i], geom[j])
      counter <- counter + 1
    }
    setTxtProgressBar(pb1, i)
  }
  close(pb1)
  
  cat("\n    Compiling spatial network...\n")
  lines_sfc <- do.call(c, geom_list)
  lines_sf <- st_sf(id1 = id1_vec, id2 = id2_vec, base_cost = as.numeric(st_length(lines_sfc)), geometry = lines_sfc, crs = 32651)
  
  cat("    Evaluating strict water barriers (Memory-Safe Chunking)...\n")
  chunk_size <- 50000
  num_chunks <- ceiling(nrow(lines_sf) / chunk_size)
  crosses_water <- logical(nrow(lines_sf))
  
  pb2 <- txtProgressBar(min = 0, max = num_chunks, style = 3, width = 50)
  for (k in 1:num_chunks) {
    start_idx <- (k - 1) * chunk_size + 1
    end_idx <- min(k * chunk_size, nrow(lines_sf))
    chunk_lines <- lines_sf[start_idx:end_idx, ]
    intersect_list <- st_intersects(chunk_lines, water_layer, sparse = TRUE)
    crosses_water[start_idx:end_idx] <- lengths(intersect_list) > 0
    setTxtProgressBar(pb2, k)
  }
  close(pb2)
  
  lines_sf$cost <- lines_sf$base_cost 
  if (any(crosses_water)) {
    lines_sf$cost[crosses_water] <- lines_sf$base_cost[crosses_water] + 500
  }
  return(lines_sf)
}

for (name in names(tiered_scenario)) {
  cat(sprintf("\nProcessing scenario: %s...\n", name))
  scenario_patches <- tiered_scenario[[name]]
  all_corridors <- generate_corridors_with_water(scenario_patches, max_dist = 1000, water_layer = water_sf_proj)
  
  if (nrow(all_corridors) > 0) {
    st_write(all_corridors %>% filter(cost <= 100) %>% st_transform(4326), sprintf("D:/Data/Corridors_%s_100.geojson", name), driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
    st_write(all_corridors %>% filter(cost <= 500) %>% st_transform(4326), sprintf("D:/Data/Corridors_%s_500.geojson", name), driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
    st_write(all_corridors %>% filter(cost <= 1000) %>% st_transform(4326), sprintf("D:/Data/Corridors_%s_1000.geojson", name), driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  }
}


# ------------------------------------------------------------------------------
# PART 6: DATA DIAGNOSTICS 
# ------------------------------------------------------------------------------
cat("\n--- PART 6: DATA DIAGNOSTICS ---\n")
diag_patches <- readRDS("D:/Data/All_MetroManila_Patches.rds")
cat("CRS EPSG (Should be 4326):", st_crs(diag_patches)$epsg, "\n")
cat("Bounding Box:\n"); print(st_bbox(diag_patches))
cat("Data Structure Preview:\n"); print(head(diag_patches, 3))


# ------------------------------------------------------------------------------
# PART 7: INTERACTIVE SHINY DASHBOARD
# Note: Running this will launch the local web server.
# ------------------------------------------------------------------------------
cat("\n--- PART 7: LAUNCHING DASHBOARD ---\n")
web_patches <- readRDS("D:/Data/All_MetroManila_Patches.rds") %>% filter(total_pixel_count >= 10)

ui <- dashboardPage(
  dashboardHeader(title = "NCR Green Metrics"),
  dashboardSidebar(
    sidebarMenu(
      radioButtons("filter_type", "Show Patches:",
                   choices = c("All Significant Patches" = "all", 
                               "Only Core Areas" = "core", 
                               "Only Islets" = "islet"))
    )
  ),
  dashboardBody(
    leafletOutput("map", height = "80vh")
  )
)

server <- function(input, output, session) {
  filtered_data <- reactive({
    if (input$filter_type == "core") {
      web_patches %>% filter(has_core == TRUE)
    } else if (input$filter_type == "islet") {
      web_patches %>% filter(is_islet == TRUE)
    } else {
      web_patches
    }
  })
  
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>% 
      setView(lng = 121.03, lat = 14.6, zoom = 11)
  })
  
  observe({
    req(filtered_data()) 
    leafletProxy("map") %>%
      clearShapes() %>%
      addPolygons(
        data = filtered_data(), 
        color = "#27ae60", 
        weight = 1, 
        fillOpacity = 0.7,
        popup = ~paste("<b>Patch ID:</b>", patch_id, "<br><b>Pixel Count:</b>", total_pixel_count) 
      )
  })
}

shinyApp(ui, server)