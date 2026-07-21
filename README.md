# NCR Urban Connectivity: Metro Manila Green Spaces
Mapping and analyzing the landscape connectivity of urban green spaces in Metro Manila using spatial network models and least-cost path analysis in R.

## Project Overview
This project extracts, classifies, and models the structural connectivity of urban green spaces across the National Capital Region (NCR) of the Philippines. By utilizing high-resolution EVI (Enhanced Vegetation Index) multi-class rasters, the pipeline identifies ecologically significant habitat patches and routes dispersal corridors while accounting for major urban barriers like concrete infrastructure and water bodies.

## Core Methodology
1. **Spatial Extraction & Classification:** Isolates dense and sparse vegetation from multi-class GeoTIFFs using `terra`, applying Rook's (4-way) contiguity to delineate distinct habitat patches.
2. **Habitat Metrics:** Calculates core area, edge area, and total physical patch size (in square meters) to evaluate the ecological value of over 77,000 distinct patches.
3. **Least-Cost Path (LCP) Routing:** Utilizes the `grainscape` package to model ideal dispersal networks across the urban matrix, using weighted cost surfaces to penalize movement through concrete and sparse vegetation.
4. **Memory-Safe Euclidean Barrier Modeling:** Implements a custom, chunk-based spatial intersection algorithm via the `sf` package to safely process hundreds of thousands of distance vectors, applying strict dispersal penalties (+500m) to any ecological corridor crossing a water body.
5. **Interactive Visualization:** Packages the filtered, mathematically projected (WGS84) dataset into a lightweight `.rds` format for deployment on an interactive Shiny/Leaflet web dashboard.

## Built With
* **R**
* `sf` (Spatial vector handling)
* `terra` (High-performance raster processing)
* `grainscape` (Landscape connectivity modeling)
* `Shiny` & `Leaflet` (Interactive web deployment)
