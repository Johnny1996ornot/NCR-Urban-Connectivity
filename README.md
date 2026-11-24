# NCR Urban Connectivity: Metro Manila Green Spaces

### 👉 [**Open the Live Interactive Web Map**](https://johnny1996ornot.github.io/NCR-Urban-Connectivity/)

Interactive spatial network visualization of dense urban vegetation patches and multi-threshold functional dispersal corridors across the 17 local government units (LGUs) of Metro Manila, Philippines.

---

## Network Connectivity Summary

| Dispersal Threshold | Maximum Link Weight (30 m pixels) | Total Functional Links | Total Network Components | Ecological Condition |
| :--- | :---: | :---: | :---: | :--- |
| **1,000 m Threshold** | 33.33 | 11,574 | 15 | Highly connected |
| **400 m Threshold** | 13.33 | 9,616 | 215 | Fragmented |
| **100 m Threshold** | 3.33 | 3,420 | 2,290 | Severely fragmented |

---

## Map Features & Workflow

1. **Land Cover & Patch Extraction:** Dense vegetation patches extracted from Landsat 8/9 imagery via supervised Random Forest classification in Google Earth Engine and analyzed at 30 m spatial resolution.
2. **Spatial Graph Modeling:** Functional patch connectivity evaluated across `100 m`, `400 m`, and `1,000 m` avian/wildlife dispersal thresholds using `grainscape` in R and exported to WGS 84 (`EPSG:4326`) GeoJSON via PyQGIS.
3. **Interactive Web Visualization:** Built with Leaflet.js Canvas rendering featuring label-free light street and satellite basemap toggles, threshold switching, and click-to-inspect patch and corridor metrics.
