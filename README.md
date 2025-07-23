# STAC_downloader
A set of scripts used to download STAC data for free. This project is still under development. The goal of this project to make a series of easy to use scripts to download data from STAC providers. 

## Overview
STAC_downloader is a command-line tool for researchers, GIS specialists, and Earth observation enthusiasts to easily download **free satellite imagery** worldwide.  
- ✅ **Sentinel‑2** (optical imagery) downloads are fully supported
- ✅ **Landsat (optical imagery)** downloads are fully supported  
  - Currently supports **Microsoft Planetary Computer** only 
- 🚧 **Sentinel‑1** (radar imagery) download functionality is planned for upcoming releases

### Core Concept
Ideally the result of this project will be a single easy to use script that will let users pick from a list of STAC providers, choose any of the supported catalogues and download tiffs for free.

<img width="808" height="598" alt="STAC selection drawio" src="https://github.com/user-attachments/assets/d6813431-a70c-4992-90e4-0bea01293f90" />


## Installation
Install STAC_downloader via [pixi](https://pixi.sh/v0.47.0/?utm_source=chatgpt.com), a fast, reproducible package manager:

💡 You only need to install the `pixi` CLI if it's not already on your system. See [pixi.sh](https://pixi.sh/) for installation instructions.
   ```bash
   curl -fsSL https://pixi.sh/install.sh | sh
  ```
1. From your STAC_downloader project directory (with pixi.toml or pyproject.toml), run:
  ```bash
  pixi install
  pixi shell
  ```
  ▶️ Make sure both commands are executed inside the main project repo.
  - `pixi install`: Installs all the project dependencies defined in the `pixi.toml` file.
  - `pixi shell`: Activates a virtual environment with all dependencies available.

## Available Scripts

| Script name                   | Description                                                                                               |
| ----------------------------- | --------------------------------------------------------------------------------------------------------- |
| `download_S2_STAC_Imagery.py` | Downloads Sentinel-2 L2A imagery (RGB and SWIR) from Microsoft Planetary Computer or AWS Earth Search     |
| `download_landsat_STAC.py`    | Downloads Landsat surface reflectance imagery from STAC endpoints ( Microsoft PC) |


### 📥 download_STAC_S2_imagery.py: Download Sentinel-2 Imagery from STAC APIs
This script downloads Sentinel-2 L2A imagery (RGB and SWIR bands) from STAC-compliant APIs such as Microsoft Planetary Computer or AWS Earth Search.

It supports:
   
   - ROI-based filtering (GeoJSON)
   - Date range selection
   - Cloud cover filtering via the Scene Classification Layer (SCL)
   - Output as clipped GeoTIFFs for RGB (MS/) and SWIR (SWIR/) bands

#### ✅ Example Usage

```bash
python download_STAC_S2_imagery.py --roi examples/rois.geojson --start-date 2022-06-01 --end-date 2022-06-30 --cloud-cover 30 --output-dir output_june_2022 --source earth-search
```
All outputs are clipped to your ROI and cloud-masked.
Example Output
```bash
PS C:\STAC_downloads> python download_STAC_S2_imagery.py --roi examples/rois.geojson --start-date 2022-06-01 --end-date 2022-06-30 --cloud-cover 30 --output-dir output_june_2022 --source earth-search
   Found 11 items in the STAC catalog. Do you want to proceed with downloading? (y/n): y
   Saved imagery for S2B_10SEF_20220628_0_L2A
   Saved imagery for S2B_10SFF_20220628_0_L2A
   Skipping S2A_10SFF_20220623_0_L2A due to cloud cover 38.07%
```

#### 🔧 Arguments
| Flag            | Description                                                              |
|-----------------|---------------------------------------------------------------------------|
| `--roi`         | Path to ROI GeoJSON file (**required**)                                   |
| `--start-date`  | Start date in format `YYYY-MM-DD` (**required**)                          |
| `--end-date`    | End date in format `YYYY-MM-DD` (**required**)                            |
| `--cloud-cover` | Max allowed cloud cover (%) — default: `50`                               |
| `--output-dir`  | Output directory — default: `output/`                                     |
| `--source`      | STAC provider: `planetary` or `earth-search` — default: `planetary`       |


#### 📁 Output Structure
All outputs are clipped to your ROI and cloud-masked.
 - output_dir/MS/ → RGB imagery
- output_dir/SWIR/ → SWIR imagery
```
    output/
    ├── MS/
    │   ├── <item_id>_MS.tif
    └── SWIR/
        ├── <item_id>_SWIR.tif
```
Example Image: The multispectral tiff (aka RGB) `S2B_10SEF_20220618_1_L2A_MS.tif` rendered in QGIS with the ROI geometry in red. 

![image](https://github.com/user-attachments/assets/590f2378-91ca-4aae-953d-4d5de811a065)




### 📥 download_landsat_STAC.py.py: Download Landsat Imagery from STAC APIs
This script downloads Landsat from the STAC-compliant API Microsoft Planetary Computer.

It supports:
   - ROI-based filtering (GeoJSON)
   - Date range selection
   - Cloud cover filtering via the Scene Classification Layer (SCL)
   - Output as clipped GeoTIFFs with all the bands saved into a single tiff. See example below
<img width="555" height="496" alt="all_in_one_tiff_diagran" src="https://github.com/user-attachments/assets/71abcc1c-b20f-44e6-80e0-b2be08bb933c" />




#### ✅ Example Usage
### Parameters
The examples below will showcase how to use these parameters to effectively search for and download the landsat data from the Planetary computer endpoint.
Note: Downloads may be slow or metered as this is a free endpoint.
| Param                   | What it does                                                      | Example                      |
| ----------------------- | ----------------------------------------------------------------- | ---------------------------- |
| `start_date`,`end_date` | ISO strings (`YYYY-MM-DD`) defining the time window               | `"2025-05-01"`               |
| `roi_geojson`           | Path to polygon(s) in EPSG:4326                                   | `"examples/rois.geojson"`    |
| `missions`              | Landsat mission numbers; confirm availability in the STAC catalog | `[8, 9]`                     |
| `collection_id`         | STAC collection ID (`landsat-c2-l1`, `landsat-c2-l2`, etc.)       | `"landsat-c2-l2"`            |
| `max_cloud`             | % cloud over the entire scene (0–100)                             | `90`                         |
| `bands`                 | List of band names; use helper or custom list                     | `default_landsat_bands(...)` |
| `interactive_mode`      | Prompt before downloading                                         | `True` / `False`             |


### Example #1
This example shows how to download data from Landsat 8 and 9 for the landsat collection 2 level 2 (surface reflectance).
- You can read more about collection landsat-c2-l2 on planetary computer's [website](https://planetarycomputer.microsoft.com/dataset/landsat-c2-l1)
- Before the download begins you will be shown the number of available images and asked if you want to proceed, enter 'y' to start the download process.

<img width="651" height="121" alt="image" src="https://github.com/user-attachments/assets/f8f8603c-0a59-4fe2-bd01-e73490d91660" />

```python
if __name__ == "__main__":
    # 1. Choose collection & missions (check the STAC catalog first)
    default_bands = default_landsat_bands(collection_id="landsat-c2-l2")

    config = LandsatConfig(
        missions=[8, 9],               # e.g. [5], [7], or [8, 9]
        provider="planetary",          # currently only "planetary"
        collection_id="landsat-c2-l2", # Level-2 SR. Use "landsat-c2-l1" for Level-1
        max_cloud=90,                  # % over the FULL SCENE tile (not clipped to ROI)
        bands=default_bands,
        interactive_mode=True,         # ask before download
    )

    # 2. Define ROI and date range
    roi_geojson = os.path.join(os.getcwd(), "examples", "rois.geojson")
    start_date  = "2025-05-01"
    end_date    = "2025-05-30"

    # 3. Download
    download_landsat_with_config(
        config=config,
        roi_geojson=roi_geojson,
        start_date=start_date,
        end_date=end_date,
        dest_dir="output_landsat8_9_pc",
    )


```

### Example #2
```
if __name__ == "__main__":

    default_bands = default_landsat_bands(collection_id="landsat-c2-l1")

    # Example usage with LandsatConfig
    config_basic = LandsatConfig(
        missions=[5],
        provider="planetary",  # or options are: "planetary"
        collection_id="landsat-c2-l1",  # Surface reflectance
        max_cloud=90,
        bands=default_bands,
        interactive_mode=True,  # Set this to True to be asked to confirm the download before starting
    )

    # Usage examples:
    roi_geojson_path = os.path.join(os.getcwd(), "examples", "rois.geojson")
    start_date = "2012-06-01"
    end_date = "2012-07-30"

    # Download with config
    download_landsat_with_config(
        config=config_basic,
        roi_geojson=roi_geojson_path,
        start_date=start_date,
        end_date=end_date,
        dest_dir="output_landsat_5_pc",  # Change to your desired output directory
    )


```

#### Output Format

## 🚀 Future Work
  - Add Sentinel‑1 download support
  - Combined all the download scripts into a single script
  - Improve cloud filtering for landsat
  - Add support for more STAC providers
  
## Contribution
Contributions are welcome! Open an issue or PR for bug fixes, feature requests, or enhancements.
