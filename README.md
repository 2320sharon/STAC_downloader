# STAC_downloader
A set of scripts used to download STAC data for free

## Overview
STAC_downloader is a command-line tool for researchers, GIS specialists, and Earth observation enthusiasts to easily download **free satellite imagery** worldwide.  
- ✅ **Sentinel‑2** (optical imagery) downloads are fully supported  
- 🚧 **Sentinel‑1** (radar imagery) and **Landsat** download functionality are planned for upcoming releases

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

| Script name                   | Description                                         |
|------------------------------|-----------------------------------------------------|
| `download_S2_STAC_Imagery.py`     |   Downloads Sentinel-2 L2A imagery (RGB and SWIR) from Microsoft Planetary Computer or AWS Earth Search                                               |

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




## 🚀 Future Work
  - Add Sentinel‑1 download support
  
  - Add Landsat integration
  
## Contribution
Contributions are welcome! Open an issue or PR for bug fixes, feature requests, or enhancements.
