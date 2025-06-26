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

### Sentinel-2 Downloader

Download Sentinel-2 L2A imagery (RGB and SWIR) from Microsoft Planetary Computer or AWS Earth Search.

  #### Features
  
  - ROI-based search using GeoJSON
  - Filter by date and cloud cover
  - Saves RGB and SWIR bands as GeoTIFFs
    - Saves the 3 band RGB tiffs to `output_folder/MS`  
    - Saves the 1 band SWIR tiffs to `output_folder/SWIR`


## 🚀 Future Work
  - Add Sentinel‑1 download support
  
  - Add Landsat integration
  
## Contribution
Contributions are welcome! Open an issue or PR for bug fixes, feature requests, or enhancements.
