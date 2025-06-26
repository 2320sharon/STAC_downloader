import os
import json
import geopandas as gpd
from typing import Tuple, List, Dict, Union
from shapely.geometry.base import BaseGeometry
from shapely.geometry import shape, mapping
from pystac_client import Client
from planetary_computer import sign
import stackstac
import rioxarray
from dataclasses import dataclass
import traceback
from pystac import Item

"""
Download Sentinel-2 L2A imagery (RGB + SWIR) from Planetary Computer or Earth Search.

Features:
- Load ROI from GeoJSON.
- Query images by date, location, and cloud cover.
- Clip to ROI and filter by cloud-free pixels using SCL.
- Save RGB (MS) and SWIR bands as GeoTIFFs.

Usage:

    Provide the script with the following parameters:
    1. `--roi`: Path to the ROI GeoJSON file.
    2. `--start-date`: Start date in 'YYYY-MM-DD' format.
    3. `--end-date`: End date in 'YYYY-MM-DD' format.
    4. `--cloud-cover`: Maximum cloud cover percentage (default: 50).
    5. `--output-dir`: Directory to save downloaded imagery (default: 'output').

    python download_STAC_S2_imagery.py --roi examples/rois.geojson --start-date 2022-06-01 --end-date 2022-06-30 --cloud-cover 30 --output-dir output_june_2022 --source earth-search

    Output Structure:
    output/
    ├── MS/
    │   ├── <item_id>_MS.tif
    └── SWIR/
        ├── <item_id>_SWIR.tif

        
Example at bottom of script.
"""


@dataclass
class DownloadSettings:
    """
    Configuration settings for downloading Sentinel-2 imagery.

    Attributes:
        roi_geojson_path (str): Path to ROI GeoJSON file.
        start_date (str): Start date in 'YYYY-MM-DD' format.
        end_date (str): End date in 'YYYY-MM-DD' format.
        max_cloud_cover (float): Maximum cloud cover allowed in percent.
        output_dir (str): Output directory for downloaded images.
        source (str): STAC source: "planetary" or "earth-search".
        additional_options (dict, optional): Additional parameters (currently unused).
    """

    roi_geojson_path: str
    start_date: str
    end_date: str
    max_cloud_cover: float
    output_dir: str = "output"
    source: str = "planetary"
    additional_options: dict = None


def load_roi(roi_geojson_path: str) -> Tuple[BaseGeometry, gpd.GeoDataFrame]:
    """
    Load a Region of Interest (ROI) from a GeoJSON file.

    Args:
        roi_geojson_path (str): Path to the GeoJSON file containing ROI geometry.

    Returns:
        Tuple[BaseGeometry, gpd.GeoDataFrame]: A shapely geometry (first feature) and a GeoDataFrame
        containing all geometries, both in EPSG:4326.
    """
    with open(roi_geojson_path) as f:
        roi_geojson = json.load(f)

    roi_geom = shape(roi_geojson["features"][0]["geometry"])
    gdf = gpd.GeoDataFrame.from_features(roi_geojson["features"])
    gdf.set_crs("EPSG:4326", inplace=True)
    return roi_geom, gdf


def get_stac_items(
    roi_geom: BaseGeometry, start_date: str, end_date: str, source: str
) -> List[Item]:
    """
    Query a STAC API for Sentinel-2 items that intersect the ROI and meet the date and cloud cover criteria.

    Args:
        roi_geom (BaseGeometry): Geometry representing the ROI.
        start_date (str): Start date in 'YYYY-MM-DD' format.
        end_date (str): End date in 'YYYY-MM-DD' format.
        source (str): Either "planetary" (Microsoft Planetary Computer) or "earth-search" (AWS).

    Returns:
        List[Item]: List of signed STAC Items matching the query.
    """
    """Search STAC API for Sentinel-2 items from the specified provider."""
    if source == "planetary":
        catalog = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
        search = catalog.search(
            collections=["sentinel-2-l2a"],
            intersects=roi_geom,
            datetime=f"{start_date}/{end_date}",
            query={"eo:cloud_cover": {"lt": 80}},
        )
        items = [sign(item) for item in search.items()]
    elif source == "earth-search":
        catalog = Client.open("https://earth-search.aws.element84.com/v1")
        search = catalog.search(
            collections=["sentinel-2-l2a"],
            intersects=roi_geom,
            datetime=f"{start_date}/{end_date}",
            query={"eo:cloud_cover": {"lt": 80}},
        )
        items = list(search.items())
    else:
        raise ValueError(f"Unsupported source: {source}")
    return items


def create_output_dirs(output_dir: str) -> Tuple[str, str]:
    """
    Create output subdirectories for RGB (MS) and SWIR bands.

    Args:
        output_dir (str): Base directory to store downloaded data.

    Returns:
        Tuple[str, str]: Paths to MS and SWIR output folders.
    """
    ms_dir = os.path.join(output_dir, "MS")
    swir_dir = os.path.join(output_dir, "SWIR")
    os.makedirs(ms_dir, exist_ok=True)
    os.makedirs(swir_dir, exist_ok=True)
    return ms_dir, swir_dir


def get_bands_by_source(source: str) -> Dict[str, List[str]]:
    """
    Get the appropriate band names for RGB, SWIR, and Scene Classification Layer (SCL) based on the STAC data source.

    Args:
         source (str): STAC source ("planetary" or "earth-search").
    Returns:
        Dict[str, Union[List[str], str]]: Band names grouped as ms_bands, swir_band, and scl_band.

    Information on the earth search colleciton:
       https://earth-search.aws.element84.com/v1/collections/sentinel-2-l2a
    information on the planetary computer colleciton:
       https://planetarycomputer.microsoft.com/api/stac/v1/collections/sentinel-2-l2a

    """
    if source == "planetary":
        return {
            "ms_bands": ["B04", "B03", "B02"],  # RGB bands
            "scl_band": "SCL",  # Scene Classification Layer
            "swir_band": "B11",  # SWIR band
        }
    elif source == "earth-search":
        return {
            "ms_bands": ["red", "blue", "green"],
            "scl_band": "scl",
            "swir_band": "swir16",
        }
    else:
        raise ValueError(f"Unsupported source: {source}")


def download_s2_imagery(settings: Union[DownloadSettings, dict]):
    """
    Main function to download and process Sentinel-2 imagery from STAC catalogs.

    Args:
        settings (Union[DownloadSettings, dict]): Configuration settings.
    """
    if isinstance(settings, dict):
        settings = DownloadSettings(**settings)

    roi_geom, gdf = load_roi(settings.roi_geojson_path)

    # Create output directories
    ms_dir, swir_dir = create_output_dirs(settings.output_dir)
    # Get STAC items from selected catalog
    items = get_stac_items(
        roi_geom, settings.start_date, settings.end_date, settings.source
    )

    if not items:
        print(
            "No items found in the specified date range and ROI. Please adjust your query. Exiting."
        )
        return

    # Ask the user if they want to proceed with the download
    proceed = input(
        f"Found {len(items)} items in the STAC catalog. Do you want to proceed with downloading? (y/n): "
    )

    if proceed.lower() != "y":
        print("Download canceled.")
        return

    # Get bands based on the source each source has different band names
    # e.g. for planetary computer: ms_bands = ["B04", "B03", "B02"]  # RGB
    bands = get_bands_by_source(settings.source)
    ms_bands = bands["ms_bands"]
    swir_band = bands["swir_band"]
    scl_band = bands["scl_band"]

    # For each item returned by the STAC search check if all required assets exist then download each item
    for item in items:
        try:
            # Ensure all required assets exist
            required_assets = ms_bands + [swir_band, scl_band]
            missing_assets = [
                band for band in required_assets if band not in item.assets
            ]
            if missing_assets:
                print(f"Skipping {item.id}: missing assets {missing_assets}")
                continue

            # Load RGB and SWIR bands using stackstac
            # Note if you provide epsg=4326, it will reproject the data to EPSG:4326 and if you specify resolution= 10 it will be relative to the EPSG:4326 CRS aka NOT 10m resolution but in degrees
            data = stackstac.stack(
                [item],
                assets=ms_bands + [swir_band],
                bounds=roi_geom.bounds,
                epsg=4326,
                chunksize=2048,
            )

            if data.shape[0] == 0:
                print(f"Skipping {item.id}: no data returned by stackstac.")
                continue

            data = data.squeeze("time")

            # convert CRS to EPSG:4326 for consistency
            data = data.rio.write_crs("EPSG:4326")

            # Clip to exact ROI geometry
            data = data.rio.clip(gdf.geometry.apply(mapping), gdf.crs)

            # Mask out clouds using SCL (Scene Classification Layer)
            scl = rioxarray.open_rasterio(
                item.assets[scl_band].href, masked=True
            ).squeeze()
            scl_clipped = scl.rio.clip(gdf.geometry.apply(mapping), gdf.crs)

            clear_mask = scl_clipped.isin([4, 5, 6, 7])  # vegetation, bare, water, etc.
            clear_pixels = clear_mask.sum().item()
            total_pixels = scl_clipped.size
            clear_ratio = clear_pixels / total_pixels

            if clear_ratio < (1 - settings.max_cloud_cover / 100):
                print(
                    f"Skipping {item.id} due to cloud cover {100 - clear_ratio * 100:.2f}%"
                )
                continue

            # Save MS (RGB)
            rgb = data.sel(band=ms_bands)
            rgb_path = os.path.join(ms_dir, f"{item.id}_MS.tif")
            rgb.rio.to_raster(rgb_path)

            # Save SWIR
            swir = data.sel(band=swir_band)
            swir_path = os.path.join(swir_dir, f"{item.id}_SWIR.tif")
            swir.rio.to_raster(swir_path)

            print(f"Saved imagery for {item.id}")

        except Exception as e:
            print(f"Failed to process {item.id}: {e}")
            print(traceback.format_exc())


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Download Sentinel-2 imagery (RGB & SWIR) from Planetary Computer or Earth Search."
    )

    parser.add_argument(
        "--roi", type=str, required=True, help="Path to the ROI GeoJSON file."
    )
    parser.add_argument(
        "--start-date", type=str, required=True, help="Start date (YYYY-MM-DD)."
    )
    parser.add_argument(
        "--end-date", type=str, required=True, help="End date (YYYY-MM-DD)."
    )
    parser.add_argument(
        "--cloud-cover",
        type=float,
        default=50.0,
        help="Maximum cloud cover percentage (default: 50).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output",
        help="Directory to save downloaded imagery (default: 'output').",
    )
    parser.add_argument(
        "--source",
        type=str,
        choices=["planetary", "earth-search"],
        default="planetary",
        help="STAC API source: 'planetary' (Microsoft) or 'earth-search' (AWS). Default is 'planetary'.",
    )

    args = parser.parse_args()

    settings = DownloadSettings(
        roi_geojson_path=args.roi,
        start_date=args.start_date,
        end_date=args.end_date,
        max_cloud_cover=args.cloud_cover,
        output_dir=args.output_dir,
        source=args.source,
    )

    download_s2_imagery(settings)
