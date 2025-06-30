import os
import json
import geopandas as gpd
from shapely.geometry import shape, mapping
from pystac_client import Client
from planetary_computer import sign
import stackstac
import rioxarray
import xarray as xr
from datetime import datetime
import numpy as np
from dataclasses import dataclass
from typing import Union
import traceback


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
    polarizations: list = None  # e.g., ["VH", "VV"]
    orbit_pass: Union[str, None] = None  # "ascending", "descending", or None
    resolution: float = 10  # in meters (will resample as necessary)


def load_roi(roi_geojson_path):
    """Load ROI geometry and convert to GeoDataFrame with EPSG:4326 CRS."""
    with open(roi_geojson_path) as f:
        roi_geojson = json.load(f)
    roi_geom = shape(roi_geojson["features"][0]["geometry"])
    gdf = gpd.GeoDataFrame.from_features(roi_geojson["features"])
    gdf.set_crs("EPSG:4326", inplace=True)
    return roi_geom, gdf


def get_stac_items(roi_geom, start_date, end_date, source):
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


def create_output_dirs(output_dir):
    """Create directories for MS (RGB) and SWIR outputs."""
    ms_dir = os.path.join(output_dir, "MS")
    swir_dir = os.path.join(output_dir, "SWIR")
    os.makedirs(ms_dir, exist_ok=True)
    os.makedirs(swir_dir, exist_ok=True)
    return ms_dir, swir_dir


def get_bands_by_source(source):
    """Return band names for RGB, SWIR, and SCL depending on STAC source.

    Args:
        source (str): STAC source, either "planetary" or "earth-search".
    Returns:
        dict: Dictionary with band names for MS, SWIR, and SCL.

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


def print_properties(catalog, collection_name="sentinel-1-grd"):
    try:
        collection = catalog.get_collection(collection_name)
        print(f"Collection '{collection.id}' description: {collection.description}\n")

        # Just get a sample item to inspect available properties
        search = catalog.search(collections=[collection_name], limit=1)
        item = next(search.get_items(), None)

        if item:
            print("Sample item properties:\n")
            for key in sorted(item.properties.keys()):
                print(f" - {key}")
        else:
            print("No items found in collection.")
    except Exception as e:
        print(f"Error retrieving collection or items: {e}")


def print_info_about_STAC(catalog, collection_name="sentinel-1-grd"):
    try:
        from pystac_client.conformance import ConformanceClasses

        conforms = catalog.conformance_classes

        print("Supported conformance classes:")
        for url in conforms:
            print(f" - {url}")

        if ConformanceClasses.FILTER in conforms:
            print("This STAC server supports FILTER.")
            print_queryable_parameters(catalog, collection_name)
        else:
            print("FILTER not supported. Cannot use get_queryables().")
            print_properties(catalog, collection_name)
    except Exception as e:
        print(f"Error retrieving STAC catalog information: {e}")
        print_properties(catalog, collection_name)


def print_queryable_parameters(catalog, collection_name="sentinel-1-grd"):
    collection = catalog.get_collection(collection_name)

    # Then get the queryables from the collection
    queryables = collection.get_queryables()
    import pprint

    pprint.pprint(queryables, indent=3)


def download_s1_imagery(settings: Union[DownloadSettings, dict]):
    """
    Downloads Sentinel-1 GRD imagery with specified settings.
    Args:
        settings (Union[DownloadSettings, dict]): Configuration including ROI, dates, polarizations, etc.
    """
    if isinstance(settings, dict):
        settings = DownloadSettings(**settings)

    roi_geom, gdf = load_roi(settings.roi_geojson_path)
    print(f"gdf.crs: {gdf.crs}")
    print(f"roi_geom: {roi_geom}")

    catalog_url = {
        "planetary": "https://planetarycomputer.microsoft.com/api/stac/v1",
        "earth-search": "https://earth-search.aws.element84.com/v1",
    }.get(settings.source)

    if catalog_url is None:
        raise ValueError(f"Unsupported source: {settings.source}")

    catalog = Client.open(catalog_url)
    # print_info_about_STAC(catalog, "sentinel-1-grd")
    print(f"catalog: {catalog}")
    query = {
        "sar:instrument_mode": {"eq": "IW"},  # just this works fine
        # "platform": {"eq": "sentinel-1"},
    }

    # I don't think either of these will work for earth search
    # if settings.polarizations:
    #     query["polarization"] = {"in": settings.polarizations} # This does not work for earht search
    # if settings.orbit_pass:
    #     query["sat:orbit_state"] = {"eq": settings.orbit_pass}
    # print_info_about_STAC(catalog, "sentinel-1-grd")
    # exit()
    search = catalog.search(
        collections=["sentinel-1-grd"],
        intersects=roi_geom,
        datetime=f"{settings.start_date}/{settings.end_date}",
        query=query,
    )

    items = [
        sign(item) if settings.source == "planetary" else item
        for item in search.items()
    ]
    print(f"Found {len(items)} Sentinel-1 items.")
    # exit()

    if not items:
        print("No Sentinel-1 items found.")
        return

    s1_dir = os.path.join(settings.output_dir, "S1")
    os.makedirs(s1_dir, exist_ok=True)

    for item in items:
        print(f"item.assets: {item.assets}")
        # exit()
        try:
            # note this is sensitive to the capitalization of the polarization names
            # For example for earth search it used vh and vv instead of VH and VV
            for pol in settings.polarizations or ["vh"]:
                # for pol in ["vh"]:
                # if pol not in item.assets:
                #     print(f"Skipping {item.id}: {pol} polarization not available.")
                #     continue

                # asset = item.assets[pol]
                asset = item.assets[pol]
                print(f"asset: {asset}")
                da = rioxarray.open_rasterio(asset.href, masked=True).squeeze()

                from rasterio.enums import Resampling

                target_crs = gdf.crs
                da_reprojected = da.rio.reproject(
                    dst_crs=target_crs,
                    # You can specify resampling method if you like:
                    resampling=Resampling.cubic,  # or bilinear/cubic
                )

                print(f"da open rasterio: {da}")
                print("Raster EPSG:", da_reprojected.rio.crs.to_epsg())
                print("GeoDataFrame EPSG:", gdf.crs.to_epsg())

                # Clip
                clipped_ds = da_reprojected.rio.clip(
                    gdf.geometry,
                    crs=da_reprojected.rio.crs,
                    drop=True,
                    all_touched=True,
                )
                output_path = os.path.join(s1_dir, f"{item.id}_{pol}.tif")
                clipped_ds.rio.to_raster(output_path)

                # we need to reproject to it to the same CRS as the ROI GeoDataFrame
                # commented out for S1
                # da = da.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs(
                #     "EPSG:4326"
                # )

                print(f"da after setting spatial dims and CRS: {da}")
                # da_clipped = da.rio.clip(gdf.geometry.apply(mapping), gdf.crs)
                # print(f"da_clipped: {da_clipped}")
                # if settings.resolution:
                #     da_clipped = da_clipped.rio.reproject(
                #         "EPSG:4326", resolution=settings.resolution
                #     )
                # print(f"da_clipped after reprojection: {da_clipped}")

                # output_path = os.path.join(s1_dir, f"{item.id}_{pol}.tif")
                # da.rio.to_raster(output_path)
                # da_clipped.rio.to_raster(output_path)
                print(f"Saved {pol} polarization image for {item.id} at {output_path}")
                exit()

                print(f"Saved S1 {pol} image for {item.id}")

        except Exception as e:
            print(f"Failed to process {item.id}: {e}")
            print(traceback.format_exc())


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

    # Get bands based on the source each source has different band names
    # e.g. for planetary computer: ms_bands = ["B04", "B03", "B02"]  # RGB
    bands = get_bands_by_source(settings.source)
    ms_bands = bands["ms_bands"]
    swir_band = bands["swir_band"]
    scl_band = bands["scl_band"]

    print(f"Found {len(items)} items in the STAC catalog.")
    if not items:
        print("No items found in the specified date range and ROI.")
        return

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


# Example usage: Search and download Sentinel-2 imagery from the earth-search STAC catalog.
# start_date = "2023-07-01"
# end_date = "2023-07-30"
# STAC_source = "earth-search"  # or options are: "earth-search", "planetary"
# settings = DownloadSettings(
#     roi_geojson_path=os.path.join(os.getcwd(), "examples", "rois.geojson"),
#     start_date=start_date,
#     end_date=end_date,
#     max_cloud_cover=50,
#     output_dir=f"s2_from_{STAC_source}_{start_date}_{end_date}",
#     source=STAC_source,
# )
# download_s2_imagery(settings)

# Example usage: Search and download Sentinel-2 imagery from the earth-search STAC catalog.
# start_date = "2020-07-01"
# end_date = "2023-07-30"
# STAC_source = "planetary"  # or options are: "earth-search", "planetary"
# settings = DownloadSettings(
#     roi_geojson_path=os.path.join(os.getcwd(), "examples", "rois.geojson"),
#     start_date=start_date,
#     end_date=end_date,
#     max_cloud_cover=50,
#     output_dir=f"s2_from_{STAC_source}_{start_date}_{end_date}",
#     source=STAC_source,
# )
# download_s2_imagery(settings)

# Example usage: Search and download Sentinel-1 imagery from the planetary computer STAC catalog.
start_date = "2021-06-01"
end_date = "2021-06-15"
settings = DownloadSettings(
    roi_geojson_path=os.path.join(os.getcwd(), "examples", "rois.geojson"),
    start_date=start_date,
    end_date=end_date,
    max_cloud_cover=100,  # Not used for S1 but required by class
    output_dir="S1_data_output_attempt_2",
    source="planetary",  # or "planetary",earth-search"
    polarizations=["vh", "vv"],
    orbit_pass="None",  # or "descending", or None for both
    resolution=10,  # meters
)
download_s1_imagery(settings)

# -----------------------------------------
# verdict for S1:
# Earth search does not work for S1 because:
# - cannot download from earth-search for free
# - the earth search catalog is a static catelogue and does NOT allow for dynamic queries
# - it uses polarizations vh and vv instead of VH and VV
# ----- ----- ----  ----- ----- -----
# Planetary computer works fine for S1:
# - search works!
# download seems to hang after printing the assets
# item.assets: {'vh': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/measurement/iw-vh.tiff?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'vv': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/measurement/iw-vv.tiff?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'thumbnail': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/preview/quick-look.png?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'safe-manifest': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/manifest.safe?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'schema-noise-vh': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/annotation/calibration/noise-iw-vh.xml?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'schema-noise-vv': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/annotation/calibration/noise-iw-vv.xml?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'schema-product-vh': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/annotation/iw-vh.xml?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'schema-product-vv': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/annotation/iw-vv.xml?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'schema-calibration-vh': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/annotation/calibration/calibration-iw-vh.xml?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'schema-calibration-vv': <Asset href=https://sentinel1euwest.blob.core.windows.net/s1-grd/GRD/2021/6/15/IW/DV/S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB_206A/annotation/calibration/calibration-iw-vv.xml?st=2025-06-11T15%3A43%3A48Z&se=2025-06-12T16%3A28%3A48Z&sp=rl&sv=2024-05-04&sr=c&skoid=9c8ff44a-6a2c-4dfb-b298-1c9212f64d9a&sktid=72f988bf-86f1-41af-91ab-2d7cd011db47&skt=2025-06-12T04%3A24%3A45Z&ske=2025-06-19T04%3A24%3A45Z&sks=b&skv=2024-05-04&sig=Kp5QJT8uVxUkf7xZ6%2Bhkxv5XlN3iUBPHgH2ZuWs5tmM%3D>, 'tilejson': <Asset href=https://planetarycomputer.microsoft.com/api/data/v1/item/tilejson.json?collection=sentinel-1-grd&item=S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB&assets=vv&assets=vh&expression=vv%3Bvh%3Bvv%2Fvh&rescale=0%2C600&rescale=0%2C270&rescale=0%2C9&asset_as_band=True&tile_format=png&format=png>, 'rendered_preview': <Asset href=https://planetarycomputer.microsoft.com/api/data/v1/item/preview.png?collection=sentinel-1-grd&item=S1B_IW_GRDH_1SDV_20210615T140744_20210615T140809_027368_0344CB&assets=vv&assets=vh&expression=vv%3Bvh%3Bvv%2Fvh&rescale=0%2C600&rescale=0%2C270&rescale=0%2C9&asset_as_band=True&tile_format=png&format=png>}
# -----------------------------------------
