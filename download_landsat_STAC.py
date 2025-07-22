import os
from typing import Iterable, Mapping, Sequence, Tuple, List, Any
import warnings
import xarray as xr
import geopandas as gpd
from shapely.geometry import shape
import planetary_computer
import pystac_client
import rioxarray  # keep this import to perform raster clip operations
from rasterio.errors import RasterioIOError
import stackstac
from typing import Union, Sequence
from shapely.geometry import mapping
from dataclasses import dataclass, field
from collections.abc import Iterable
import tqdm
import warnings

# this is to ensure that warnings are raised as errors
# @dev only
# warnings.filterwarnings("error", category=RuntimeWarning)


@dataclass
class LandsatConfig:
    missions: Union[int, Sequence[int]]
    provider: str = "pc"
    collection_id: str = "landsat-c2-l2"
    max_cloud: float = 90
    bands: Sequence[str] = field(default_factory=list)
    interactive_mode: bool = False


# Helper: common routine for searching Landsat STAC (Collection 2)
def search_landsat_collection(
    missions: Union[int, Sequence[int]],
    provider: str,
    roi_geom: dict,
    start_date: str,
    end_date: str,
    max_cloud: float,
    collection_id: str = "landsat-c2-l2",  # default to Collection 2 Level-2 Landsat data
):
    """
    Internal helper to search the Landsat Collection 2 STAC for a given mission.
    Returns the STAC Item with lowest cloud cover within criteria, or None.

    Parameters:
        missions: Single Landsat mission number (5, 7, 8, 9) or a sequence of missions.
        provider: Name of STAC provider (e.g., 'pc' for Planetary Computer)
            Note: Only planetary computer is supported currently.
        roi_geom: GeoJSON-like dict defining the ROI polygon geometry that the search will intersect.
        Results filtered to only those intersecting the geometry.
        start_date : Date to start searching for imagery (YYYY-MM-DD).
        end_date: Date to stop searching for imagery (YYYY-MM-DD).
        max_cloud: Maximum acceptable cloud cover percentage (0-100).
        collection_id: STAC collection ID to search within (default is Collection 2 Level-2).

    Returns:
        pystac_client.item_search.ItemSearch: Search result containing matching Landsat items.

    Raises:
        ValueError: If no items are found matching the search criteria.
        ValueError: If the provider does not have the specified collection.
        Exception: If an error occurs while searching the catalog.

    """

    # Validate that the requested provider has the specified collection
    if not validate_provider_has_collection(provider, collection_id):
        raise ValueError(
            f"Provider {provider} does not have collection {collection_id}"
        )
    # Open the STAC API catalog
    catalog = get_provider_catalog(provider)
    print(f"Searching Landsat {', '.join(map(str, missions))} imagery from {catalog}")

    # Build date-time range string
    date_range = (
        f"{start_date}/{end_date}" if end_date else f"{start_date}/{start_date}"
    )

    # Prepare query filters. Use STAC "query" to filter by platform (satellite) and cloud cover.
    # 'eo:cloud_cover' is a common metadata field for cloud percentage:contentReference[oaicite:9]{index=9}.
    query = {
        "platform": create_landsat_platform_filter(missions),
        "eo:cloud_cover": {"lte": max_cloud},
    }

    print(f"Searching for items within date range: {date_range}")

    # Perform the search (geometry filter via intersects for exact ROI, or bbox)
    search = catalog.search(
        collections=[collection_id],
        intersects=roi_geom,
        datetime=date_range,
        query=query,
    )

    # Check if any items were found
    if not search:
        print("No items found.")
        raise ValueError("No items found matching criteria")

    return search


def get_items_from_search(search: pystac_client.item_search.ItemSearch):
    """
    Retrieves items from a STAC ItemSearch object.

    Args:
        search (pystac_client.item_search.ItemSearch): The STAC ItemSearch object to retrieve items from.

    Returns:
        pystac.item_collection.ItemCollection : Collection of items found in the search.

    Raises:
        ValueError: If no items are found matching the search criteria.
    """
    items = search.item_collection()
    print(f"Found {len(items)} items")
    if not items:
        raise ValueError("No imagery found matching criteria")
    return items


def determine_provider(provider: str) -> str:
    """
    Determine the provider string for STAC search.
    Accepts variations of 'pc' for Planetary Computer or 'earth-search' for AWS Earth Search.
    Note: Only Planetary Computer is supported currently.
    """
    if provider.lower() in [
        "pc",
        "planetary",
        "planetarycomputer",
        "planetary-computer",
    ]:
        return "pc"
    elif provider.lower() in [
        "aws",
        "earth-search",
    ]:  # not available currently because landsat is NOT free to download
        return "earth-search"
    else:
        raise ValueError(
            f"Unknown provider '{provider}'. Use 'pc' for Planetary Computer or 'earth-search' for AWS."
        )


def try_get_item_assets(
    catalog: Any, collection_id: str
) -> Iterable[str] | Mapping[str, Any] | str | None:
    """Return item_assets or None if missing anywhere."""
    try:
        collection = catalog.get_child(collection_id)
        return getattr(collection, "item_assets", None)
    except AttributeError:
        return None


def filter_assets(
    requested: Sequence[str],
    item_assets: Iterable[str] | Mapping[str, Any] | str | None,
) -> Tuple[list[str], list[str], bool]:
    """
    Returns (available, unavailable, used_fallback).
    If item_assets is None, we can't filter → return requested, [] and True.
    """
    if item_assets is None:
        return list(requested), [], True

    if isinstance(item_assets, Mapping):
        available_set = set(item_assets.keys())
    elif isinstance(item_assets, str):
        available_set = {item_assets}
    else:
        available_set = set(item_assets)

    available = [a for a in requested if a in available_set]
    unavailable = [a for a in requested if a not in available_set]
    return available, unavailable, False


def choose_assets(catalog, collection_id, requested, is_interactive: bool = True):
    avail, unavail, fallback = filter_assets(
        requested, try_get_item_assets(catalog, collection_id)
    )

    if fallback:
        print("Collection has no item_assets; using your original request.")
        return avail

    if not avail:
        print("None of your requested assets are available.")
        return []

    if is_interactive:
        print("Available:", ", ".join(avail))
        if unavail:
            print("Unavailable:", ", ".join(unavail))

        if input("Proceed with available assets? (y/n): ").strip().lower() != "y":
            print("Cancelled.")
            return []

    return avail


def get_provider_catalog(provider: str) -> pystac_client.Client:
    """
    Get the appropriate STAC client based on the provider.
    Returns a pystac_client.Client instance for the specified provider.
    """
    if determine_provider(provider) == "pc":
        return pystac_client.Client.open(
            "https://planetarycomputer.microsoft.com/api/stac/v1",
            modifier=planetary_computer.sign_inplace,
        )


def convert_gdf_geojson(gdf: gpd.GeoDataFrame) -> dict:
    """
    Converts a GeoPandas GeoDataFrame to a GeoJSON Polygon geometry dictionary.

    Parameters:
        gdf (gpd.GeoDataFrame): The GeoDataFrame containing geometries to be unified.

    Returns:
        dict: A GeoJSON dictionary representing the exterior coordinates of the unified polygon.

    Note:
        The function unions all geometries in the GeoDataFrame and extracts the exterior coordinates
        of the resulting polygon.
    """
    roi_shape = gdf.union_all()
    roi_geom = shape(roi_shape)
    roi_geojson_geom = {
        "type": "Polygon",
        "coordinates": [list(roi_geom.exterior.coords)],
    }
    return roi_geojson_geom


def create_landsat_cloud_mask(data, cloud_band: str = "qa_pixel"):
    """
    Create a cloud mask for Landsat imagery based on the specified cloud quality band.

    Args:
        data (xr.DataArray): The Landsat data array containing quality assessment bands.
        cloud_band (str, optional): The name of the cloud quality band to use for masking. Defaults to "qa_pixel".

    Returns:
        xr.DataArray: A boolean mask where True indicates clear pixels (no cloud, cirrus, or dilated cloud detected).
    """
    # Get the bitmask for bad weather conditions
    bad_bits = get_mask_bits()
    # Create an xarray of the cloud band ( unsigned int so we can convert to binary for binary ops)
    qa = data.sel(band=cloud_band).astype("uint16")  # shape: (time, y, x)

    cloudy = (qa & bad_bits) != 0  # same shape, lazy bool
    return cloudy


def add_cloud_mask_to_data(
    data, cloud_mask: xr.DataArray, band_name: str = "cloud_mask"
):
    """
    Add a cloud mask to the Landsat data array.

    Args:
        data (xr.DataArray): The Landsat data array.
        cloud_mask (xr.DataArray): The cloud mask to add.
        band_name (str, optional): The name of the new band for the cloud mask. Defaults to "cloud_mask".

    Returns:
        xr.DataArray: The updated data array with the cloud mask added as a new band.
    """
    # Ensure the cloud mask has the correct dimensions
    cloud_mask = cloud_mask.expand_dims(band=[band_name])
    # Concatenate the cloud mask along the existing band dimension
    return xr.concat([data, cloud_mask], dim="band")


def get_mask_bits(bits: set = (1, 2, 3)):
    """
    Get the bitmask for the specified bad weather conditions.
    Parameters:
        bits (set, optional): Set of bits to consider as "bad weather" (clouds, cirrus, etc.). Defaults to {1, 2, 3}.
            The bits correspond to:
            - 0 : Fill
                - 0: Image data
                - 1: Fill Data
            - 1 : Dilated cloud
            - 2 : Cirrus
            - 3 : Cloud
            - 4 : Cloud shadow
            - 5 : Snow
            - 6 : Clear
            - 7 : Water
            - 8-9 : Cloud confidence
                0: No confidence level set
                1: Low confidence
                2: Medium confidence
                3: High confidence
            - Bits 10-11: Cloud Shadow Confidence
                0: No confidence level set
                1: Low confidence
                2: Reserved
                3: High confidence
            - Bits 12-13: Snow / Ice Confidence
                0: No confidence level set
                1: Low confidence
                2: Reserved
                3: High confidence
            - Bits 14-15: Cirrus Confidence
                0: No confidence level set
                1: Low confidence
                2: Reserved
                3: High confidence


        See https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_TOA

    Returns:
        xr.DataArray: A boolean mask where True indicates clear pixels (no cloud, cirrus, or dilated cloud detected).
    """
    # Define the “bad-weather” bits to reject dilated cloud = bit 1,cirrus = bit 2,cloud = bit 3
    bits = set(bits)
    BAD_BITS = sum(1 << bit for bit in bits)  # Combine bits into a single integer mask
    return BAD_BITS


def meets_clear_threshold(clear_mask, threshold=50):
    """
    Determines if the proportion of clear pixels in a mask meets a specified threshold.

    Args:
        clear_mask (numpy.ndarray or similar): A mask array where clear pixels are indicated as True or 1 (typically as boolean or binary values).
        threshold (int, optional): The minimum percentage of clear pixels required to meet the threshold. Defaults to 50.

    Returns:
        bool: True if the proportion of clear pixels exceeds the threshold, False otherwise.

    Prints:
        A message indicating whether the item is being saved or skipped based on the clear pixel proportion.
    """
    if (clear_mask.values.sum() / clear_mask.values.size) > threshold / 100:
        print("More than 50% clear pixels, saving item.")
        return True
    if clear_mask.sum() / clear_mask.size < 0.5:
        print("Less than 50% clear pixels, skipping item.")
        return False


def get_item_id(item):
    """
    Get the item ID from a STAC item.

    Args:
        item (pystac.Item): The STAC item to extract the ID from.

    Returns:
        str: The ID of the item.
    """
    item_id = item.coords["id"].values  # or item_data.
    item_id_str = str(item_id)
    return item_id_str


def save_bands_to_raster(
    data,
    bands,
    time_index: int,
    output_dir: str,
    driver: str = "GTiff",
    compress: str = "LZW",
):
    """
    Saves specified bands from a time-indexed xarray dataset to a raster file.
    Parameters:
        data (xarray.DataArray): The input data array containing bands and time dimensions.
        bands (list or array-like): List of band names or indices to save.
        time_index (int): The index of the time dimension to select.
        output_dir (str): Directory where the raster file will be saved.
        driver (str, optional): Raster file format driver (default is "GTiff").
        compress (str, optional): Compression type for the raster file (default is "LZW").
    Raises:
        RasterioIOError: If there is an error writing the raster file.
        Exception: For any other errors encountered during saving.
    Notes:
        The output filename is generated from the item ID of the selected time slice.
        If an error occurs during saving, a warning is issued and the error is printed.
    """
    item_data = data.isel(time=time_index)
    item_id_str = get_item_id(item_data)
    filename = f"{item_id_str}.tif"

    data_for_bands = item_data.sel(band=bands)
    raster_path = os.path.join(output_dir, filename)
    try:
        data_for_bands.rio.to_raster(raster_path, driver=driver, compress=compress)
    except RasterioIOError as e:
        warnings.warn(f"Skipping item {item_id_str} due to RasterioIOError: {e}")
        print(f"Skipping item {item_id_str} due to RasterioIOError: {e}")
    except Exception as e:
        print(f"Error saving {bands} bands to {raster_path}: {e}")
    else:
        print(f"Saved {bands} bands to {raster_path}")


def filter_bands(assets, bands):
    """
    Filters the provided list of bands, returning only those that are present in the assets.

    Args:
        assets (dict or list): The available assets, typically a dictionary or list of band names.
        bands (list): The list of band names to filter.

    Returns:
        list: A list containing only the bands that are present in the assets.
    """
    # if the band is not in the assets, then we will not include it
    bands = [band for band in bands if band in assets]
    return bands


def clean_bands(bands: list = None) -> list:
    """
    Normalize any iterable/string of band names:
    - treat a single string as one item
    - strip whitespace, lowercase
    - preserve order, drop duplicates/empties
    """
    if bands is None or bands == "" or (isinstance(bands, Iterable) and not bands):
        return []
    # Allow a lone string instead of an iterable
    if isinstance(bands, str) or not isinstance(bands, Iterable):
        bands = [bands]

    cleaned, seen = [], set()
    for raw in bands:
        name = str(raw).strip().lower()
        if name and name not in seen:
            cleaned.append(name)
            seen.add(name)
    return cleaned


def default_landsat_bands():
    """Return a list of the  default Landsat bands.

    The default bands are:
    - "red"
    - "green"
    - "blue"
    - "swir16" (Shortwave Infrared 1)
    - "nir08" (Near Infrared)
    - "qa_pixel" (Quality Assessment Pixel)

    """
    DEFAULT_LANDSAT_BANDS = ["red", "green", "blue", "swir16", "nir08", "qa_pixel"]
    return DEFAULT_LANDSAT_BANDS.copy()


def resolve_landsat_bands(bands=None):
    """
    Given a list of band names, this function returns a list of valid Landsat band names.

    If no bands are provided, it returns a default set of Landsat bands:
    - "red"
    - "green"
    - "blue"
    - "swir16" (Shortwave Infrared 1)
    - "nir08" (Near Infrared)
    - "qa_pixel" (Quality Assessment Pixel)

    Parameters
    ----------
    bands : list or str, optional
        A list of band names to clean and return. If None, an empty string, or an empty iterable is provided,
        the default set of bands will be returned.


    Returns
    -------
    list[str]
        A list of cleaned, unique band names, ensuring they are in lowercase and stripped of whitespace.
    """
    cleaned = clean_bands(bands)
    return cleaned if cleaned else default_landsat_bands()


def validate_provider_has_collection(provider: str, collection_id: str) -> bool:
    """Check if the specified provider has the given Landsat collection.
    Args:
        provider (str): Provider name (e.g., 'pc').
        collection_id (str): Collection ID to check.
    Returns:
        bool: True if the collection exists for the provider, False otherwise.
    """
    try:
        catalog = get_provider_catalog(provider)
        collection = catalog.get_collection(collection_id)
        return collection is not None
    except Exception as e:
        print(
            f"Error checking collection '{collection_id}' for provider '{provider}': {e}"
        )
        return False


def prompt_user_for_confirmation(num_items: int, missions: list[str]) -> bool:
    """
    Prompts the user to confirm whether to proceed with downloading items.

    Args:
        num_items (int): The number of items matching the criteria.
        missions (list[str]): List of mission names associated with the items.

    Returns:
        bool: True if the user confirms to continue with the download, otherwise exits the program.
    """
    # ask user if they want to continue with the download
    print(
        f"Found {num_items} items matching criteria for missions: {', '.join(missions)}."
    )
    # ask user if they want to continue
    print(
        f"Do you want to continue with the download? (y/n) [default: y]: ",
    )
    user_input = input().strip().lower()
    if user_input not in ["y", "yes", ""]:
        print("Download cancelled by user.")
        exit()


def get_stack_of_items(
    search,
    bands: list,
    missions: list[str],
    interactive_mode: bool = False,
):
    """
        Retrieves a stack of items from the Landsat collection using a STAC search.
    Args:
        search: The STAC search object used to query items.
        bands (list): List of band names to include in the stack.
        missions (list[str]): List of Landsat missions to filter or display.
        interactive_mode (bool, optional): If True, prompts the user for confirmation before proceeding. Defaults to False.
    Returns:
        stackstac.Stack: A stack of items as a Dask array with the specified bands and EPSG:4326 projection.
    Raises:
        ValueError: If no items are found in the search.
        Exception: For errors during band filtering or stacking.
    Notes:
        - The function filters bands to those available in the assets of the first item.
        - Debug print statements are included for items and bands.
        - Signing of collection items may not work for the planetary computer.
    """
    # get the number of items found
    items = get_items_from_search(search)
    # ask user if they want to continue with the download
    if interactive_mode:
        prompt_user_for_confirmation(len(items), missions)

    bands = filter_bands(
        items[0].assets, bands
    )  # remove any bands that are not available in the assets

    return stackstac.stack(
        items,
        assets=bands,
        epsg=4326,
        chunksize=2048,
    )


def create_landsat_platform_filter(missions: Union[int, Sequence[int]]) -> dict:
    """
    Create a platform filter for Landsat missions.
    Args:
        missions (int or Sequence[int]): Single Landsat mission number (5, 7, 8, 9) or a sequence of missions.
    Returns:
        dict: A dictionary representing the platform filter for STAC search.
    Raises:
        TypeError: If missions is not an int or a sequence of ints.
    Notes:
        - The function converts single mission numbers to a list for consistency.
        - If multiple missions are provided, the filter will use "in" to match any of them.
        - If a single mission is provided, it uses "eq" for exact match.
    """
    if isinstance(missions, int):
        missions = [missions]
    elif not isinstance(missions, Sequence):
        raise TypeError("missions must be an int or a sequence of ints")

    platform_names = [f"landsat-{m}" for m in missions]
    return (
        {"in": platform_names} if len(platform_names) > 1 else {"eq": platform_names[0]}
    )


def get_item_assets(catalog, collection_id):
    """Get item assets for a given collection."""
    try:
        collection = catalog.get_child(collection_id)
    except AttributeError as e:
        print(f"Collection {collection_id} not found in catalog.")
        return []
    try:
        item_assets = collection.item_assets
    except AttributeError as e:
        print(
            f"Item assetsproperty not found for collection {collection_id}. Cannot retrieve item assets (aka bands)."
        )
        return []

    return item_assets


def download_landsat_with_config(
    config: LandsatConfig,
    roi_geojson: str,
    start_date: str,
    end_date: str,
    dest_dir: str,
):
    """
    Download Landsat imagery using a configuration object.

    Parameters:
        config: LandsatConfig object with download settings
        roi_geojson: Path to ROI GeoJSON file
        start_date: Start date (YYYY-MM-DD)
        end_date: End date (YYYY-MM-DD)
        dest_dir: Output directory
    """
    gdf = gpd.read_file(roi_geojson)
    roi_geojson_geom = convert_gdf_geojson(gdf)

    search = search_landsat_collection(
        missions=config.missions,
        provider=config.provider,
        roi_geom=roi_geojson_geom,
        start_date=start_date,
        end_date=end_date,
        max_cloud=config.max_cloud,
        collection_id=config.collection_id,
    )

    # once we have a search verify that the requested bands are available
    catalog = get_provider_catalog(config.provider)
    bands = clean_bands(config.bands)
    filtered_bands = choose_assets(catalog, config.collection_id, bands)

    data = get_stack_of_items(
        search,
        bands=filtered_bands,
        missions=[f"Landsat-{m}" for m in config.missions],
        interactive_mode=config.interactive_mode,
    )

    # Clip to the ROI geometry
    data = data.rio.clip(gdf.geometry.apply(mapping), gdf.crs)

    cloudy_array = create_landsat_cloud_mask(data, cloud_band="qa_pixel")
    # add the cloud mask band
    data = add_cloud_mask_to_data(
        data=data,
        cloud_mask=cloudy_array,
        band_name="cloud_mask",
    )

    # here is where we can filter data by the cloud mask
    # @todo create new filtered_data that only contains items with less than MAX_SCENE_CLOUD_COVER%

    out_dir = dest_dir
    os.makedirs(out_dir, exist_ok=True)

    # Loop through each time slice in the data and save the bands to separate files

    for i in tqdm.tqdm(range(len(data.time)), desc="Downloading"):
        # Save each satellite to different directories
        satellite = data.isel(time=i).coords["platform"].item()
        satellite_dir = os.path.join(out_dir, satellite)
        print(f"Saving data for {satellite} to {satellite_dir}")
        os.makedirs(satellite_dir, exist_ok=True)
        print(f"bands: {data.band.values.tolist()}")
        # Save each time slice to a separate file
        save_bands_to_raster(
            data=data,
            bands=data.band.values.tolist(),  # this makes call to get data from api
            time_index=i,
            output_dir=satellite_dir,
            driver="GTiff",
            compress="LZW",
        )


# @todo : Right now this controls the cloud cover for the entire tile (aka not the cloud cover with the clipped ROI)
#       maybe we should make this filter the cloud cover within the ROI
# @todo easy: remove old todo prints

# Make a script to validate the downloaded from each of these sources works properly
# It should validate the script works, the files downloaded and they have the bands they are supposed to


# Example usage with LandsatConfig

if __name__ == "__main__":

    default_bands = default_landsat_bands()

    # Example usage with LandsatConfig
    config_basic = LandsatConfig(
        missions=[8, 9],
        provider="planetary",  # or options are: "planetary"
        collection_id="landsat-c2-l2",  # Surface reflectance
        max_cloud=90,
        bands=default_bands,
        interactive_mode=True,  # Set this to True to be asked to confirm the download before starting
    )

    # Usage examples:
    roi_geojson_path = os.path.join(os.getcwd(), "examples", "rois.geojson")
    start_date = "2025-05-01"
    end_date = "2025-05-30"

    # Download with config
    download_landsat_with_config(
        config=config_basic,
        roi_geojson=roi_geojson_path,
        start_date=start_date,
        end_date=end_date,
        dest_dir="output_landsat8_9_pc",  # Change to your desired output directory
    )
