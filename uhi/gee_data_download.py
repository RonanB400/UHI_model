"""
Google Earth Engine and Planetary Computer Script for downloading and processing
urban heat island detection data
"""

import ee
import geemap
import os
import requests
import numpy as np
import planetary_computer as pc
import pystac_client
import rasterio
from rasterio.warp import reproject, Resampling
from pathlib import Path
import geopandas as gpd
import rioxarray
from shapely.geometry import box
import datetime
import getpass
from oauthlib.oauth2 import BackendApplicationClient
from requests_oauthlib import OAuth2Session
import netCDF4 as nc

# Initialize Earth Engine
try:
    ee.Initialize()
except Exception:
    ee.Authenticate()
    ee.Initialize()

# Define Lyon Area of Interest
# Lyon, France coordinates
u_lon = 4.834838
u_lat = 45.769467
u_poi = ee.Geometry.Point(u_lon, u_lat)

# Define a region of interest with a buffer zone around Le Wagon, Lyon.
# Buffer increased to 5km for a better UHI study area
lyon_geometry = u_poi.buffer(5000)  # Buffer in meters

# Define time period (summer months often show strongest UHI effects)
start_date = '2022-06-01'
end_date = '2022-08-31'

# Common projection for all datasets
target_crs = 'EPSG:4326'
lst_resolution = 100  # LST resolution in meters

output_dir = os.path.join(os.path.dirname(os.getcwd()), 'data')


def download_landsat_data(output_dir):
    """Download Landsat 8/9 data for Lyon at 100m resolution"""
    print("Processing Landsat data...")
    
    # Get Landsat 8/9 Level 2 Collection 2 Tier 1 
    # This dataset includes surface temperature (ST) bands
    landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
        .filterBounds(lyon_geometry) \
        .filterDate(start_date, end_date) \
        .filter(ee.Filter.lt('CLOUD_COVER', 20))  # Filter for low cloud cover

    print(f"Number of Landsat scenes: {landsat.size().getInfo()}")

    # Select a sample image (the least cloudy)
    landsat_image = ee.Image(landsat.sort('CLOUD_COVER').first())
    
    # Get the image date
    image_date = ee.Date(landsat_image.get('system:time_start')).format('yyyy-MM-dd').getInfo()

    # Scale factors for Landsat Collection 2
    def apply_scale_factors(image):
        optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
        thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
        return image.addBands(optical_bands, None, True) \
            .addBands(thermal_bands, None, True)

    landsat_image = apply_scale_factors(landsat_image)

    # Extract Land Surface Temperature (LST)
    lst = landsat_image.select('ST_B10')
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Get download URL at 100m resolution (modified from original 30m)
        url = lst.getDownloadURL({
            'scale': lst_resolution,
            'crs': target_crs,
            'region': lyon_geometry,
            'format': 'GEO_TIFF'
        })
        
        # Download the file with date and source in filename
        response = requests.get(url)
        if response.status_code == 200:
            lst_file = os.path.join(output_dir, f'Lyon_LST_Landsat_{image_date}.tif')
            with open(lst_file, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded LST to {lst_file}")
            return lst_file
        else:
            print(f"Error downloading LST: HTTP {response.status_code}")
            print("The file might be too large for direct download.")
    except Exception as e:
        print(f"Error: {e}")
        print("File might be too large. Try a smaller region.")
    
    # Function to extract LST from Landsat collection
    def process_landsat_collection(collection):
        # Apply scale factors for LST
        def apply_thermal_scale(image):
            date = ee.Date(image.get('system:time_start'))
            thermal = image.select('ST_B10').multiply(0.00341802).add(149.0)
            return thermal.rename('LST').set('date', date.format('YYYY-MM-dd'))
        
        return collection.map(apply_thermal_scale)

    # Get LST time series
    lst_collection = process_landsat_collection(landsat)

    # Export LST time series
    dates = lst_collection.aggregate_array('date').getInfo()
    
    # Create time series directory
    ts_dir = os.path.join(output_dir, 'TimeSeries')
    os.makedirs(ts_dir, exist_ok=True)
    
    # Note: We'll only download the first image as example
    # (downloading all dates might exceed API limits)
    if dates:
        try:
            image = ee.Image(lst_collection.filter(
                ee.Filter.eq('date', dates[0])).first())
            
            url = image.getDownloadURL({
                'scale': lst_resolution,
                'crs': target_crs,
                'region': lyon_geometry,
                'format': 'GEO_TIFF'
            })
            
            response = requests.get(url)
            if response.status_code == 200:
                ts_file = os.path.join(ts_dir, f'Lyon_LST_Landsat_{dates[0]}.tif')
                with open(ts_file, 'wb') as f:
                    f.write(response.content)
                print(f"Downloaded LST for {dates[0]} to {ts_dir}")
                print("Note: Only downloaded 1 date as example.")
            else:
                print("Error downloading time series: HTTP "
                      f"{response.status_code}")
        except Exception as e:
            print(f"Error downloading time series: {e}")
    return None


def download_sentinel_data(output_dir):
    """Download Sentinel-2 data for Lyon at 10m resolution"""
    print("Processing Sentinel-2 data...")
    
    # Get Sentinel-2 Level 2A data
    sentinel = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
        .filterBounds(lyon_geometry) \
        .filterDate(start_date, end_date) \
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))

    print(f"Number of Sentinel-2 scenes: {sentinel.size().getInfo()}")

    # Select a sample image
    sentinel_image = ee.Image(sentinel.sort('CLOUDY_PIXEL_PERCENTAGE').first())
    
    # Get the image date
    image_date = ee.Date(sentinel_image.get('system:time_start')).format('yyyy-MM-dd').getInfo()

    # Calculate NDVI (Normalized Difference Vegetation Index)
    ndvi = sentinel_image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Download NDVI at 10m resolution
        url = ndvi.getDownloadURL({
            'scale': 10,  # 10m resolution
            'crs': target_crs,
            'region': lyon_geometry,
            'format': 'GEO_TIFF'
        })
        
        response = requests.get(url)
        if response.status_code == 200:
            ndvi_file = os.path.join(output_dir, f'Lyon_NDVI_Sentinel_{image_date}.tif')
            with open(ndvi_file, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded NDVI to {ndvi_file}")
            return ndvi_file
        else:
            print(f"Error downloading NDVI: HTTP {response.status_code}")
    except Exception as e:
        print(f"Error downloading NDVI: {e}")
    
    return None


def download_impervious_surface(output_dir):
    """Download impervious surface data from Global Human Settlement Layer"""
    print("Processing impervious surface data...")
    
    # Use the JRC Global Human Settlement Layer data
    imperviousness = ee.ImageCollection("JRC/GHSL/P2023A/GHS_BUILT_C") \
        .filterDate('2018-01-01', '2018-12-31') \
        .first() \
        .select('built_characteristics')
    
    # Get the image date
    image_date = ee.Date(imperviousness.get('system:time_start')).format('yyyy-MM-dd').getInfo()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Download impervious surface at 10m resolution
        url = imperviousness.getDownloadURL({
            'scale': 10,
            'crs': target_crs,
            'region': lyon_geometry,
            'format': 'GEO_TIFF'
        })
        
        response = requests.get(url)
        if response.status_code == 200:
            impervious_file = os.path.join(output_dir, f'Lyon_Impervious_GHSL_{image_date}.tif')
            with open(impervious_file, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded impervious surface to {impervious_file}")
            return impervious_file
        else:
            print("Error downloading impervious data: HTTP "
                  f"{response.status_code}")
    except Exception as e:
        print(f"Error downloading impervious data: {e}")
    
    return None


def download_land_cover(output_dir):
    """Download ESA WorldCover land cover data using Google Earth Engine"""
    print("Processing ESA WorldCover land cover data...")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Get ESA WorldCover from Earth Engine
        worldcover = ee.ImageCollection("ESA/WorldCover/v200") \
            .filterBounds(lyon_geometry) \
            .first() \
            .select('Map')
        
        # Get the image date
        image_date = "2021"  # ESA WorldCover 2021 by default
        if worldcover.get('system:time_start'):
            image_date = ee.Date(worldcover.get('system:time_start')).format('yyyy-MM-dd').getInfo()
        
        # Download land cover at 10m resolution
        url = worldcover.getDownloadURL({
            'scale': 10,  # 10m resolution
            'crs': target_crs,
            'region': lyon_geometry,
            'format': 'GEO_TIFF'
        })
        
        response = requests.get(url)
        if response.status_code == 200:
            lc_file = os.path.join(output_dir, f'Lyon_LandCover_ESA_{image_date}.tif')
            with open(lc_file, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded land cover data to {lc_file}")
            return lc_file
        else:
            print(f"Error downloading land cover: HTTP {response.status_code}")
    except Exception as e:
        print(f"Error downloading land cover data: {e}")
    
    return None


def align_and_resample(output_dir, lst_file, ndvi_file, impervious_file, 
                       land_cover_file):
    """Align all layers to a common spatial extent and CRS, preserving native resolutions"""
    print("Aligning data layers...")
    
    aligned_dir = os.path.join(output_dir, 'Aligned')
    os.makedirs(aligned_dir, exist_ok=True)
    
    # Create full paths to input files
    lst_path = os.path.join(output_dir, lst_file) if lst_file else None
    ndvi_path = os.path.join(output_dir, ndvi_file) if ndvi_file else None
    impervious_path = os.path.join(output_dir, impervious_file) if impervious_file else None
    land_cover_path = os.path.join(output_dir, land_cover_file) if land_cover_file else None
    
    # First, determine the common bounding box from LST (our reference)
    if lst_path and os.path.exists(lst_path):
        with rasterio.open(lst_path) as src:
            # Get the bounding box in the CRS coordinates
            dst_crs = src.crs
            bounds = src.bounds
            
            # Prepare LST aligned filename
            lst_aligned_file = os.path.splitext(lst_file)[0] + '_aligned' + os.path.splitext(lst_file)[1]
            
            # Copy LST to aligned directory
            dst_lst = os.path.join(aligned_dir, lst_aligned_file)
            with rasterio.open(dst_lst, 'w', **src.meta.copy()) as dst:
                dst.write(src.read())
            
            print(f"Aligned LST saved to {dst_lst}")
            
            # Align NDVI
            if ndvi_path and os.path.exists(ndvi_path):
                ndvi_aligned_file = os.path.splitext(ndvi_file)[0] + '_aligned' + os.path.splitext(ndvi_file)[1]
                align_to_bbox(ndvi_path, 
                             os.path.join(aligned_dir, ndvi_aligned_file),
                             bounds, dst_crs)
            
            # Align impervious surface
            if impervious_path and os.path.exists(impervious_path):
                impervious_aligned_file = os.path.splitext(impervious_file)[0] + '_aligned' + os.path.splitext(impervious_file)[1]
                align_to_bbox(impervious_path,
                             os.path.join(aligned_dir, impervious_aligned_file),
                             bounds, dst_crs)
            
            # Align land cover
            if land_cover_path and os.path.exists(land_cover_path):
                landcover_aligned_file = os.path.splitext(land_cover_file)[0] + '_aligned' + os.path.splitext(land_cover_file)[1]
                align_to_bbox(land_cover_path,
                             os.path.join(aligned_dir, landcover_aligned_file),
                             bounds, dst_crs)
            
            print(f"All layers aligned to the same extent and saved to {aligned_dir}")
            return aligned_dir
    else:
        print("LST file not found. Cannot align layers.")
    
    return None


def align_to_bbox(src_path, dst_path, bounds, dst_crs):
    """Align a raster to match a bounding box and CRS without changing its resolution"""
    try:
        with rasterio.open(src_path) as src:
            # Create a copy of the source metadata
            kwargs = src.meta.copy()
            
            # Update the CRS to match the destination
            kwargs.update({
                'crs': dst_crs,
            })
            
            # Calculate the new geotransform to match the bounds while keeping pixel size
            # This only adjusts the origin, not the resolution
            src_res_x, src_res_y = src.res
            
            # Create a dataset with the same resolution but new bounds
            with rasterio.open(dst_path, 'w', **kwargs) as dst:
                # Read the source data
                data = src.read()
                
                # Reproject to new CRS but maintain resolution
                dst_transform = rasterio.transform.from_bounds(
                    bounds.left, bounds.bottom, bounds.right, bounds.top, 
                    round((bounds.right - bounds.left) / src_res_x),
                    round((bounds.top - bounds.bottom) / src_res_y)
                )
                
                reproject(
                    source=data,
                    destination=np.zeros(
                        (src.count, 
                         round((bounds.top - bounds.bottom) / src_res_y),
                         round((bounds.right - bounds.left) / src_res_x))),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest
                )[0]
            
            print(f"Aligned {src_path} to {dst_path} (preserved resolution: {src_res_x}m)")
            return dst_path
    except Exception as e:
        print(f"Error aligning {src_path}: {e}")
    
    return None


def visualize_data():
    """Visualize data on a map"""
    # Create a map
    map_lyon = geemap.Map()
    map_lyon.centerObject(lyon_geometry, 11)
    map_lyon.addLayer(lyon_geometry, {}, 'Lyon Area')
    
    # Get Landsat data
    landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
        .filterBounds(lyon_geometry) \
        .filterDate(start_date, end_date) \
        .filter(ee.Filter.lt('CLOUD_COVER', 20))
    
    landsat_image = ee.Image(landsat.sort('CLOUD_COVER').first())
    
    # Apply scale factors
    def apply_scale_factors(image):
        optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
        thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
        return image.addBands(optical_bands, None, True) \
            .addBands(thermal_bands, None, True)

    landsat_image = apply_scale_factors(landsat_image)
    
    # Display the Land Surface Temperature (LST)
    lst = landsat_image.select('ST_B10')
    lst_vis = {
        'min': 20, 
        'max': 40, 
        'palette': ['blue', 'yellow', 'orange', 'red']
    }
    map_lyon.addLayer(lst.clip(lyon_geometry), lst_vis, 
                     'Land Surface Temperature')

    # True color composite
    rgb_vis = {'bands': ['SR_B4', 'SR_B3', 'SR_B2'], 'min': 0, 'max': 0.3}
    map_lyon.addLayer(landsat_image.clip(lyon_geometry), rgb_vis, 'True Color')
    
    # Get Sentinel-2 data
    sentinel = ee.ImageCollection('COPERNICUS/S2_SR') \
        .filterBounds(lyon_geometry) \
        .filterDate(start_date, end_date) \
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    
    sentinel_image = ee.Image(sentinel.sort('CLOUDY_PIXEL_PERCENTAGE').first())
    
    # Calculate NDVI
    ndvi = sentinel_image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    ndvi_vis = {
        'min': -0.2, 
        'max': 0.8, 
        'palette': ['brown', 'yellow', 'green', 'darkgreen']
    }
    map_lyon.addLayer(ndvi.clip(lyon_geometry), ndvi_vis, 'NDVI')

    # Get impervious surface data
    imperviousness = ee.ImageCollection("JRC/GHSL/P2023A/GHS_BUILT_C") \
        .filterDate('2018-01-01', '2018-12-31') \
        .first() \
        .select('built_characteristics')
    
    imp_vis = {
        'min': 0,
        'max': 25,
        'palette': ['black', 'grey', 'white']
    }
    map_lyon.addLayer(imperviousness.clip(lyon_geometry), imp_vis, 
                     'Impervious Surface')
    
    # Sentinel-2 true color
    s2_rgb_vis = {'bands': ['B4', 'B3', 'B2'], 'min': 0, 'max': 3000}
    map_lyon.addLayer(
        sentinel_image.clip(lyon_geometry), 
        s2_rgb_vis, 
        'Sentinel-2 True Color'
    )
    
    return map_lyon


def download_ecostress_data(output_dir):
    """Download NASA ECOSTRESS LST data for Lyon using the EARTHDATA API"""
    print("Processing ECOSTRESS LST data...")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # NASA EARTHDATA authentication
    def get_earthdata_credentials():
        """Get NASA EARTHDATA login credentials"""
        print("Please enter your NASA EARTHDATA Login credentials")
        username = input("Username: ")
        password = getpass.getpass("Password: ")
        return username, password
    
    def get_earthdata_token(username, password):
        """Get authentication token for EARTHDATA API"""
        # NASA EARTHDATA Authentication
        urs_url = "https://urs.earthdata.nasa.gov"
        app_client_id = "earthdata-search-client"
        app_client_secret = "none"
        
        # Create OAuth2 session
        client = BackendApplicationClient(client_id=app_client_id)
        oauth = OAuth2Session(client=client)
        
        # Get token
        token = oauth.fetch_token(
            token_url=f"{urs_url}/oauth/token",
            username=username,
            password=password,
            client_id=app_client_id,
            client_secret=app_client_secret,
        )
        return token
    
    def search_ecostress_lst_data(token, bbox, start_date, end_date):
        """Search for ECOSTRESS LST data using CMR API"""
        # CMR API endpoint
        cmr_url = "https://cmr.earthdata.nasa.gov/search"
        
        # Convert bbox to CMR format
        # CMR expects: [west, south, east, north]
        bbox_str = f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}"
        
        # Format dates for CMR
        start_dt = datetime.datetime.strptime(start_date, "%Y-%m-%d")
        end_dt = datetime.datetime.strptime(end_date, "%Y-%m-%d")
        start_date_str = start_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
        end_date_str = end_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
        
        # ECOSTRESS LST collection short name
        collection = "ECO2LSTE.001"  # ECOSTRESS L2 Land Surface Temperature and Emissivity
        
        # Set up parameters for granule search
        params = {
            "collection_concept_id": collection,
            "temporal": f"{start_date_str},{end_date_str}",
            "bounding_box": bbox_str,
            "page_size": 100,
            "sort_key": "-start_date"  # Sort by start date descending
        }
        
        # Construct URL
        headers = {"Authorization": f"Bearer {token['access_token']}"}
        granules_url = f"{cmr_url}/granules.json"
        
        # Make request to CMR API
        response = requests.get(granules_url, params=params, headers=headers)
        
        if response.status_code == 200:
            results = response.json()['feed']['entry']
            print(f"Found {len(results)} ECOSTRESS LST granules")
            return results
        else:
            print(f"Error searching for data: {response.status_code}")
            print(response.text)
            return []
    
    def download_ecostress_file(url, token, output_path):
        """Download ECOSTRESS file with authentication"""
        headers = {"Authorization": f"Bearer {token['access_token']}"}
        response = requests.get(url, headers=headers, stream=True)
        
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=1024*1024):
                    f.write(chunk)
            print(f"Downloaded file to {output_path}")
            return output_path
        else:
            print(f"Error downloading file: {response.status_code}")
            return None
    
    def process_ecostress_netcdf(file_path, output_dir):
        """Process ECOSTRESS NetCDF/HDF file and convert to GeoTIFF"""
        try:
            # Open the NetCDF/HDF file
            dataset = nc.Dataset(file_path, 'r')
            
            # Extract LST data
            lst_data = dataset.variables['LST'][:]
            
            # Get geolocation information
            lat = dataset.variables['lat'][:]
            lon = dataset.variables['lon'][:]
            
            # Get metadata
            acquisition_date = dataset.time_coverage_start
            
            # Close the dataset
            dataset.close()
            
            # Convert to simple datetime string for filename
            date_str = acquisition_date.split('T')[0]
            
            # Create output filename
            output_file = os.path.join(
                output_dir, 
                f'Lyon_LST_ECOSTRESS_{date_str}.tif'
            )
            
            # Create a GeoTIFF file with proper georeference
            # Note: This is a simplified version, actual implementation would need
            # proper transformation of coordinates
            transform = rasterio.transform.from_bounds(
                lon.min(), lat.min(), lon.max(), lat.max(), 
                lst_data.shape[1], lst_data.shape[0]
            )
            
            with rasterio.open(
                output_file,
                'w',
                driver='GTiff',
                height=lst_data.shape[0],
                width=lst_data.shape[1],
                count=1,
                dtype=lst_data.dtype,
                crs=target_crs,
                transform=transform,
                nodata=0
            ) as dst:
                dst.write(lst_data, 1)
                
            print(f"Converted ECOSTRESS data to GeoTIFF: {output_file}")
            return output_file
        
        except Exception as e:
            print(f"Error processing ECOSTRESS file: {e}")
            return None
    
    try:
        # Get bounding box for Lyon directly from coordinates
        # Instead of using EE geometry that causes issues
        bbox = [
            u_lon - 0.05,  # west - approx 5km buffer
            u_lat - 0.05,  # south - approx 5km buffer
            u_lon + 0.05,  # east - approx 5km buffer
            u_lat + 0.05   # north - approx 5km buffer
        ]
        
        # Get credentials
        username, password = get_earthdata_credentials()
        
        # Get auth token
        token = get_earthdata_token(username, password)
        
        # Search for ECOSTRESS data
        granules = search_ecostress_lst_data(token, bbox, start_date, end_date)
        
        if not granules:
            print("No ECOSTRESS LST data found for the specified area and time period.")
            return None
        
        # Process only the first granule for demonstration
        granule = granules[0]
        
        # Get download URL from granule links
        download_url = None
        for link in granule.get('links', []):
            if link.get('rel') == 'download':
                download_url = link.get('href')
                break
        
        if not download_url:
            print("Could not find download URL for the selected granule.")
            return None
        
        # Create filename from granule ID
        granule_id = granule.get('id', 'unknown')
        download_path = os.path.join(output_dir, f"{granule_id}.h5")
        
        # Download the file
        downloaded_file = download_ecostress_file(download_url, token, download_path)
        
        if downloaded_file:
            # Process the file and convert to GeoTIFF
            ecostress_tif = process_ecostress_netcdf(downloaded_file, output_dir)
            return os.path.basename(ecostress_tif) if ecostress_tif else None
        
    except Exception as e:
        print(f"Error downloading ECOSTRESS data: {e}")
    
    return None


if __name__ == "__main__":
    # Output directory
    output_directory = './data'
    use_google_drive = False  # Set to True to export to Google Drive
    
    print("Downloading data locally to:", output_directory)
    
    # Download all datasets
    lst_file = download_landsat_data(output_dir=output_directory)
    ndvi_file = download_sentinel_data(output_dir=output_directory)
    impervious_file = download_impervious_surface(output_dir=output_directory)
    land_cover_file = download_land_cover(output_dir=output_directory)
    
    # Download ECOSTRESS LST data
    ecostress_file = download_ecostress_data(output_dir=output_directory)
    
    # Extract just the filenames from the paths
    lst_filename = os.path.basename(lst_file) if lst_file else None
    ndvi_filename = os.path.basename(ndvi_file) if ndvi_file else None
    impervious_filename = os.path.basename(impervious_file) if impervious_file else None
    land_cover_filename = os.path.basename(land_cover_file) if land_cover_file else None
    
    # Align and resample all layers
    align_and_resample(output_directory, lst_filename, ndvi_filename, 
                       impervious_filename, land_cover_filename)
    
    # Visualize data
    map_lyon = visualize_data()
    try:
        map_lyon.to_streamlit(height=600)
    except Exception:
        # If not in streamlit environment, save as HTML
        html_file = os.path.join(output_directory, 'map.html')
        map_lyon.save(html_file)
        print(f"Map saved to {html_file}")
    
    print("Data processing complete.")
    if use_google_drive:
        print("Check your Google Drive for exported files.")
    else:
        print(f"Check {output_directory} for downloaded files.") 