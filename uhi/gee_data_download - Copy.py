"""
Google Earth Engine Script for downloading Landsat and Sentinel-2 data for Lyon
"""

import ee
import geemap
import os
import requests
from pathlib import Path

# Initialize Earth Engine
try:
    ee.Initialize()
except Exception:
    ee.Authenticate()
    ee.Initialize()

# Define Lyon Area of Interest using the simplest possible format
# Lyon, France coordinates (small test area)
#bbox = [4.82, 45.77, 4.84, 45.78]  # [west, south, east, north]

# Define the location of interest (Lyon, France).
u_lon = 4.834838
u_lat = 45.769467
u_poi = ee.Geometry.Point(u_lon, u_lat)

# Define a region of interest with a buffer zone of 1 km around Le Wagon, Lyon.
lyon_geometry = u_poi.buffer(1000)  # Buffer in meters

# Define time period (summer months often show strongest UHI effects)
start_date = '2022-06-01'
end_date = '2022-08-31'


def download_landsat_data(output_dir):
    """Download Landsat 8/9 data for Lyon"""
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
        # Get download URL (limited to 32MB)
        url = lst.getDownloadURL({
            'scale': 30,
            'crs': 'EPSG:4326',
            'region': lyon_geometry,
            'format': 'GEO_TIFF'
        })
        
        # Download the file
        response = requests.get(url)
        if response.status_code == 200:
            with open(os.path.join(output_dir, 'Lyon_LST.tif'), 'wb') as f:
                f.write(response.content)
            print(f"Downloaded LST to {output_dir}/Lyon_LST.tif")
        else:
            print(f"Error downloading LST: HTTP {response.status_code}")
            print("The file might be too large for direct download.")
            print("Consider using Google Drive export instead.")
    except Exception as e:
        print(f"Error: {e}")
        print("File might be too large. Try a smaller region or use Google Drive.")
    
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
                'scale': 30,
                'crs': 'EPSG:4326',
                'region': lyon_geometry,
                'format': 'GEO_TIFF'
            })
            
            response = requests.get(url)
            if response.status_code == 200:
                with open(os.path.join(ts_dir, f'Lyon_LST_{dates[0]}.tif'), 'wb') as f:
                    f.write(response.content)
                print(f"Downloaded LST for {dates[0]} to {ts_dir}")
                print("Note: Only downloaded 1 date as example.")
                print("Too many dates might exceed download limits.")
            else:
                print(f"Error downloading time series: HTTP {response.status_code}")
        except Exception as e:
            print(f"Error downloading time series: {e}")
            print("Consider using Google Drive export for time series data.")


def download_sentinel_data(output_dir):
    """Download Sentinel-2 data for Lyon"""
    print("Processing Sentinel-2 data...")
    
    # Get Sentinel-2 Level 2A data
    sentinel = ee.ImageCollection('COPERNICUS/S2_SR') \
        .filterBounds(lyon_geometry) \
        .filterDate(start_date, end_date) \
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))

    print(f"Number of Sentinel-2 scenes: {sentinel.size().getInfo()}")

    # Select a sample image
    sentinel_image = ee.Image(sentinel.sort('CLOUDY_PIXEL_PERCENTAGE').first())

    # Calculate NDVI (Normalized Difference Vegetation Index)
    ndvi = sentinel_image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Download NDVI
        url = ndvi.getDownloadURL({
            'scale': 10,
            'crs': 'EPSG:4326',
            'region': lyon_geometry,
            'format': 'GEO_TIFF'
        })
        
        response = requests.get(url)
        if response.status_code == 200:
            with open(os.path.join(output_dir, 'Lyon_NDVI.tif'), 'wb') as f:
                f.write(response.content)
            print(f"Downloaded NDVI to {output_dir}/Lyon_NDVI.tif")
        else:
            print(f"Error downloading NDVI: HTTP {response.status_code}")
            print("File might be too large for direct download.")
    except Exception as e:
        print(f"Error downloading NDVI: {e}")
        print("Consider using Google Drive export instead.")
        
    try:
        # Download RGB (true color)
        rgb_bands = sentinel_image.select(['B4', 'B3', 'B2'])
        url = rgb_bands.getDownloadURL({
            'scale': 10,
            'crs': 'EPSG:4326',
            'region': lyon_geometry,
            'format': 'GEO_TIFF'
        })
        
        response = requests.get(url)
        if response.status_code == 200:
            with open(os.path.join(output_dir, 'Lyon_S2_RGB.tif'), 'wb') as f:
                f.write(response.content)
            print(f"Downloaded RGB to {output_dir}/Lyon_S2_RGB.tif")
        else:
            print(f"Error downloading RGB: HTTP {response.status_code}")
            print("File might be too large for direct download.")
    except Exception as e:
        print(f"Error downloading RGB: {e}")
        print("Consider using Google Drive export instead.")


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
    map_lyon.addLayer(lst.clip(lyon_geometry), lst_vis, 'Land Surface Temperature')

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

    # Sentinel-2 true color
    s2_rgb_vis = {'bands': ['B4', 'B3', 'B2'], 'min': 0, 'max': 3000}
    map_lyon.addLayer(
        sentinel_image.clip(lyon_geometry), 
        s2_rgb_vis, 
        'Sentinel-2 True Color'
    )
    
    return map_lyon


if __name__ == "__main__":
    # Set to False to download directly to local machine instead of Google Drive
    output_directory = './data'
    
    print("Downloading data locally to:", output_directory)
    print("Note: Files may be too large for direct download.")
    print("If downloads fail, try using Google Drive (set use_google_drive=True)")
    
    download_landsat_data(output_dir=output_directory)
    download_sentinel_data(output_dir=output_directory)
    
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