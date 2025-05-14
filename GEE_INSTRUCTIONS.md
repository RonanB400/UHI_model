# Using Google Earth Engine to Download Satellite Data for Lyon UHI Project

This guide explains how to use Google Earth Engine (GEE) to download and process Landsat 8/9 and Sentinel-2 data for the Urban Heat Island (UHI) study in Lyon, France.

## Prerequisites

1. **Google Earth Engine Account**
   - Sign up at [https://earthengine.google.com/](https://earthengine.google.com/)
   - Request access if you don't already have it (approval usually takes 1-2 days)

2. **Python Environment Setup**
   - Install required packages:
     ```bash
     pip install earthengine-api geemap requests
     ```

3. **Google Earth Engine Authentication**
   - The first time you run the script, you'll need to authenticate:
     ```python
     import ee
     ee.Authenticate()  # This will open a browser window
     ee.Initialize()
     ```

## Using the Provided Script

We've created `gee_data_download.py` that handles downloading and processing both Landsat and Sentinel-2 data.

### Running the Script

```bash
python gee_data_download.py
```

This will:
1. Download Landsat 8/9 data (including Land Surface Temperature)
2. Download Sentinel-2 data (including NDVI calculations)
3. Export the data to your Google Drive (default) or locally
4. Display a map visualization (if run in a notebook environment)

### Script Features

- **Time Period**: Default is summer 2022 (June-August) to maximize UHI detection
- **Area of Interest**: Lyon, France (approximately 4.7°E, 45.65°N to 5.0°E, 45.85°N)
- **Cloud Filtering**: Only scenes with <20% cloud cover
- **Data Exports**:
  - Land Surface Temperature (LST) from Landsat
  - NDVI (vegetation index) from Sentinel-2
  - True color composites
  - Time series data (all available dates)

## Downloading Data Locally vs. to Google Drive

The script now supports two download methods:

### Google Drive Download (Default)
- Handles large files (no size limits)
- Processing happens on Google's servers
- Results appear in your Google Drive when complete
- Better for full-scale analysis

### Local Download
- Limited to files under 32MB by the Earth Engine API
- May fail for large areas or high-resolution data
- Faster for small areas
- No need to download from Google Drive afterward

To switch to local downloads, edit this section in the script:

```python
# Change to False for local downloads
use_google_drive = False
output_directory = './data'  # Change to your preferred directory
```

**Important Note on File Sizes:**
- Landsat data (30m resolution) for all of Lyon may exceed the 32MB direct download limit
- Sentinel-2 data (10m resolution) will almost certainly exceed the limit
- Options for local downloads:
  1. Reduce the area size (e.g., select a smaller neighborhood)
  2. Reduce the output resolution
  3. Use Google Drive for the full dataset

## Customizing the Script

You can modify these parameters in the script:

```python
# Change the area of interest
lyon_geometry = ee.Geometry.Rectangle([4.7, 45.65, 5.0, 45.85])

# Change the time period
start_date = '2022-06-01'
end_date = '2022-08-31'
```

For a specific neighborhood, you can define a more precise geometry:

```python
# Example for a specific neighborhood
neighborhood = ee.Geometry.Polygon([
    [4.82, 45.75], 
    [4.82, 45.77],
    [4.84, 45.77],
    [4.84, 45.75]
])
```

## Accessing Downloaded Data

### For Google Drive Downloads:
1. Go to your Google Drive
2. Look for folders named `UHI_Lyon` and `UHI_Lyon_TimeSeries`
3. GeoTIFF files will be available for each dataset

### For Local Downloads:
1. Check the specified output directory (default is `./data`)
2. Files will be saved as GeoTIFF with appropriate names
3. A map.html file will also be generated for visualization

## Processing Tips

- **Land Surface Temperature (LST)**:
  - Landsat LST is in Kelvin, subtract 273.15 to convert to Celsius
  - For UHI detection, look for temperature differences between urban and rural areas

- **NDVI Analysis**:
  - NDVI ranges from -1 to 1
  - Higher values (>0.5) indicate dense vegetation
  - Lower values (<0.2) in urban areas typically correlate with higher temperatures

## Troubleshooting

- **Authentication Issues**: If you see authentication errors, run `ee.Authenticate()` separately
- **Export Timeouts**: Large exports might time out; try reducing the area or resolution
- **Memory Errors**: If you encounter memory errors, try processing smaller areas or fewer dates
- **Local Download Failures**: If you see "File too large" errors when downloading locally, switch to Google Drive export

## Additional Resources

- [Google Earth Engine Documentation](https://developers.google.com/earth-engine/guides)
- [geemap Documentation](https://geemap.org/) 