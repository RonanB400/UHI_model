"""
Interactive map visualization for Urban Heat Island data
"""

import os
import glob
import rasterio
import numpy as np
import folium
from folium import plugins
import branca.colormap as cm
import argparse


def load_raster(filepath):
    """Load a raster file and return the data with metadata"""
    with rasterio.open(filepath) as src:
        data = src.read(1)  # Read the first band
        bounds = src.bounds
        transform = src.transform
        crs = src.crs
        nodata = src.nodata
        
        # Mask no data values
        if nodata is not None:
            data = np.ma.masked_equal(data, nodata)
            # Convert masked array to regular array with NaN values
            data = data.filled(np.nan)
            
    return data, bounds, transform, crs


def create_map(center_lat=45.769467, center_lon=4.834838, zoom=12):
    """Create a base map centered on Lyon"""
    m = folium.Map(location=[center_lat, center_lon], 
                  zoom_start=zoom,
                  tiles='CartoDB positron')
    
    # Add layer control
    folium.LayerControl().add_to(m)
    
    # Add scale - use this method that works in current folium versions
    folium.TileLayer(
        tiles='CartoDB positron',
        attr='CartoDB',
        name='CartoDB Positron',
        overlay=True,
        control=True,
    ).add_to(m)
    
    # We'll skip the scale control since it's problematic
    # and the map will show scale automatically
    
    # Add fullscreen button
    plugins.Fullscreen().add_to(m)
    
    # Add mouse position
    plugins.MousePosition().add_to(m)
    
    return m


def create_simple_visualization(data_dir, output_file):
    """
    Create a simplified interactive map visualization 
    that avoids the ImageOverlay issue
    
    Args:
        data_dir: Directory containing the aligned data files
        output_file: Output HTML file path 
    
    Returns:
        Folium map object
    """
    # Check if the directory exists
    if not os.path.exists(data_dir):
        raise ValueError(f"Directory does not exist: {data_dir}")
    
    # Find aligned files
    lst_files = glob.glob(os.path.join(data_dir, '*LST*aligned.tif'))
    ndvi_files = glob.glob(os.path.join(data_dir, '*NDVI*aligned.tif'))
    impervious_files = glob.glob(os.path.join(data_dir, '*Impervious*aligned.tif'))
    landcover_files = glob.glob(os.path.join(data_dir, '*LandCover*aligned.tif'))
    
    if not any([lst_files, ndvi_files, impervious_files, landcover_files]):
        raise ValueError(f"No aligned files found in {data_dir}")
    
    # Create the map
    m = create_map()
    
    # Add simplified markers for data locations
    if lst_files:
        with rasterio.open(lst_files[0]) as src:
            bounds = src.bounds
            center_lat = (bounds.bottom + bounds.top) / 2
            center_lon = (bounds.left + bounds.right) / 2
            # Add circle marker at the center with popup showing file info
            folium.CircleMarker(
                location=[center_lat, center_lon],
                radius=15,
                popup=f"LST Data: {os.path.basename(lst_files[0])}",
                color='red',
                fill=True,
                fill_color='red',
                fill_opacity=0.6,
            ).add_to(m)
            print(f"Added LST marker from {lst_files[0]}")
            # Recenter map on LST data
            m.location = [center_lat, center_lon]
    
    if ndvi_files:
        with rasterio.open(ndvi_files[0]) as src:
            bounds = src.bounds
            center_lat = (bounds.bottom + bounds.top) / 2
            center_lon = (bounds.left + bounds.right) / 2
            folium.CircleMarker(
                location=[center_lat, center_lon],
                radius=15,
                popup=f"NDVI Data: {os.path.basename(ndvi_files[0])}",
                color='green',
                fill=True,
                fill_color='green',
                fill_opacity=0.6,
            ).add_to(m)
            print(f"Added NDVI marker from {ndvi_files[0]}")
    
    if impervious_files:
        with rasterio.open(impervious_files[0]) as src:
            bounds = src.bounds
            center_lat = (bounds.bottom + bounds.top) / 2
            center_lon = (bounds.left + bounds.right) / 2
            folium.CircleMarker(
                location=[center_lat, center_lon],
                radius=15,
                popup=f"Impervious Data: {os.path.basename(impervious_files[0])}",
                color='gray',
                fill=True,
                fill_color='gray',
                fill_opacity=0.6,
            ).add_to(m)
            print(f"Added impervious surface marker from {impervious_files[0]}")
    
    if landcover_files:
        with rasterio.open(landcover_files[0]) as src:
            bounds = src.bounds
            center_lat = (bounds.bottom + bounds.top) / 2
            center_lon = (bounds.left + bounds.right) / 2
            folium.CircleMarker(
                location=[center_lat, center_lon],
                radius=15,
                popup=f"Land Cover Data: {os.path.basename(landcover_files[0])}",
                color='blue',
                fill=True,
                fill_color='blue',
                fill_opacity=0.6,
            ).add_to(m)
            print(f"Added land cover marker from {landcover_files[0]}")
    
    # Add rectangle to show the data coverage area
    if lst_files:
        with rasterio.open(lst_files[0]) as src:
            bounds = src.bounds
            folium.Rectangle(
                bounds=[[bounds.bottom, bounds.left], [bounds.top, bounds.right]],
                color='red',
                weight=2,
                fill=False,
                popup='LST Data Coverage'
            ).add_to(m)
    
    # Add explanation text
    folium.Marker(
        location=[45.77, 4.85],
        popup="""
        <h3>UHI Map</h3>
        <p>This interactive map shows the location of the data files used 
        for Urban Heat Island analysis.</p>
        <p>Data files are located at:</p>
        <p><code>{}</code></p>
        <p>For full visualization, please use specialized GIS software 
        like QGIS to view the original GeoTIFF files.</p>
        """.format(data_dir),
        icon=folium.Icon(icon='info-sign')
    ).add_to(m)
    
    # Save the map
    m.save(output_file)
    print(f"Map saved to {output_file}")
    print("Open this HTML file in your browser to view the interactive map")
    print("For the full visualization, please use QGIS to view the GeoTIFF files")
    
    return m


def main(output_dir):
    """Create visualization from a specified directory"""
    # Check if the directory exists
    if not os.path.exists(output_dir):
        print(f"Error: Directory does not exist: {output_dir}")
        return None
        
    # Get the aligned directory
    aligned_dir = os.path.join(output_dir, 'Aligned')
    if not os.path.exists(aligned_dir):
        print(f"Error: Aligned directory not found: {aligned_dir}")
        print("Available directories in the output_dir:")
        for item in os.listdir(output_dir):
            if os.path.isdir(os.path.join(output_dir, item)):
                print(f"- {item}")
        return None
    
    # Print the files in the aligned directory for debugging
    print(f"Looking for files in: {aligned_dir}")
    all_files = glob.glob(os.path.join(aligned_dir, '*.*'))
    if all_files:
        print("Files found in Aligned directory:")
        for f in all_files:
            print(f"- {os.path.basename(f)}")
    else:
        print("No files found in the Aligned directory")
        
    # Output path
    output_file = os.path.join(output_dir, 'uhi_map.html')
    
    try:
        # Create the simplified visualization instead
        return create_simple_visualization(aligned_dir, output_file)
    except ValueError as e:
        print(f"Error creating visualization: {e}")
        return None


# Command-line interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Create an interactive map for Urban Heat Island data'
    )
    parser.add_argument(
        'data_dir', 
        type=str,
        help='Directory containing the data (with Aligned subdirectory)'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        help='Output HTML file path (default: data_dir/uhi_map.html)',
        default=None
    )
    
    args = parser.parse_args()
    
    # Check command-line mode vs direct call
    if args.output:
        # Use specified output file
        create_simple_visualization(os.path.join(args.data_dir, 'Aligned'), args.output)
    else:
        # Use default behavior
        main(args.data_dir) 