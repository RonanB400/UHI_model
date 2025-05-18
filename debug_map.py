"""
Debug script for the UHI Map visualization
"""

import os
import sys
from pathlib import Path
from uhi import visualize_map

def check_path(path):
    """Check if a path exists and print its contents if it does"""
    print(f"\nChecking path: {path}")
    if os.path.exists(path):
        print(f"✓ Path exists")
        if os.path.isdir(path):
            files = os.listdir(path)
            if files:
                print(f"Contents ({len(files)} items):")
                for item in files:
                    item_path = os.path.join(path, item)
                    type_marker = "[DIR]" if os.path.isdir(item_path) else "[FILE]"
                    print(f"  {type_marker} {item}")
            else:
                print("Directory is empty")
        else:
            print("Path is a file, not a directory")
    else:
        print("✗ Path does not exist")
    print("-" * 50)

# Script runs from the Project directory
current_dir = os.getcwd()
print(f"Current working directory: {current_dir}")

# Check various potential data paths
possible_paths = [
    # Path in notebook
    os.path.join(current_dir, 'UHI_model', 'data'),
    # Path with potential error
    os.path.join(current_dir, 'data'),
    # Other possible locations
    os.path.join(current_dir, 'UHI_model', 'data', 'Aligned'),
    os.path.join(current_dir, 'data', 'Aligned'),
]

for path in possible_paths:
    check_path(path)

# Try to find any aligned files anywhere in the project
print("\nSearching for aligned files in the project directory...")
for root, dirs, files in os.walk(current_dir):
    aligned_files = [f for f in files if 'aligned.tif' in f.lower()]
    if aligned_files:
        print(f"Found aligned files in: {root}")
        for f in aligned_files:
            print(f"  {f}")

print("\n\nAttempting visualization with different paths...")

# Try the path from the notebook
uhi_model_data_path = os.path.join(current_dir, 'UHI_model', 'data')
if os.path.exists(uhi_model_data_path):
    print(f"\nTrying visualization with: {uhi_model_data_path}")
    visualize_map.main(uhi_model_data_path)

# If that fails, try to run with whatever path contains aligned files
print("\nDone debugging.") 