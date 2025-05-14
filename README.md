# Detection and Modeling of Urban Heat Islands in Lyon Using Open Data and Satellite Imagery

## **Project Title:**

**Detection and Modeling of Urban Heat Islands in Lyon Using Open Data and Satellite Imagery**

## **Objective:**

To identify and analyze urban heat islands in Lyon, France, by combining remote sensing data, open meteorological datasets, and urban morphology. The goal is to highlight heat-prone areas and suggest mitigation strategies.

## **Scope:**

A **proof of concept** focusing on:

* A **specific neighborhood or pilot area** in Lyon
* **Data-driven mapping** of surface temperatures and vegetation
* Exploratory analysis of the relationship between temperature, urban structure, and vegetation

---

## **Data Sources:**

### **Satellite Data:**

* **Landsat:** 30 m resolution, surface temperature
* **Sentinel-2:** 10 m resolution, vegetation indices (NDVI), images, morphology


### **Open Data:**

* **Data Grand Lyon:** Hourly meteorological data from local stations
* Urban morphology: building footprints, heights, street width, etc.
* Land use and surface types (asphalt, vegetation, water, etc.)

---

## **Methodology (Simplified Pipeline):**

### 1. **Data Collection:**

* Download satellite imagery (Landsat, Sentinel-2)
* Retrieve meteorological and urban morphology data from open sources

### 2. **Data Preprocessing:**

* Extract land surface temperature (LST)
* Calculate NDVI for vegetation density
* Harmonize spatial resolution and align datasets

### 3. **Exploratory Data Analysis (EDA):**

* Visualize temperature and vegetation distribution
* Identify hotspots (UHI zones)
* Study correlations with built environment (height, density, albedo)

### 4. **Output:**

* A model able to give a more precise measure of the surface temperature, and able to predict UHI
A map of UHIs in the study area
* Visual dashboard or report with insights and suggestions

---

## **Deliverables:**

* Jupyter notebook or Python script with full pipeline
* Visuals: maps, graphs (temperature, NDVI, correlation)
* Brief report or slide deck with findings, methodology, and next steps
