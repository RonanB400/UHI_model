# 🌡️ Urban Heat Island Downscaling & Mitigation Project

## 🎯 Project Goal

This project aims to produce high-resolution (10 m) Urban Heat Island (UHI) maps for cities like Lyon by combining ECOSTRESS thermal satellite data with high-resolution predictors (NDVI, impervious surfaces, land cover). The final objective is to recommend optimal zones for tree planting to mitigate extreme urban heat, especially during summer heatwaves.

---

## 🧭 Step 1 — Downscale ECOSTRESS LST from 1 km (MODIS) to 70 m

### 🔹 1.1 Data Collection

- **ECOSTRESS LST** (70 m) via Google Earth Engine  
- **MODIS LST** (1 km) as baseline for comparison  
- **Sentinel-2 (10 m)**: compute NDVI, NDBI, NDWI  
- **GHSL Impervious Surface** (10–30 m)  
- **ESA WorldCover 10 m** land cover map  
- Optional: Microsoft Tree Canopy, VIIRS Nightlights, slope/elevation

### 🔹 1.2 Preprocessing

- Clip to city (e.g., Lyon), align all datasets to common grid/CRS  
- Resample predictors and align with ECOSTRESS  
- Filter ECOSTRESS by quality flags (QC)  
- Build training samples:  
  `X = [NDVI, canopy, imperviousness, land use, etc.]`  
  `y = ECOSTRESS LST (°C)`

### 🔹 1.3 Modeling & Evaluation

- Train models: `RandomForest`, `XGBoost`, `LightGBM`  
- Evaluate with R², RMSE, spatial residuals  
- Compare MODIS vs ECOSTRESS vs model predictions  
- Output: trained model (.pkl), performance plots, residual maps

---

## 🧭 Step 2 — Downscale ECOSTRESS from 70 m to 10 m

### 🔹 2.1 Prepare 10 m Inference Grid

- Use high-res predictors (Sentinel-2 NDVI, impervious, canopy, etc.)  
- Build tabular or raster dataset for inference

### 🔹 2.2 Predict Fine-Scale LST

- Load model trained in Step 1  
- Predict LST at 10 m: `LST_10m = model.predict(X_10m)`  
- Save predicted map as `UHI_LST_10m.tif`

### 🔹 2.3 Visualization

- Use `folium`, `QGIS`, or `kepler.gl` to map heat hotspots  
- Overlay population, parks, or vulnerable zones for context

---

## 🧭 Step 3 — Recommend Tree Planting to Reduce UHI

### 🔹 3.1 Cooling Potential Mapping

- Define tree-planting priority zones:
  - High LST
  - Low NDVI
  - Available land cover (non-road, non-building)

- Compute a **“cooling potential score”** per pixel

### 🔹 3.2 Optimization

- Use simple heuristics or greedy search to rank top planting spots  
- Optional: include cost surfaces or buffers for parks, roads

### 🔹 3.3 Visual Outputs

- Maps of recommended planting zones  
- “What-if” visualizations: cooling impact of 1,000 new trees  
- GeoJSON or shapefile of prioritized sites


### **Deliverables:**

* Web Dashboard (Streamlit or Mapbox)
* Visuals: maps, graphs (temperature, NDVI, correlation)
* Brief report or slide deck with findings, methodology, and next steps

---

## 📁 Project Folder Structure
│
├── data/
│ ├── raw/ # ECOSTRESS, Sentinel, GHSL, etc.
│ └── processed/ # Aligned and resampled data
│
├── notebooks/
│ ├── 01_download_data.ipynb
│ ├── 02_preprocessing.ipynb
│ ├── 03_train_model.ipynb
│ ├── 04_predict_highres_LST.ipynb
│ ├── 05_tree_recommendation.ipynb
│
├── models/ # Trained models (.pkl, .joblib)
├── outputs/ # Maps and visual results
│ ├── rasters/
│ └── visualizations/
├── scripts/ # Utility functions
│ └── utils.py
├── requirements.txt
└── README.md


---

## 🚀 Ideas to Extend the Project

### 🔄 1. Generalize to Other Cities

Apply the full pipeline to Paris, Marseille, or international cities.

### 📊 2. Integrate Air Temperature

Use in-situ air temperature (e.g., Netatmo, Météo-France) to calibrate LST and assess human heat exposure.

### 🛠 3. Scenario Simulator

Create a dashboard to test:
- “What if we increase NDVI by 20% here?”
- “What if we remove half the parking lots?”

### 🧑‍⚖️ 4. Add Social Equity Layer

Overlay demographic or income data to target green interventions where people are most vulnerable.

### 🌐 5. Public Dashboard

Build a Streamlit or Mapbox interface to let stakeholders explore UHI patterns and recommendations.

---



