# Sharpening Sentinel-2 with Planetscope
## Description
Image sharpening has been developed to get the best of both spectral and spatial resolution by integrating different image sources to improve the information and quality of multispectral images. Satellites with no Pan band can use other satellites’ high-resolution bands to emulate a Pan band, provided the pixel resolution and the sharpening algorithm are suited to the specific process. In this context, Li et al. (2020) tested methods to combine Sentinel-2 with high-resolution Planetscope-0 imagery for Earth Observation studies. 

In this repository, we provide an automated and slightly modified version of the RGB-NIR sharpening and SWIR-1(2) synthesizing routine proposed by Li et al. (2020). This is a significant upgrade that makes it possible to easily and quickly sharpen many images or time series and removes the human factor, which can introduce errors. The algorithm workflow is presented in the diagram below.

![diagram](workflow_diagram.jpeg)

## Repository Structure
```
Sharpening-S2_Planetscope
│── README.md
│── requirements.txt
│── workflow_diagram.jpeg
├── code/
│   ├── __init__.py
│   ├── sharpen_classes.py
├── data/
│   ├── pan_test.tif
│   ├── target_swir_test.tif
│   └── target_vnir_test.tif
└── tests/
│   ├── sharpening_test.ipynb
```
## How to cite
Tayer T.C., Douglas M.M., Cordeiro M.C.R., Tayer A.D.N., Callow J.N., Beesley L. & McFarlane D. (2023) Improving the accuracy of the Water Detect algorithm using Sentinel-2, Planetscope and sharpened imagery: a case study in an intermittent river, GIScience & Remote Sensing, 60:1, DOI: 10.1080/15481603.2023.2168676

## Tutorial


## Supported Formats
-- Pan 
* Planetscope - 4 bands (B, G, R, NIR)

-- Target

If only sharpening 4 bands:
* Sentinel-2 - 4 bands (B, G, R, NIR)

If sharpening and synthesizing (SWIR-1 or SWIR-2):
* Sentinel-2 - 4 bands (B,G,R,NIR) + Sentinel-2 - 2 bands (SWIR-1, SWIR-2)
* Sentinel-2 - 6 bands (B, G, R, NIR, SWIR-1, SWIR-2).
## Dependencies
The required libraries are (conda install):
```
arosics>=1.2.6
geoarray>=0.12.3
geopandas>=0.10.2
ipython>=7.29.0
numpy>=1.21.2
pandas>=1.2.2
rasterio>=1.2.0
rioxarray>=0.4.0
scikit-learn>=1.0.1
shapely>=1.7.1
xarray>=0.17.0
opencv>=4.5.3
```
Note that GDAL is required to run Rasterio and others. Thus, we highly recommend using conda to install packages either from:

1. Creating a new conda virtual environment with requirements using:

`conda create --name YourVenvName python=x.x anaconda --file .\requirements.txt`

Or 2. installing requirements in your current environment using conda:

`conda install .\requirements.txt`

## References

Li, Z.; Zhang, H.K.; Roy, D.P.; Yan, L.; Huang, H. Sharpening the Sentinel-2 10 and 20 m Bands to Planetscope-0 3 m Resolution. Remote Sens. 2020, 12, 2406. https://doi.org/10.3390/rs12152406

Tayer T.C., Douglas M.M., Cordeiro M.C.R., Tayer A.D.N., Callow J.N., Beesley L. & McFarlane D. (2023) Improving the accuracy of the Water Detect algorithm using Sentinel-2, Planetscope and sharpened imagery: a case study in an intermittent river, GIScience & Remote Sensing, 60:1, DOI: 10.1080/15481603.2023.2168676
