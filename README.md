# Sharpening Sentinel-2 with Planetscope
## Synopsis
Image sharpening has been developed to get the best of both spectral and spatial resolution by integrating different image sources to improve the information and quality of multispectral images. Satellites with no Pan band can use other satellites’ high-resolution bands to emulate a Pan band, provided the pixel resolution and the sharpening algorithm are suited to the specific process. In this context, Li et al. (2020) tested methods to combine Sentinel-2 with high-resolution Planetscope-0 imagery for Earth Observation studies. 

In this repository we provide an automated and slightly modified version for the RGB-NIR sharpening and SWIR-1(2) sythesizing routine proposed by Li et al. (2020). This is a significant upgrade that makes it possible to easily and quickly sharpen many images or time series and removes the human factor, which can introduce error. The algorithm workflow is presented in the diagram below.

![diagram](workflow_diagram.jpeg)

All the details and tests are described in the article ------ .

## Repository Structure
```
Sharpening-S2_Planetscope
│── README.md
│── requirements.txt
│── workflow_diagram.jpeg
├── code/
│   ├── __init__.py
│   ├── sharpen_classes.py
│   ├── main.py
├── data/
│   ├── pan_test.tif
│   ├── target_swir_test.tif
│   └── target_vnir_test.tif
└── tests/
│   ├── img1/
│   ├── img2/
```
## How to cite
Tayer et al 2022......

## Tutorial


## Supported Formats


## Dependencies
The required libraries are:
```
arosics==1.2.6
geoarray==0.12.3
geopandas==0.10.2
ipython==7.30.0
numpy==1.20.1
pandas==1.2.2
rasterio==1.2.0
rioxarray==0.4.0
scikit-learn==1.0.1
Shapely==1.7.1
xarray==0.17.0
opencv==4.5.3
```

Create a new conda virtual environment with requirements using:

`conda create --name YourVenvName --file .\requirements.txt`

or, install requirements in your current environment using conda:

`conda install .\requirements.txt`
## References

Li, Z.; Zhang, H.K.; Roy, D.P.; Yan, L.; Huang, H. Sharpening the Sentinel-2 10 and 20 m Bands to Planetscope-0 3 m Resolution. Remote Sens. 2020, 12, 2406. https://doi.org/10.3390/rs12152406
