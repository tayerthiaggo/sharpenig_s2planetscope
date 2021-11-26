# Sharpening Sentinel-2 with Planetscope
## Synopsis
Image fusion, also called image sharpening or multisensory merging, has been developed to get the best of both spectral and spatial resolution by integrating different image sources to improve the information and quality of multispectral images. Satellites with no Pan band can use other satellites’ high-resolution bands to emulate a Pan band, provided the pixel resolution and the sharpening algorithm are suited to the specific process. For example, to take advantage of Planetscope’s four finest spectral resolution bands, Li et al. (2020) tested methods to combine spectrally-corrected Sentinel-2 imagery with high-resolution but uncorrected Planetscope-0 imagery for Earth Observation studies. Li et al. (2020) conclude that it is feasible to generate 3m VNIR (Visible – Blue/Green/Red, and Near-Infrared), red-edge, and SWIR reflectance on days that Sentinel-2 and Planetscope-0 are spatially overlapping, which is most beneficial for image classification approaches, and, hence, opens the possibility to the application of water detection algorithms. 

The automatization of the sharpening routine proposed by Li et al. (2020) is a significant upgrade that makes it possible to easily and quickly sharpen many images or time series and removes the human factor, which can introduce error. 

This study sharpened a series of multispectral Sentinel-2 with Planetscope imagery by automating and slightly modifying the algorithm suggested by Li et al. (2020).


## How to cite


## Tutorial


## Supported Formats


## Dependencies

```


```

## Reference

Li, Z.; Zhang, H.K.; Roy, D.P.; Yan, L.; Huang, H. Sharpening the Sentinel-2 10 and 20 m Bands to Planetscope-0 3 m Resolution. Remote Sens. 2020, 12, 2406. https://doi.org/10.3390/rs12152406
