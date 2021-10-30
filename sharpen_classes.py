import os
import geopandas as gpd
import json
import rasterio
from rasterio.mask import mask
from rasterio.transform import from_origin
from shapely.geometry import box, mapping
import shapely
from geoarray import GeoArray
from arosics import COREG
import numpy as np
import math
from itertools import product
import pandas as pd
from sklearn import linear_model
import xarray as xr
import rioxarray as rxr
import cv2 as cv
import numpy as np


def stack_sharpen_vnir_swir (sharp_vnir, sharp_swir, out_path_name):

    minx, miny, maxx, maxy = sharp_vnir.rio.bounds()
    sharp_vnir_bbox = box(minx, miny, maxx, maxy)
    gdf= gpd.GeoDataFrame({'geometry': [sharp_vnir_bbox]}, crs= sharp_vnir.rio.crs)
    sharp_swir = sharp_swir.rio.clip(gdf.geometry.apply(mapping), gdf.crs)

    minx, miny, maxx, maxy = sharp_swir.rio.bounds()
    sharp_swir_bbox = box(minx, miny, maxx, maxy)
    gdf= gpd.GeoDataFrame({'geometry': [sharp_swir_bbox]}, crs= sharp_vnir.rio.crs)
    sharp_vnir = sharp_vnir.rio.clip(gdf.geometry.apply(mapping), gdf.crs)
    
    stack = xr.concat([sharp_vnir, sharp_swir], dim = 'band')
    stack.rio.to_raster(out_path_name, compress='lzw')
    
    return stack

class Sharpen():
    
    def __init__(self, ms_pan, ms_target, out_path, synth=False):
        self.ms_pan = ms_pan
        self.ms_target = ms_target
        self.out_path = out_path
        self.synth = synth
        
        _ms_pan = rxr.open_rasterio(ms_pan, chunks= (1,500, 500))
        _ms_target = rxr.open_rasterio(ms_target, chunks= (1,500, 500))
        if synth == True:
            target_transform= _ms_target.rio.transform()
            assert target_transform[0] == 20 and _ms_target.shape[0] == 2, 'Please provide S2-SWIR 1 and SWIR 2 stacked image or change synth parameter to False'
        else:
            assert _ms_pan.shape[0] == _ms_target.shape[0], 'Both rasters need to have VNIR bands'
        
        scratch_path = os.path.join(self.out_path, 'scratch_path')
        try:
            os.mkdir(scratch_path)
        except:
            pass
        self.scratch_path = scratch_path

    def check_projection_and_reproject(self): 
        """Checks if the inputs images have equal projection system. If not, the target image will be reprojected."""
        #open input images as dask arrays
        _ms_pan = rxr.open_rasterio(self.ms_pan, chunks = (1,500, 500))
        _ms_target = rxr.open_rasterio(self.ms_target, chunks = (1,500, 500))
        #check if projections are equal 
        if _ms_pan.rio.crs == _ms_target.rio.crs:
            return self.ms_target
        #if not, reproject target using pan
        else:
            _ms_target_name = os.path.basename(self.ms_target)
            ms_target = os.path.join(self.scratch_path, _ms_target_name[:-4] + '_proj.tif')
            _ms_target_affine = _ms_target.rio.transform()
            _ms_target = _ms_target.rio.reproject(_ms_pan.rio.crs, resolution =_ms_target_affine[0])
            _ms_target = _ms_target.where(_ms_target != -9999, drop=True).astype(dtype='uint16')
            _ms_target = _ms_target.rio.write_nodata(0)
            _ms_target.rio.to_raster(ms_target, compress='lzw')
            return ms_target

    def check_pan_extent (self):
        """Checks if the pan image footprint is within the target image extent. If the pan image footprint is larger than the target, the pan image will be clipped to fit the target extent"""
        #input images
        ms_pan = rasterio.open(self.ms_pan, chunks= (1,500, 500))
        ms_target = self.check_projection_and_reproject()
        ms_target_r = rasterio.open(ms_target, chunks= (1,500, 500))
        
        #retrieve pan and target boundaries
        minx, miny, maxx, maxy = ms_pan.bounds
        ms_pan_bbox = box(minx, miny, maxx, maxy)
        minx, miny, maxx, maxy = ms_target_r.bounds
        ms_target_r_bbox = box(minx, miny, maxx, maxy)
        
        #check if pan is larger than target. Else clip pan.
        diff_area= ms_pan_bbox.difference(ms_target_r_bbox)

        if diff_area.area == 0:
            return self.ms_pan, ms_target
        else:
            geo_clip = gpd.GeoDataFrame({'geometry': ms_target_r_bbox}, index=[0], crs= ms_target_r.crs)
            coords = [json.loads(geo_clip.to_json())['features'][0]['geometry']]
            out_img, out_transform = mask(ms_pan, coords, crop=True)
            out_meta = ms_pan.meta.copy()
            out_meta.update({"driver": "GTiff",
                             "height": out_img.shape[1],
                             "width": out_img.shape[2],
                             "transform": out_transform,
                             "crs": ms_pan.crs})
            raster_name= os.path.basename(self.ms_pan)[:-4]
            if self.synth == True:
                ms_pan= os.path.join(self.scratch_path, raster_name + '_clip_to_target_synth.tif')
            else:
                ms_pan= os.path.join(self.scratch_path, raster_name + '_clip_to_target.tif')
            
            #Export clipped image to scratch path
            with rasterio.open(ms_pan, "w", **out_meta) as dest:
                dest.write(out_img)
                dest.close()
            return ms_pan, ms_target
        
    def adjust_and_coregister_ms_target (self):
        """Uses AROSICS to coregister and resample target image and adjust arrays to same number of rows and columns"""
        ms_pan, ms_target = self.check_pan_extent()
        ms_pan = xr.open_rasterio(ms_pan, chunks= (1,500, 500))
        ms_target = xr.open_rasterio(ms_target, chunks= (1,500, 500))
        
        #retrieve pan image transform values
        pan_transform = [ms_pan.transform[2], ms_pan.transform[0], ms_pan.transform[1], ms_pan.transform[5], ms_pan.transform[3], ms_pan.transform[4]]
        pan_projection = ms_pan.rio.crs.wkt
        
        #retrieve and buffer pan boundaries
        bbox = minx, miny, maxx, maxy= ms_pan.rio.bounds() 
        polygon = shapely.geometry.box(*bbox)
        buffer_poly = polygon.buffer(120, join_style = 2)
        gdf = gpd.GeoDataFrame({'geometry': [buffer_poly]}, crs= ms_pan.rio.crs)
        
        #get and transpose array values
        ms_pan_arr = ms_pan.values
        ms_pan_arr = np.nan_to_num(ms_pan_arr)
        ms_pan_arr= ms_pan_arr.transpose((1, 2, 0))
        
        #convert ms_pan to Geoarray
        pan_garr = GeoArray(ms_pan_arr, pan_transform, pan_projection)
        
        #clip ms_target using buffered pan boundaries
        clipped_ms_target = ms_target.rio.clip(gdf.geometry.apply(mapping), gdf.crs)
        
        #get and transpose array values
        target_transform = [clipped_ms_target.rio.bounds()[0], ms_target.transform[0], ms_target.transform[1], clipped_ms_target.rio.bounds()[3], ms_target.transform[3], ms_target.transform[4]]
        ms_target_arr= clipped_ms_target.values
        ms_target_arr = np.nan_to_num(ms_target_arr)
        ms_target_arr= ms_target_arr.transpose((1, 2, 0))
        
        #convert ms_target to Geoarray
        target_garr =  GeoArray(ms_target_arr, target_transform, pan_projection)
        
        #set AROSICS parameters
        kwargs = {
                'ws'  : (64,64),
                'align_grids'  : True,
                'match_gsd'    : True,
                'mask_baddata_ref' : None,
                 'q'            : True,
                'nodata': [0,0],
                }
        #run coregistration
        CRL = COREG(pan_garr, target_garr,**kwargs)
        result= CRL.correct_shifts()
        
        ## result to xarray
        #retrieve raster corner coordinates and pixel resolution
        res = result['updated geotransform'][1]
        xmin = result['updated geotransform'][0]
        ymax = result['updated geotransform'][3]

        #store x and y size (horizontal and vertical array length)
        xsize = result['arr_shifted'].shape[0]
        ysize = result['arr_shifted'].shape[1]

        #get num of bands
        bands= [x+1 for x in range(result['arr_shifted'].shape[2])]

        #create an array that stores all possible coordinates
        x = np.arange(xmin, xmin + xsize * res, res)
        y = np.arange(ymax, ymax - ysize * res, -res)

        #create meshgrids
        x_mesh, y_mesh = np.meshgrid(x, y, indexing='ij')

        #data into dataset
        coreg_target = xr.DataArray(
            data= result['arr_shifted'],
            dims= ['y', 'x', 'band'],
            coords=dict(
                y_m=(["y", "x"], y_mesh),
                x_m=(["y", "x"], x_mesh),
                band=bands,
            )
        )

        coreg_target=coreg_target.transpose('band', 'y', 'x')
        coreg_target.attrs['crs']= ms_pan.attrs['crs']
        coreg_target.attrs['transform'] = (result['updated geotransform'][1], result['updated geotransform'][2],\
                                            result['updated geotransform'][0], result['updated geotransform'][4],\
                                            result['updated geotransform'][5], result['updated geotransform'][3])
        
        ## adjust arrays to same number of rows and columns
        coreg_target= coreg_target.chunk(chunks= (1,500, 500))
        minx, miny, maxx, maxy = ms_pan.rio.bounds()
        ms_pan_bbox = box(minx, miny, maxx, maxy)
        gdf= gpd.GeoDataFrame({'geometry': [ms_pan_bbox]}, crs= ms_pan.rio.crs)
        coreg_target = coreg_target.rio.clip(gdf.geometry.apply(mapping), gdf.crs)
        
        minx, miny, maxx, maxy = coreg_target.rio.bounds()
        coreg_target_bbox = box(minx, miny, maxx, maxy)
        gdf= gpd.GeoDataFrame({'geometry': [coreg_target_bbox]}, crs= ms_pan.rio.crs)
        ms_pan = ms_pan.rio.clip(gdf.geometry.apply(mapping), gdf.crs)
        
        #set no data
        for num, band in enumerate(coreg_target):
            band.values = xr.where(ms_pan[num] ==0 , 0, band.values)
        coreg_target.attrs['_FillValue']= ms_pan.attrs['_FillValue']

        return ms_pan, ms_target, coreg_target

    def spatial_degradation_math (self, ms_pan, target_transform, _dict_planet_pan):
        #SRTM_MTF values from ESA reports
        _dict_SRTM_MTF = {'blue': 0.29, 'green': 0.28, 'red': 0.265, 'nir': 0.23, 'swir1': 0.19 , 'swir2': 0.24}
        
        def mtf (_dict_SRTM_MTF, band, target_transform):
            #set f values depending on pixel spatial res
            if target_transform[0] == 10:
                f = 1/20
            elif target_transform[0] == 20:
                f = 1/40
            numerad = (f) **2
            denumerad = 2* math.log(_dict_SRTM_MTF[band])
            S = math.sqrt(-(numerad/denumerad))
            
            return S
        #convolution matrix
        dict_siglam = {}
        _dict_planet_degraded = {}
        for band in _dict_planet_pan:
            S = mtf (_dict_SRTM_MTF, band, target_transform)
            dict_siglam[band] = S
            sig = 1/((dict_siglam[band])*2*math.pi)

            #Generate 41x41 Matrix
            amin = -20
            amax = 21
            list_m = list(product(range(amin, amax), repeat=2))

            #Matrix to df
            df= pd.DataFrame(list_m)
            df.rename(columns={0: 'i', 1:'j'}, inplace=True)
            df.sort_values(['i','j'], ascending=[True, False], inplace=True)

            df['cst'] = ( 1 /(2*math.pi*(sig/ms_pan.rio.transform()[0])**2))
            df['numerad'] = df['i']**2 + df['j']**2
            df['denum'] = 2*(sig/ms_pan.rio.transform()[0])**2
            df['exp'] = np.exp(-(df['numerad']/df['denum']))
            df['result'] = df['cst']*df['exp']
            #normalize matrix
            df['norm'] = df['result'] /(df['result'].sum())

            #define matrix array
            x = np.array(list(df['norm'])).reshape(41, 41)

            #define kernel for convolution
            kernel = 1/41 * x
            conv = cv.filter2D(_dict_planet_pan[band], -1, kernel)
            mask = (_dict_planet_pan[band] > 0)
            if np.any(~mask):
                scaling_vals = cv.filter2D(mask.astype(np.float32), -1, kernel)
                conv = conv.astype(np.float32)
                conv[mask] /= scaling_vals[mask]
                conv[~mask] = 0
                scaling_vals = None
                mask = None

            _dict_planet_degraded[band] = conv
        
        #degraded values to xarray
        array= []
        for band in _dict_planet_degraded:
            array.append(np.array(_dict_planet_degraded[band]))
        array= np.nan_to_num(array)
        degraded_pan= ms_pan.copy()
        degraded_pan.values= array

        degraded_pan= degraded_pan.chunk(chunks= (1,500, 500))
        return degraded_pan
    
    def spatial_degradation (self):
        """Spatial degradation based on MTF values"""
        
        ms_pan, ms_target, coreg_target = self.adjust_and_coregister_ms_target()
        target_transform= ms_target.rio.transform()
        _dict_planet_pan = {'blue': ms_pan.values[0], 'green': ms_pan.values[1], 'red': ms_pan.values[2], 'nir': ms_pan.values[3]}
        degraded_pan = self.spatial_degradation_math(ms_pan, target_transform, _dict_planet_pan)
        
        return ms_pan, ms_target, coreg_target, degraded_pan

    def synthetic_and_degraded_swir (self):
        """Synthetic SWIR 1 and 2 and Spatial degradation"""
        
        ms_pan, ms_target, coreg_target, degraded_pan = self.spatial_degradation()
        
        df_deg= pd.DataFrame()
        df_res= pd.DataFrame()
        df_target= pd.DataFrame()
        
        #flatten array values and compute linear regression
        array= []
        for num, band in enumerate(degraded_pan, start= 1):
            df_deg[num] = band.values.flatten()
        for num, band in enumerate(ms_pan, start= 1):
            df_res[num] = band.values.flatten()
        for num, band in enumerate(coreg_target):
            df_deg[num+5] = coreg_target.values[num].flatten()
            x = df_deg[[1,2,3,4]]
            y = df_deg[num+5]
            regr = linear_model.LinearRegression()
            regr.fit(x, y)
            coef_arr = regr.coef_
            coef_list = coef_arr.tolist()
            df_res[num+5] = coef_list[0]* df_res[1]  +  coef_list[1]* df_res[2]  +  coef_list[2]* df_res[3]  + coef_list[3]* df_res[4]
            array_synth = df_res[num+5].to_numpy().reshape(coreg_target[num].shape)
            array.append(array_synth)

        array= np.nan_to_num(array)
        synth_pan= coreg_target.copy()
        synth_pan.values= array
        
        #degrade synthetic band
        synth_pan= synth_pan.chunk(chunks= (1,500, 500))
        _dict_planet_pan = {'swir1': synth_pan.values[0], 'swir2': synth_pan.values[1]}
        target_transform= ms_target.rio.transform()
        degraded_synth = self.spatial_degradation_math(synth_pan, target_transform, _dict_planet_pan)
        
        return synth_pan, coreg_target, degraded_synth
    
    def hpf_sharpening (self):
        
        if self.synth == True:
            synth_pan, coreg_target, degraded_synth = self.synthetic_and_degraded_swir()
            sharp = (coreg_target/degraded_synth)*synth_pan
        else:
            ms_pan, ms_target, coreg_target, degraded_pan = self.spatial_degradation()
            sharp = (coreg_target/degraded_pan)*ms_pan

        sharp= sharp.fillna(0)
        sharp = xr.where(sharp < 0, 0, sharp)
        sharp.attrs['_FillValue']= 0
   
        return sharp