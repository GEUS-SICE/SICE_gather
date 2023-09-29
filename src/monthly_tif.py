# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 11:51:12 2023

@author: rabni
"""

import rasterio
import os
import glob
import numpy as np
import datetime
from pyproj import CRS
import xarray as xr
from rasterio.transform import Affine


def date_delta(first_day,delta):
    return (first_day + datetime.timedelta(days=delta))

def multimaps(area,month,var,base_folder):
    
    #logging.info("Processing: " + month)
    print(f"Processing: {month}")
    
    output_folder = base_folder + os.sep + 'monthlymaps'

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
   
    
    first_day = datetime.date(int(month[:4]), int(month[-2:]), 1)
    last_day = datetime.date(int(month[:4]), int(month[-2:]) + 1, 1) - datetime.timedelta(days=1)
    days = [date_delta(first_day,n).strftime("%Y_%m_%d") for n in range(32) if date_delta(first_day,n)<=last_day]
    month = month.replace('-','_')
    
    d_ref = days[0]
    ds_ref = xr.open_dataset(f'https://thredds.geus.dk/thredds/dodsC/SICEvEDC_500m/{area}/sice_500_{d_ref}.nc')
    list_of_dsvar = list(ds_ref.keys())
    
    print('loading files, and saving into matrix')
    for v in var:
        if v in list_of_dsvar:
            start = 0
            for i,d in enumerate(days): 
                
                ref = f'sice_500_{d}.nc'
                
                ds = xr.open_dataset(f'https://thredds.geus.dk/thredds/dodsC/SICEvEDC_500m/{area}/{ref}')
                z = np.array(ds[v])
                
                if start == 0: 
                    data = np.tile(z * np.nan, (len(days), 1, 1))
                    start = 1
                
                data[i,:,:] = z  
        else:
            print(f'{v} does not exist as a SICE variable, please check on the SICE_gather GitHub for varable names')
            
        print('merging....')
        mergemedian = np.nanmedian(data, axis=0)        
        name_median = f"{month}_{v}_sice_monthlymedian.tif"
        
        try: 
            x,y = np.meshgrid(np.array(ds_ref[v].x),np.array(ds_ref[v].y))
        except: 
            x,y = np.meshgrid(np.array(ds_ref[v].x2),np.array(ds_ref[v].y2))
        
        exporttiff(x,y,mergemedian,CRS.from_string("+init=EPSG:3413"),output_folder,name_median)
        
        print(f"{month} for {v} has been exported")
        #logging.info(f"{month} has been exported")
    

def exporttiff(x,y,z,crs,path,filename):
    
    "Input: xgrid,ygrid, data paramater, the data projection, export path, name of tif file"
    
    resx = (x[0,1] - x[0,0])
    resy = (y[1,0] - y[0,0])
    transform = Affine.translation((x.ravel()[0]),(y.ravel()[0])) * Affine.scale(resx, resy)
    
    if resx == 0:
        resx = (x[0,0] - x[1,0])
        resy = (y[0,0] - y[0,1])
        transform = Affine.translation((y.ravel()[0]),(x.ravel()[0])) * Affine.scale(resx, resy)
    
    with rasterio.open(
    path + os.sep + filename,
    'w',
    driver='GTiff',
    height=z.shape[0],
    width=z.shape[1],
    count=1,
    dtype=z.dtype,
    crs=crs,
    transform=transform,
    ) as dst:
        dst.write(z, 1)
    
    dst.close()
    
    return None 

if __name__ == "__main__":
    
    
    months = ['2018-08'] # list of months to compute, in format "yyyy-mm"
    var = ['grain_diameter']
    base = os.path.abspath('..')
    area = 'Greenland'
    multimaps(area,months,var,base)
    