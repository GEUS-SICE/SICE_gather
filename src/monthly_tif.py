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


def date_delta(first_day, delta):
    
    """
    Calculate a new date by adding a delta to a given date.

    Args:
        first_day (datetime.date): The starting date.
        delta (int): Number of days to add.

    Returns:
        datetime.date: The resulting date.
    """
    return (first_day + datetime.timedelta(days=delta))


def multimaps(area, month, var, base_folder):
    
    """
    Generate monthly maps for a specified area and variables.

    Args:
        area (str): The area you want to map.
        month (list): List of months to compute in the format "yyyy-mm".
        var (list): List of variables creating monthly maps.
        base_folder (str): Path to the repository.

    This function loads data for the given area and variables, processes it, and exports monthly maps.
    """
    
    for m in month:
        
        print(f"Processing: {m}")
    
        output_folder = os.path.join(base_folder, 'monthlymaps')
    
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        
        first_day = datetime.date(int(m[:4]), int(m[-2:]), 1)
        last_day = datetime.date(int(m[:4]), int(m[-2:]) + 1, 1) - datetime.timedelta(days=1)
        days = [date_delta(first_day, n).strftime("%Y_%m_%d") for n in range(32) if date_delta(first_day, n) <= last_day]
        mm = m.replace('-', '_')
    
        d_ref = days[0]
        ds_ref = xr.open_dataset(f'https://thredds.geus.dk/thredds/dodsC/SICEvEDC_500m/{area}/sice_500_{d_ref}.nc')
        list_of_dsvar = list(ds_ref.keys())
    
        print('Loading files and saving into a matrix')
        for v in var:
            if v in list_of_dsvar:
                start = 0
                for i, d in enumerate(days):
    
                    ref = f'sice_500_{d}.nc'
    
                    ds = xr.open_dataset(f'https://thredds.geus.dk/thredds/dodsC/SICEvEDC_500m/{area}/{ref}')
                    z = np.array(ds[v])
    
                    if start == 0:
                        data = np.tile(z * np.nan, (len(days), 1, 1))
                        start = 1
    
                    data[i, :, :] = z
            else:
                print(f'{v} does not exist as a SICE variable. Please check on the SICE_gather GitHub for variable names.')
    
            print('Merging data...')
            mergemedian = np.nanmedian(data, axis=0)
            name_median = f"{mm}_{v}_sice_monthlymedian.tif"
    
            try:
                x, y = np.meshgrid(np.array(ds_ref[v].x), np.array(ds_ref[v].y))
            except:
                x, y = np.meshgrid(np.array(ds_ref[v].x2), np.array(ds_ref[v].y2))
    
            exporttiff(x, y, mergemedian, CRS.from_string("+init=EPSG:3413"), output_folder, name_median)
    
            print(f"{mm} for {v} has been exported")


def exporttiff(x, y, z, crs, path, filename):
    
    """
    Export data as a GeoTIFF file.

    Args:
        x (ndarray): x grid.
        y (ndarray): y grid.
        z (ndarray): Data parameter.
        crs (CRS): The data projection.
        path (str): Export path.
        filename (str): Name of the TIFF file.
    """
    resx = (x[0, 1] - x[0, 0])
    resy = (y[1, 0] - y[0, 0])
    transform = Affine.translation((x.ravel()[0]), (y.ravel()[0])) * Affine.scale(resx, resy)

    if resx == 0:
        resx = (x[0, 0] - x[1, 0])
        resy = (y[0, 0] - y[0, 1])
        transform = Affine.translation((y.ravel()[0]), (x.ravel()[0])) * Affine.scale(resx, resy)

    with rasterio.open(
    os.path.join(path, filename),
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


if __name__ == "__main__":
    months = ['2018-09']  # List of months to compute, in format "yyyy-mm"
    var = ['grain_diameter']  # List of variables creating monthly maps
    base = os.path.abspath('..')  # Path to repo
    area = 'Greenland'  # Area you want to map
    multimaps(area, months, var, base)