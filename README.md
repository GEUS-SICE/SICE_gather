# SICE_gather

@author: Rasmus Nielsen (rabni@geus.dk) and Jason Box (jeb@geus.dk)

to read from Thredds server at https://thredds.geus.dk/thredds/catalog/SICE_Greenland_500m/catalog.html
do:
    pip install opendap-protocol
    conda update --all
and
    conda install -c conda-forge netcdf4

variables available, can take just what are needed instead of gathering the heavy product

['ANG', 'AOD_550', 'O3_SICE', 'al', 'albedo_bb_planar_sw', 'albedo_bb_spherical_sw', 
 'albedo_spectral_planar_01'.. 'albedo_spectral_planar_21',
 'cloud_mask', 'crs', 'cv1', 'cv2', 'factor', 'grain_diameter', 'isnow', 'lat', 'lon', 'r0', 
 'rBRR_01', .. 'rBRR_21', 
 'r_TOA_01', .. 'r_TOA_21', 'saa', 'snow_specific_surface_area', 'sza',
 'threshold', 'vaa', 'vza']
