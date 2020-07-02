# mhm-preproc

This is a collection of Python 3 scripts for preprocessing datafor use in mesoscale Hydrologic Model

modis_lai.py

This script processes Leaf Area Index (LAI) data by MODIS (MOD15A2H)  to a netCDF file.

LAI data (500 m) by MODIS is provided in 8-daily HDF files. You can use this script to clip and reproject the raw data and merge it into a netCDF file.

Run the main() function of the script and give the directory containing the raw HDF files (data_dir) and a reference raster (reference_raster) as arguments.

For each calendar year of data your data directory should have seperate sub directories. 
The data directory should follow this structure:

    LAI/LAI_2007
    LAI/LAI_2008
    LAI/LAI_2009
    ...

The reference raster should be a geoTIFF, the script will extract the geographical information like extent and cell width and will apply it on the output LAI data. For this task bilinear resampling is used, which should may be changed to nearest neighbour method.


-------------------------------------------------------------------------------------------------------------------------------------------------


soilgrids_lut.py

This script generates a Lookup-Table soil_classedefinition.txt based on a dbf table derived by k-mean-clustering soilgrids data


-------------------------------------------------------------------------------------------------------------------------------------------------

References

[1] Myneni, R., Knyazikhin, Y., Park, T. (2015). MOD15A2H MODIS/Terra Leaf Area Index/FPAR 8-Day L4 Global 500m SIN Grid V006 [Data set]. NASA EOSDIS Land Processes DAAC. Accessed 2020-01-22 from https://doi.org/10.5067/MODIS/MOD15A2H.006

