import os
import glob
from datetime import date, datetime, timedelta
import numpy as np
import pandas as pd
import xarray as xr
import geoprocess as gproc
from pyproj import Proj
from osgeo import gdal, gdalconst


def match_raster(f_input, f_reference):
    """This function takes a LAI .hdf file as input and returns a reprojected and resampled gdal dataset fitting the
    geographic parameters of a reference raster"""

    # open the .hdf file and read raster with GDAL
    f_input = gdal.Open(f_input, gdalconst.GA_ReadOnly)
    f_input = gdal.Open(f_input.GetSubDatasets()[1][0], gdal.GA_ReadOnly)

    # define the MODIS sinusoidal projection according to product user guide
    input_proj = 'PROJCS["Sinusoidal", ' \
                 'GEOGCS["GCS_unnamed ellipse",' \
                 'DATUM["D_unknown", ' \
                 'SPHEROID["Unknown",6371007.181,0]],' \
                 'PRIMEM["Greenwich",0], ' \
                 'UNIT["Degree",0.017453292519943295]],' \
                 'PROJECTION["Sinusoidal"], ' \
                 'PARAMETER["central_meridian",0],' \
                 'PARAMETER["false_easting",0], ' \
                 'PARAMETER["false_northing",0],' \
                 'UNIT["Meter",1]]'

    # get attributes of LAI map
    input_params = f_input.GetGeoTransform()
    input_band = f_input.GetRasterBand(1)
    input_x = f_input.RasterXSize
    input_y = f_input.RasterYSize
    input_data = f_input.ReadAsArray()
    input_data = input_data.astype(np.float)

    # set all fill values according to MOD15A2H manual to NaN, multiply data with scale factor
    scale_fact = 0.1
    fill_values = list(range(249, 256))

    for fvalue in fill_values:
        input_data[input_data == fvalue] = np.nan
    input_data = np.where(input_data == np.nan, input_data, (input_data * scale_fact))

    # create new gdal dataset
    driver = gdal.GetDriverByName('MEM')
    input_new = driver.Create('', input_x, input_y, 1, gdal.GDT_Float32)
    input_new.SetGeoTransform(input_params)
    input_new.SetProjection(input_proj)
    input_new.GetRasterBand(1).WriteArray(input_data)
    f_input = input_new

    # reproject MODIS to UTM33N
    reprojected = gdal.Warp('', f_input, srcSRS=input_proj, dstSRS='EPSG:25833', format='VRT',
                  outputType=gdal.GDT_Float32, xRes=500, yRes=500)

    # get meta data from reference raster
    reference = gdal.Open(f_reference, gdalconst.GA_ReadOnly)
    reference_trans = reference.GetGeoTransform()
    reference_proj = reference.GetProjection()
    bandreference = reference.GetRasterBand(1)
    bandreference.SetNoDataValue(np.nan)
    x = reference.RasterXSize
    y = reference.RasterYSize

    # clip to reference raster extent and cell width
    clipped = driver.Create('', x, y, 1, bandreference.DataType)
    clipped.SetGeoTransform(reference_trans)
    clipped.SetProjection(reference_proj)
    gdal.ReprojectImage(reprojected, clipped, reference_proj, reference_proj, gdalconst.GRA_Bilinear)
    del reprojected

    return clipped


def modis_dates(files):
    """This functions reads the julian dates from the .hdf file names and converts them to pandas datetime index
    for time dimension of xarray"""

    timerow = []
    for file in files:
        julian_date = file.split('\\')[1].split('.')[1]
        year, d = julian_date[1:5], julian_date[5:]
        timerow.append(pd.to_datetime(date(int(year), 1, 1) + timedelta(days=int(d) - 1)))
    return timerow


def create_dataset(directory, year, ref):
    """Creates xarray dataset for a year"""

    print(f'Creating dataset for {year}...')
    filelist = sorted(glob.glob(directory + "/LAI_" + year + "/*.hdf"))
    time_dim = modis_dates(filelist)

    ds_year = []
    for file in filelist:
        print(f'processing {file}')
        ds_year.append(match_raster(file, ref))

    # convert gdal dataset to numpy array
    ds_year = gproc.gdal_to_array(ds_year)

    # stack np arrays over time
    ds_init, ds_year = ds_year[0], ds_year[1:]
    for ds in ds_year:
        ds_init = np.append(ds_init, np.atleast_3d(ds), axis=0)

    # flip raster along the vertical axis
    new_ds = np.empty((len(ds_year) + 1, ds_year[0].shape[1], ds_year[0].shape[2]))
    for i in range(0, ds_init.shape[0]):
        new_ds[i][new_ds[i] == 0] = np.nan
        new_ds[i] = np.flip(ds_init[i], axis=0)

    # create xarray with time coord
    lai_array = xr.DataArray(new_ds,
                             dims=['time', 'y', 'x'],
                             coords={'time': time_dim})

    # create new Dataset
    xr_ds = xr.Dataset({'lai': (['time', 'y', 'x'], lai_array)},
                       coords={'time': time_dim})

    # upsample to monthly average, fill na/zero with e-10. This is necessary because mHM does neither accept NaN, -9999
    # nor 0.0 (this is probably because of floating point uncertainty in FORTRAN from e-20)
    xr_ds = xr_ds.resample(time='M').mean()
    xr_ds = xr_ds.fillna(0.0000000001)
    xr_ds = xr_ds.where(xr_ds != 0, 0.0000000001)

    # flip the dataset along the y-axis, don't know why it's upside down in the first place
    # later note: you flipped it in line 118, so what did you expect, a plasma tv?
    xr_ds = xr_ds.sortby('y', ascending=False)
    return xr_ds


def main(data_dir, reference_raster):
    """Takes the data directory containing LAI data and reference raster as input
    gives netCDF file with merged data as output

    Structure of LAI data directory:
    LAI/LAI_2007
    LAI/LAI_2008
    ...
    """
    dirs = glob.glob(data_dir + '/*')
    years = [x[-4:] for x in dirs]
    output = []
    for year in years:
        output.append(create_dataset(data_dir, year, reference_raster))
    print('Merging datasets...')
    xr_ds = xr.merge(output, compat='override')
    xr_ds.lai.attrs = {'long_name': 'LAI from MODIS'}

    xr_ds.time.attrs = {'standard_name': 'time'}

    xr_ds.encoding = {'unlimited_dims': ['time']}

    # set global attributes
    current_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    xr_ds.attrs = {'title': 'LAI from MODIS ',
                   'institution': 'Helmholtz Centre for Environmental Research - UFZ, Dpt. ENVINF, '
                                  'processed by Marco Hannemann',
                   'source': 'MODIS MOD15A2H',
                   'conventions': 'CF-1.7',
                   'references': 'https://www.github.com/marcoHannemann/mhm-preproc',
                   'history': f'{current_dt} Python 3.7'}

    xr_ds.to_netcdf('lai_temp.nc', encoding={'lai': {'dtype': 'double'},
                                              'time': {'dtype': 'int'}})
    print('Reopen netcdf and add attributes...')
    ds = xr.open_dataset('lai_temp.nc')

    # get coordinates from reference raster file
    latlon = create_latlon(reference_raster)
    easting, northing, latgrid, longrid = latlon

    la = ds.lai.values
    dt = ds.time.values

    ds = xr.Dataset({'lai': (['time', 'northing', 'easting'], la),
                     'time': dt,
                     'northing': ('northing', northing),
                     'easting': ('easting', easting),
                     'lat': (['northing', 'easting'], latgrid),
                     'lon': (['northing', 'easting'], longrid)})

    ds.lai.attrs = {'long_name': 'monthly LAI',
                    '_FillValue': -9999.0,
                    'coordinates': 'lon lat'}

    ds.easting.attrs = {'long_name': 'x-coordinate in cartesian coordinates UTM33N',
                        '_FillValue': False,
                        'units': 'm'}

    ds.northing.attrs = {'long_name': 'y-coordinate in cartesian coordinates UTM33N',
                         '_FillValue': False,
                         'units': 'm'}

    ds.lat.attrs = {'long_name': 'latitude',
                    '_FillValue': False,
                    'units': 'degrees_north'}

    ds.lon.attrs = {'long_name': 'longitude',
                    '_FillValue': False,
                    'units': 'degrees_east'}

    ds.encoding = {'unlimited_dims': ['time']}
    current_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # set global attributes
    ds.attrs = {'title': 'LAI from MODIS',
                'institution': 'Helmholtz Centre for Environmental Research - UFZ, Dpt. ENVINF'
                'processed by Marco Hannemann',
                'source': 'radiospectrometer',
                'conventions': 'CF-1.7',
                'references': 'https://www.github.com/marcoHannemann/mhm-preproc',
                'history': '{} Python 3.7'.format(current_dt)}

    ds = ds.fillna(0.0000000001)
    print('Creating netCDF...')
    ds.to_netcdf('lai_250.nc',
                 encoding={'time': {'dtype': 'int'},
                           'lai': {'dtype': 'double'}},
                 unlimited_dims='time')

	ds.close()
	os.remove('lai_temp.nc')
main('data/LAI/', 'data/grid250.tif')
