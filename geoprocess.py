from osgeo import osr, gdal, gdalconst
import numpy as np
import wradlib as wrl
from pyproj import Proj


def gdal_to_array(datasets):

    """converts gdal.dataset to numpy array"""
    arrs = []
    for ds in datasets:
        arr = np.array(ds.GetRasterBand(1).ReadAsArray())
        arr = np.expand_dims(arr, axis=0)
        arrs.append(arr)
    return arrs

def create_latlon(f_reference):
    """This function creates a lon and lat array for the netCDF variables"""

    p = Proj(proj='utm', zone=33, ellps='WGS84')
    f_reference = gdal.Open(f_reference, gdalconst.GA_ReadOnly)
    header = f_reference.GetGeoTransform()

    ncols = f_reference.RasterXSize
    nrows = f_reference.RasterYSize

    xx = np.linspace(header[0] + (header[1] / 2), header[0] + ncols * header[1] - (header[1] / 2), ncols)
    yy = np.linspace(header[3] + (header[5] / 2), header[3] + nrows * header[5] - (header[5] / 2), nrows)

    a, b = np.meshgrid(xx, yy)
    easting, northing = xx, yy
    longrid, latgrid = p(a, b, inverse=True)

    return easting, northing, latgrid, longrid