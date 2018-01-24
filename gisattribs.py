from osgeo import gdal
import numpy as np

class GISAttribs(object):
    def __init__(self):
        pass
    
    @classmethod
    def fromFile(cls,fname):
        gd = gdal.Open(fname, gdal.GA_ReadOnly)
        c=cls()
        c.filename=fname
        c.datatype = gdal.GetDataTypeName(gd.GetRasterBand(1).DataType)
        c.ndv = gd.GetRasterBand(1).GetNoDataValue()
        c.xsize = gd.RasterXSize
        c.ysize = gd.RasterYSize
        c.geotransform = gd.GetGeoTransform()
        c.projection = gd.GetProjection()
        c.projectionref = gd.GetProjectionRef()
        c.xscale=np.abs(c.geotransform[1])
        c.yscale=np.abs(c.geotransform[5])
        c.xorigin=np.abs(c.geotransform[0])
        c.yorigin=np.abs(c.geotransform[3])
        return c

    @classmethod
    def fromScenario(cls,s):
        c=cls()
        v=s._v_attrs
        if v.__contains__("filename"):
            c.filename=v.filename
        else:
            c.filename=v.rasterName
        if v.__contains__("datatype"):
            c.datatype=v.datatype
        else:
            c.datatype=None
        if v.__contains__("ndv"):
            c.ndv=v.ndv
        else:
            c.ndv=0
        c.xsize=v.map_x_size
        c.ysize=v.map_y_size
        c.xscale=v.map_x_scale
        c.yscale=v.map_y_scale
        c.xorigin=v.map_x_origin
        c.yorigin=v.map_y_origin
        if v.__contains__("geotransform"):
            c.geotransform=v.geotransform
        else:
            c.geotransform=None
        if v.__contains__("projection"):
            c.projection=v.projection
        else:
            c.projection=None
        if v.__contains__("projectionref"):
            c.projectionref=v.projectionref
        else:
            c.projectionref=None
        return c
    
    def write(self,s):
        v=s._v_attrs
        v.filename=self.filename
        v.datatype=self.datatype
        v.ndv=self.ndv
        v.map_x_size=self.xsize
        v.map_y_size=self.ysize
        v.map_x_scale=self.xscale
        v.map_y_scale=self.yscale
        v.map_x_origin=self.xorigin
        v.map_y_origin=self.yorigin
        v.geotransform=self.geotransform
        v.projection=self.projection
        v.projectionref=self.projectionref

    def __repr__(self):
        s="GIS File Attributes:\n"
        s+="Filename: %s" % self.filename + '\n'
        s+="No Data Value: %s" % self.ndv + '\n'
        s+="x size: " + str(self.xsize) + '\n'
        s+="y size: " + str(self.ysize) + '\n'
        return s

