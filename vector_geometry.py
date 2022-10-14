# Uses code by A.Minaev

import os
from osgeo import gdal, ogr, osr
# from geodata import changeXY
from .temp_files import TempName, StopFromStorage, DeleteFile


def GetImageBoundaryBox(image_ds):
    global xmax, ymax, ymin, xmin, image_bbox
    try:
        image_geotransform = image_ds.GetGeoTransform()
        if image_geotransform != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
            # get resolution
            xsize = image_ds.RasterXSize
            ysize = image_ds.RasterYSize
            # get coordinates of boundary box
            xmin = image_geotransform[0]
            ymin = image_geotransform[3] + xsize * image_geotransform[4] + ysize * image_geotransform[5]
            xmax = image_geotransform[0] + xsize * image_geotransform[1] + ysize * image_geotransform[2]
            ymax = image_geotransform[3]

        else:
            # extract gcps and their coordinates
            points = image_ds.GetGCPs()
            if points:
                points_x = [point.GCPX for point in points]
                points_y = [point.GCPY for point in points]

                # get coordinates of boundary box
                xmin = min(points_x)
                ymin = min(points_y)
                xmax = max(points_x)
                ymax = max(points_y)

        # create a polygon with boundary box
        image_bbox_ring = ogr.Geometry(ogr.wkbLinearRing)
        image_bbox_ring.AddPoint(xmax, ymin)
        image_bbox_ring.AddPoint(xmax, ymax)
        image_bbox_ring.AddPoint(xmin, ymax)
        image_bbox_ring.AddPoint(xmin, ymin)
        image_bbox_ring.AddPoint(xmax, ymin)

        image_bbox = ogr.Geometry(ogr.wkbPolygon)
        image_bbox.AddGeometry(image_bbox_ring)

    except:
        pass

    return image_bbox


@StopFromStorage
def GetVectorLayerGeometry(vpath):
    vector_ds = ogr.Open(vpath, 1)
    vector_layer = vector_ds.GetLayer()

    # creating multipolygon for features in vector layer
    vector_geometry = ogr.Geometry(ogr.wkbMultiPolygon)

    for i in range(vector_layer.GetFeatureCount()):
        feature = vector_layer.GetFeature(i)
        feature_geometry = feature.GetGeometryRef()
        vector_geometry.AddGeometry(feature_geometry)

    return vector_geometry

# def VectorizeRasterBand(rpath, vpath):
#     cmd = f'python3 gdal_polygonize.py  -q -b mask "{rpath}" "{vpath}"'
#     # print(cmd)
#     os.system(cmd)
#     vector_ds = ogr.Open(vpath, 1)
#     vector_layer = vector_ds.GetLayer()
#     for feat in vector_layer:
#         if feat.GetField('DN') == 0:
#             vector_layer.DeleteFeature(feat.GetFID())
#     vector_ds = None


@StopFromStorage
def rasterDataCover(rpath, vpath):
    # mapping between gdal type and ogr field type
    type_mapping = {gdal.GDT_Byte: ogr.OFTInteger,
                    gdal.GDT_UInt16: ogr.OFTInteger,
                    gdal.GDT_Int16: ogr.OFTInteger,
                    gdal.GDT_UInt32: ogr.OFTInteger,
                    gdal.GDT_Int32: ogr.OFTInteger,
                    gdal.GDT_Float32: ogr.OFTReal,
                    gdal.GDT_Float64: ogr.OFTReal,
                    gdal.GDT_CInt16: ogr.OFTInteger,
                    gdal.GDT_CInt32: ogr.OFTInteger,
                    gdal.GDT_CFloat32: ogr.OFTReal,
                    gdal.GDT_CFloat64: ogr.OFTReal}

    ds = gdal.Open(rpath)
    prj = ds.GetProjection()
    src_band = ds.GetRasterBand(1)
    mask_band = src_band.GetMaskBand()
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(vpath)
    srs = osr.SpatialReference(wkt=prj)

    dst_layer = dst_ds.CreateLayer('', srs=srs)
    raster_field = ogr.FieldDefn('id', type_mapping[src_band.DataType])
    dst_layer.CreateField(raster_field)
    options = ['DATASET_FOR_GEOREF=' + rpath]
    gdal.Polygonize(mask_band, mask_band, dst_layer, 0, options, callback=None)
    del rpath, ds, mask_band, src_band, dst_ds, dst_layer


# Can be replaced by rasterCoverGeometry()
@StopFromStorage
def GetRasterCover(rpath):
    vpath = TempName('tmp', 'shp')
    # VectorizeRasterBand(rpath, vpath)
    rasterDataCover(rpath, vpath)
    return GetVectorLayerGeometry(vpath)


# !!! Still has problems in axis order while srs transformations -- need to solve
#checking intersection of image and shapefile
def GetIntersection(rpath, vpath, target_srs = None):

    #get extent of layers
    # image_bbox = GetImageBoundaryBox(image_ds)
    raster_geometry = GetRasterCover(rpath)
    vector_geometry = GetVectorLayerGeometry(vpath)

    image_ds = gdal.Open(rpath)
    image_srs = image_ds.GetSpatialRef()
    vector_ds = ogr.Open(vpath)
    vector_srs = vector_ds.GetLayer().GetSpatialRef()
    if target_srs is None:
        target_srs = image_srs

    if image_srs.ExportToProj4()!=vector_srs.ExportToProj4():
        raster_geometry.Transform(osr.CoordinateTransformation(image_srs, vector_srs))
    if not raster_geometry.Intersects(vector_geometry):
        invert = True
        raster_geometry = changeXY(raster_geometry)
    else:
        invert = False

    #get square of intersection between image and vector layer
    intersection = raster_geometry.Intersection(vector_geometry)

    if intersection is not None:
        if target_srs.ExportToProj4()!=vector_srs.ExportToProj4():
            intersection.Transform(osr.CoordinateTransformation(vector_srs, target_srs))
        if invert:
            intersection = changeXY(intersection)

    return intersection


@StopFromStorage
def rasterCoverGeometry(raster_path, raster_cover_path=None, srcnodata=None, spatial_reference=None):

    if raster_cover_path is None:
        raster_cover_path = TempName(ext='shp')

    if srcnodata is not None:
        rasterDataCoverSlow(raster_path, raster_cover_path, no_data_value=srcnodata)
    else:
        rasterDataCover(raster_path, vector_path=raster_cover_path)

    cover_geometry = vectorLayerGeometry(ogr.Open(raster_cover_path), geometryType=ogr.wkbMultiPolygon)

    if spatial_reference is not None:
        cover_geometry = geometryCoordinateTransformation(cover_geometry,
                                                          gdal.Open(raster_path).GetSpatialRef(),
                                                          spatial_reference)

    return cover_geometry


@StopFromStorage
def vectorLayerGeometry(data_source, iLayer=0,
                        geometryType=ogr.wkbMultiPolygon,
                        srs=None):
    layer = data_source.GetLayer(iLayer)
    if layer is not None:
        vector_geometry = ogr.Geometry(geometryType)
        for feature in layer:
            feature_geometry = feature.GetGeometryRef()
            if feature_geometry.GetGeometryType() == ogr.wkbMultiPolygon:
                vector_geometry = vector_geometry.Union(feature_geometry)
            else:
                vector_geometry.AddGeometry(feature_geometry)
        if srs is not None:
            geometryCoordinateTransformation(
                vector_geometry, layer.GetSpatialRef(), srs)
        return vector_geometry


def geometryCoordinateTransformation(geometry, start_spatial_ref, end_spatial_ref):
    if (start_spatial_ref is not None) and (end_spatial_ref is not None):
        if start_spatial_ref.ExportToProj4() != end_spatial_ref.ExportToProj4():
            geometry.Transform(osr.CoordinateTransformation(start_spatial_ref, end_spatial_ref))
    return geometry


def geometryExtent(geometry):
    xmin, xmax, ymin, ymax = geometry.GetEnvelope()
    return [xmin, ymin, xmax, ymax]


def geometryShapefile(geometry, vector_path, spatial_reference=None):
    data_source = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(vector_path)
    if data_source is not None:
        layer = data_source.CreateLayer('', srs=spatial_reference, geom_type=geometry.GetGeometryType())
        feature = ogr.Feature(ogr.FeatureDefn())
        feature.SetGeometry(geometry)
        layer.CreateFeature(feature)
        data_source = None
        with open(vector_path.replace(os.path.splitext(vector_path)[-1], '.prj'), 'w') as prj:
            prj.write(spatial_reference.ExportToWkt())


# Contains data for clipping raster with vector
class VectorClipper(object):

    def __init__(self, vector_path):
        self.path = vector_path
        self.data_source = ogr.Open(vector_path)
        if self.data_source is None:
            self.layer = None
            self.spatial_reference = None
            self.geometry = None
        else:
            self.layer = self.data_source.GetLayer()
            self.spatial_reference = self.layer.GetSpatialRef()
            self.geometry = vectorLayerGeometry(self.data_source, geometryType=ogr.wkbMultiPolygon)

    def __str__(self):
        if self.spatial_reference is None:
            return f'VectorClipper: "{self.path}" without spatial reference data'
        else:
            return f'VectorClipper: {self.path} in {self.spatial_reference.ExportToProj4()}'

    def __repr__(self):
        return self.__str__()

    def intersect(self, check_geometry):
        if self.geometry.Intersects(check_geometry):
            if check_geometry.Within(self.geometry):
                return 2
            elif self.geometry.Within(check_geometry):
                return 3
            else:
                return 1
        else:
            return 0

    def intersection(self, check_geometry):
        intersect_status = self.intersect(check_geometry)
        if intersect_status == 0:
            return None
        elif intersect_status == 1:
            return self.geometry.Intersection(check_geometry)
        elif intersect_status == 2:
            return check_geometry
        elif intersect_status == 3:
            return self.geometry
        else:
            raise Exception(f'Wrong intersect status {intersect_status}, a value between 0 and 3 is required')

    def clipRasterParameters(self, raster_path, make_raster_data_cover=True, end_spatial_ref=None, srcnodata=None):
        # print(raster_path)
        if make_raster_data_cover:
            print('making_raster_cover')
            raster_cover_path = TempName('tmp', 'shp')
            if srcnodata is not None:
                print('UsingSlow')
                rasterDataCoverSlow(raster_path, raster_cover_path, no_data_value=srcnodata)
            else:
                # rasterDataCover(raster_path, vector_path=raster_cover_path)
                rasterDataCover(raster_path, raster_cover_path)
            cover_geometry = vectorLayerGeometry(ogr.Open(raster_cover_path), geometryType=ogr.wkbMultiPolygon)
        else:
            cover_geometry = GetImageBoundaryBox(gdal.Open(raster_path))

        cover_geometry = geometryCoordinateTransformation(cover_geometry,
                                                          gdal.Open(raster_path).GetSpatialRef(),
                                                          self.spatial_reference)
        # print(self.geometry.ExportToWkt(), cover_geometry.GetEnvelope(), self.spatial_reference.ExportToProj4())
        intersect_status = self.intersect(cover_geometry)

        if intersect_status == 0:
            # No intesection
            raise Exception(f'{self.path} and {raster_path} areas do not intersect, cannot clip')

        elif intersect_status == 1:
            # Raster intersects vector
            intersection_geometry = self.geometry.Intersection(cover_geometry)
            if end_spatial_ref is not None:
                intersection_geometry = geometryCoordinateTransformation(intersection_geometry,
                                                                         self.spatial_reference, end_spatial_ref)
            clip_geometry_path = TempName('tmp', 'shp')
            geometryShapefile(intersection_geometry, clip_geometry_path, spatial_reference=self.spatial_reference)
            # print(clip_geometry_path, geometryExtent(intersection_geometry))
            return {'crop_to_cutline': True,
                    'cutline': clip_geometry_path,
                    'Extent': geometryExtent(intersection_geometry)}

        elif intersect_status == 2:
            # Raster within vector
            return {'crop_to_cutline': False,
                    'cutline': None,
                    'Extent': None}

        elif intersect_status == 3:
            # Vector within raster
            return {'crop_to_cutline': True,
                    'cutline': self.path,
                    'Extent': geometryExtent(self.geometry)}

        # Wrong intersect status
        else:
            raise Exception(f'Wrong intersect status {intersect_status}, a value between 0 and 3 is required')


@StopFromStorage
def rasterDataCoverSlow(raster_path, raster_cover_path, band_num=1, no_data_value=None):
    raster = gdal.Open(raster_path)
    band = raster.GetRasterBand(band_num)
    arr = band.ReadAsArray()
    if no_data_value is None:
        no_data_value = band.GetNoDataValue()
    mask_arr = arr != no_data_value

    mask_temp_path = TempName(ext='tif')
    mask_raster = gdal.GetDriverByName('GTiff').Create(mask_temp_path,
                            raster.RasterXSize, raster.RasterYSize, 1, 1)
    mask_raster.SetProjection(raster.GetProjection())
    mask_raster.SetGeoTransform(raster.GetGeoTransform())
    mask_band = mask_raster.GetRasterBand(1)
    mask_band.WriteArray(mask_arr)
    mask_band.SetNoDataValue(0)
    mask_raster = None

    rasterDataCover(mask_temp_path, raster_cover_path)

    DeleteFile(mask_temp_path)


def removeGeometryPits(geometry):

    geometry_type = geometry.GetGeometryType()
    wkt = geometry.ExportToWkt()

    if geometry_type == 3:
        wkt_parts = [wkt]
        ending = '))'
    elif geometry_type == 6:
        wkt_parts = wkt.split(')),((')
        ending = ')))'
    else:
        return geometry

    for i, part in enumerate(wkt_parts):
        wkt_parts[i] = part.split('),(')[0]

    wkt = ')),(('.join(wkt_parts) + ending

    return ogr.Geometry(wkt=wkt)



def removeGeometryPitsFromFile(path):
    data_source = ogr.Open(path, 1)
    layer = data_source.GetLayer()
    for feature in layer:
        geometry = feature.GetGeometryRef()
        feature.SetGeometry(removeGeometryPits(geometry))
        layer.SetFeature(feature)
    data_source = None


def getVectorMask(data_source, array_shape, geotrans):
    layer = data_source.GetLayer()
    y_res, x_res = array_shape
    dataset = gdal.GetDriverByName('MEM').Create('', x_res, y_res, 1, gdal.GDT_Byte)
    dataset.SetGeoTransform(geotrans)
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(0)
    gdal.RasterizeLayer(dataset, [1], layer, burn_values=[1])
    mask = band.ReadAsArray().astype(bool)
    return mask


def clipRasterWithMask(rpath, mpath):
    rds = gdal.Open(rpath, 1)
    raster_cover_path = TempName(ext='shp')
    rasterDataCover(mpath, raster_cover_path)
    removeGeometryPitsFromFile(raster_cover_path)
    mask_arr = getVectorMask(ogr.Open(raster_cover_path),
                             (rds.RasterYSize, rds.RasterXSize),
                             rds.GetGeoTransform())
    for bandnum in range(1, rds.RasterCount + 1):

        rband = rds.GetRasterBand(bandnum)
        new_data = rband.ReadAsArray() * mask_arr  # Raster Data NoDataValue is supposed to be zero
        rband.WriteArray(new_data)