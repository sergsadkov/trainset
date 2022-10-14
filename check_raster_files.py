import os, sys, re, math, json
import numpy as np
from copy import copy
try:
    from osgeo import gdal, ogr, osr
except:
    import gdal, ogr, osr
from .temp_files import *
from .logger import *

# Does not support the following gdal_translate options: -r -tr -expand -outsize -srcwin -projwin -projwin_srs -a_srs -a_coord_epoch -a_offset -a_ullr -colorinterp_X -gcp -sds -norat -noxmp
# Does not support the following gdal_rasterize options: -of -a_srs -te -tap -ts

__default_parameters__ = {
        'RasterSize_error': 0.01,
        'PixelSize_error': 0.01,
        'NoDataValue_error': 0.01,
        'DataMinimum_error': 0.01,
        'DataMaximum_error': 0.01,
        '__preserve_original_pixel_size': True,
    }

__gdal_data_type__ = ['Unknown', 'Byte', 'UInt16', 'Int16', 'UInt32', 'Int32', 'Float32', 'Float64', 'CInt32', 'CFloat32', 'CFloat64']
__int_parameters__ = ['RasterCount', 'RasterCountMin', 'RasterCountMax', 'RasterSize', 'RasterSizeMin', 'RasterSizeMax', 'DataType', 'EPSG']
__float_parameters__ = ['PixelSize', 'PixelSizeMin', 'PixelSizeMax', 'NoDataValue', 'NoDataValueMin', 'NoDataValueMax', 'DataMinimum', 'DataMinimumMin', 'DataMinimumMax', 'DataMaximum', 'DataMaximumMin', 'DataMaximumMax']
__string_parameters__ = ['Projection', 'Compression']

class GDALoptions:

    def __init__(self, gdalf):
        if gdalf == 'gdalwarp':
            self.single_options = ['tps', 'rpc', 'geoloc', 'tap', 'multi', 'q', 'overwrite', 'crop_to_cutline', 'nomd', 'setci', 'srcalpha', 'nosrcalpha', 'dstalpha']
            self.str_options = ['ovr', 'if', 'of', 'cutline', 'cl', 'cwhere', 'csql', 'cvmd']
            self.int_options = ['order', 'refine_gcps', 'wm', 'cblend']
            self.float_options = ['et', 'srcnodata', 'dstnodata']
        elif gdalf == 'gdal_translate':
            self.single_options = ['strict', 'unscale', 'epo', 'eco', 'nogcp', 'sds', 'q' 'stats', 'norat', 'noxmp']
            self.str_options = ['if', 'of', 'expand', 'colorinterp', 'scale']
            self.int_options = ['mask']
            self.float_options = ['exponent', 'a_scale', 'a_offset', 'a_nodata']
        elif gdalf == 'gdal_rasterize':
            self.single_options = ['i', 'at', '3d', 'add', 'q']
            self.str_options = ['a', 'l', 'sql', 'dialect', 'init', 'optim']
            self.int_options = ['burn']
            self.float_options = []
        else:
            raise Exception(f'Unknown GDAL function: {gdalf}')

# See details on https://gdal.org/drivers/raster/gtiff.html#raster-gtiff
__gdal_creation_options__ = ['NUM_THREADS', 'GEOREF_SOURCES', 'SPARSE_OK', 'TWF', 'RPB', 'RPCTXT', 'INTERLEAVE', 'TILED', 'BLOCKXSIZE', 'BLOCKYSIZE', 'NBITS', 'COMPRESS', 'PREDICTOR', 'DISCARD_LSB', 'JPEG_QUALITY', 'JPEGTABLESMODE', 'ZLEVEL', 'ZSTD_LEVEL', 'MAX_Z_ERROR', 'WEBP_LEVEL', 'WEBP_LOSSLESS', 'JXL_LOSSLESS', 'JXL_EFFORT', 'JXL_DISTANCE', 'PHOTOMETRIC', 'ALPHA', 'PROFILE', 'BIGTIFF', 'PIXELTYPE', 'COPY_SRC_OVERVIEWS', 'GEOTIFF_KEYS_FLAVOR', 'GEOTIFF_VERSION']
# See details on https://gdal.org/api/gdalwarp_cpp.html#_CPPv4N15GDALWarpOptions16papszWarpOptionsE
__gdal_warp_options__ = ['INIT_DEST', 'WRITE_FLUSH', 'SKIP_NOSOURCE', 'UNIFIED_SRC_NODATA', 'CUTLINE',
    'CUTLINE_BLEND_DIST', 'CUTLINE_ALL_TOUCHED', 'OPTIMIZE_SIZE', 'NUM_THREADS', 'STREAMABLE_OUTPUT',
    'SRC_COORD_PRECISION', 'SRC_ALPHA_MAX', 'DST_ALPHA_MAX', 'SAMPLE_GRID', 'SAMPLE_STEPS', 'SOURCE_EXTRA',
    'APPLY_VERTICAL_SHIFT', 'MULT_FACTOR_VERTICAL_SHIFT']
# See details on https://gdal.org/api/gdal_alg.html#_CPPv426GDALCreateRPCTransformerV2PK13GDALRPCInfoV2idPPc
__gdal_transform_options__ = ['RPC_HEIGHT', 'RPC_HEIGHT_SCALE', 'RPC_DEM', 'RPC_DEM_INTERPOLATION', 'RPC_DEM_MISSING_VALUE',
    'RPC_DEM_SRS', 'RPC_DEM_APPLY_VDATUM_SHIFT', 'RPC_PIXEL_ERROR_THRESHOLD', 'RPC_MAX_ITERATIONS', 'RPC_FOOTPRINT']
# See details on https://gdal.org/drivers/raster/gtiff.html
__gdal_open_options_geotiff__ = ['NUM_THREADS', 'GEOREF_SOURCES', 'SPARSE_OK']

def AsList(obj):
    new_obj = copy(obj)
    if isinstance(new_obj, (tuple, list)):
        return new_obj
    else:
        return [new_obj]

def CheckValue(variable, value_exact, value_presicion, value_min, value_max):
    if value_exact is not None:
        if hasattr(variable, '__float__') and (value_presicion is not None):
            try:
                if abs((variable - value_exact) / (variable + value_exact)) > value_presicion:
                    return 1
            except ZeroDivisionError:
                if abs(variable) > value_presicion:
                    return 1
        elif variable != value_exact:
            return 1
    if hasattr(variable, '__float__'):
        if value_min is not None:
            if variable < value_min:
                return 2
        if value_max is not None:
            if variable > value_max:
                return 3
    return 0

def SpatialReferenceFromEPSG(epsg):
    if epsg is not None:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg))
        return srs

def SpatialReferenceFromProj4(proj4):
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj4)
    return srs

def GetUTMzoneEPSG(lat, lon):
    return int(f'32{(6,7)[lat<0]}{round(((lon+180)%360)/6):02}')

def GetSrsUTM(rds):
    rsrs = rds.GetSpatialRef()
    tsrs = SpatialReferenceFromEPSG(4326)
    if rsrs.IsSame(tsrs):
        transform = None
    else:
        rsrs.SetAxisMappingStrategy(0)
        tsrs.SetAxisMappingStrategy(0)
        transform = osr.CoordinateTransformation(rsrs, tsrs)
    rgt = rds.GetGeoTransform(can_return_null = True)
    if rgt:
        wkt = 'POLYGON (({hx} {hy}, {hx} {ly}, {lx} {ly}, {lx} {hy}, {hx} {hy}))'.format(
            lx = rgt[0], hx = rgt[0] + rgt[1] * math.cos(rgt[2]) * rds.RasterXSize + rgt[5] * math.sin(rgt[4]) * rds.RasterYSize,
            ly = rgt[3], hy = rgt[3] + rgt[5] * math.cos(rgt[4]) * rds.RasterYSize + rgt[1] * math.sin(rgt[2]) * rds.RasterXSize,)
    else:
        wkt = f'MULTIPOINT ({",".join(["%f %f," % (point.GCPX, point.GCPY) for point in rds.GetGCPs()])})'
    try:
        area_geometry = ogr.Geometry(wkt = wkt)
        if transform:
            area_geometry.Transform(transform)
        point_geometry = area_geometry.Centroid()
    except Exception as e:
        print(e)
        return None
    return GetUTMzoneEPSG(point_geometry.GetY(), point_geometry.GetX())

def CheckBandsMatch(rpath, bands):
    if bands is None:
        return True
    rds = gdal.Open(rpath)
    if rds:
        return list(range(1, rds.RasterCount+1)) == bands

def unbanded(option):
    if re.search('^.+_\d+$', option):
        return option.split('_')[0]
    else:
        return option

def GetPercentile(histogram, min = 0.02, max = 0.98):
    total = sum(histogram)
    min_value = total * min
    max_value = total * max
    min_position = 0
    min_sum = histogram[0]
    while min_sum < min_value:
        min_position += 1
        min_sum += histogram[min_position]
    max_position = len(histogram) - 1
    max_sum = total - histogram[max_position]
    while max_sum > max_value:
        max_position -= 1
        max_sum -= histogram[max_position]
    return min_position, max_position

def DataTypeMinMaxRange(data_type):
    try:
        return [(0, 255, 256), (0,65535, 65536), (-32767, 32768, 65536), (0,4294967295,4294967296),
                (-2147483647,2147483648,4294967296), (-2147483647,2147483648,4294967296),
                (-9223372036854775807,9223372036854775808,18446744073709551616),
                (-2147483647,2147483648,4294967296), (-2147483647,2147483648,4294967296),
                (-9223372036854775807,9223372036854775808,18446744073709551616)][data_type-1]
    except:
        return (-32767, 32768, 65536)

class RasterParameters(dict):

    def __init__(self, *args, **kwargs):
        self.update(__default_parameters__)
        for arg in args:
            self.update(arg)
        self.update(**kwargs)

    def gget(self, *keys, default = None, return_key = False):
        for key in keys:
            if key in self:
                if return_key:
                    return key, self[key]
                else:
                    return self[key]
        if return_key:
            return None, default
        else:
            return default

    def valueList(self, *keys, default = None):
        return [self.get(key, default) for key in keys]

    def select(self, *args):
        return RasterParameters(**{key: self[key] for key in args})

    def alter(self, **kwargs):
        new_dict = copy(self)
        new_dict.update(**kwargs)
        return new_dict

    # Not ready yet
    def UpdateFromRasterFile(self, rpath):
        rds = gdal.Open(rpath)
        if rds is not None:
            self.update(RasterXSize = rds.RasterXSize, RasterYSize = rds.RasterYSize, RasterCount = rds.RasterCount,
                DataType = rds.GetRasterBand(1).DataType, SpatialReference = rds.GetSpatialRef(),
                GeoTransform = rds.GetGeoTransform(), NoDataValue = rds.GetRasterBand(1).GetNoDataValue())

    def CheckRaster(self, rpath):
        miss = Mismatch()
        rds = gdal.Open(rpath)
        if rds:
            miss.Check('RasterCount', rds.RasterCount, self.get('RasterCount'), None, self.get('RasterCountMin'), self.get('RasterCountMax'))
            miss.Check('RasterXSize', rds.RasterXSize, self.gget('RasterXSize', 'RasterSize'), None, self.gget('RasterXSizeMin', 'RasterSizeMin'), self.gget('RasterXSizeMax', 'RasterSizeMax'))
            miss.Check('RasterYSize', rds.RasterYSize, self.gget('RasterYSize', 'RasterSize'), None, self.gget('RasterYSizeMin', 'RasterSizeMin'), self.gget('RasterYSizeMax', 'RasterSizeMax'))
            rtrans = rds.GetGeoTransform(1)
            miss.Check('PixelXSize', abs(rtrans[1]), self.gget('PixelXSize', 'PixelSize'), None, self.gget('PixelXSizeMin', 'PixelSizeMin'), self.gget('PixelXSizeMax', 'PixelSizeMax'))
            miss.Check('PixelYSize', abs(rtrans[5]), self.gget('PixelYSize', 'PixelSize'), None, self.gget('PixelYSizeMin', 'PixelSizeMin'), self.gget('PixelYSizeMax', 'PixelSizeMax'))
            miss.Check('COMPRESS', rds.GetMetadata(domain='IMAGE_STRUCTURE').get('COMPRESSION', 'NONE'), self.get('COMPRESS'), None, None, None)
            miss.Check('DataType', [rds.GetRasterBand(i + 1).DataType for i in range(rds.RasterCount)], self.get('DataType'), None, None, None)
            miss.Check('NoDataValue', [rds.GetRasterBand(i + 1).GetNoDataValue() for i in range(rds.RasterCount)], self.get('NoDataValue'), self.get('NoDataValue_error'), self.get('NoDataValueMin'), self.get('NoDataValueMax'))
            check_values = self.valueList('DataMinimum', 'DataMinimumMin', 'DataMinimumMax', 'DataMaximum', 'DataMaximumMin', 'DataMaximumMax')
            for i in range(rds.RasterCount):
                if all([check_value is None for check_value in check_values]):
                    break
                rband = rds.GetRasterBand(i + 1)
                rarr_values = np.unique(rband.ReadAsArray())
                rarr_values = rarr_values[rarr_values != rband.GetNoDataValue()]
                try:
                    check_min = miss.Check('DataMinimum', np.amin(rarr_values), check_values[0], self.get('DataMinimum_error'), check_values[1], check_values[2])
                    if check_min:
                        for i in range(1, 4):
                            if i in check_min:
                                check_values[i-1] = None
                except AssertionError:
                    pass
                try:
                    check_max = miss.Check('DataMaximum', np.amax(rarr_values), check_values[3], self.get('DataMaximum_error'), check_values[4], check_values[5])
                    if check_max:
                        for i in range(1, 4):
                            if i in check_max:
                                check_values[i+2] = None
                except AssertionError:
                    pass
            tsrs = self.GetSpatialReference()
            if tsrs:
                rsrs = rds.GetSpatialRef()
                if tsrs == 'UTM':
                    if re.search('^32[67]\d\d$', rsrs.GetAttrValue('AUTHORITY',1)):
                        pass
                    else:
                        miss['Proj4'] = (f"{rsrs.ExportToProj4()}", f"{self.GetSpatialReference(rpath=rpath).ExportToProj4()}")
                elif not rsrs.IsSame(tsrs):
                    miss['Proj4'] = (f"{rsrs.ExportToProj4()}", f"{tsrs.ExportToProj4()}")
        return miss

    def GetSpatialReference(self, rpath = None, for_extent = False):
        dkey = ('', 'Extent')[for_extent]
        key, srs = self.gget(*tuple([dkey+key for key in ['SpatialReference', 'Proj4', 'ProjWkt', 'EPSG']]), return_key=True)
        if (key is not None) and (srs is not None):
            if srs == 'UTM':
                if rpath is not None:
                    return SpatialReferenceFromEPSG(GetSrsUTM(gdal.Open(rpath)))
                else:
                    return srs
            elif key == 'SpatialReference':
                return srs
            elif key == 'Proj4':
                return SpatialReferenceFromProj4(srs)
            elif key == 'ProjWkt':
                return osr.SpatialReference(wkt = srs)
            elif key == 'EPSG':
                return SpatialReferenceFromEPSG(srs)
            else:
                raise Exception('Unknown srs format:', srs)

    def GetGDALoptions(self, gdalf):
        options = []
        gdalo = GDALoptions(gdalf)
        for option in self:
            value = self.get(option)
            if value is not None:
                if unbanded(option) in gdalo.single_options:
                    options.append(f' -{option}')
                elif unbanded(option) in gdalo.str_options:
                    options.append(f' -{option} {value}')
                elif unbanded(option) in (gdalo.int_options + gdalo.float_options):
                    options.append(f' -{option} {value}')
        return options

    def GetGDALoptionsList(self, type):
        options = []
        gdal_options = {'creation': __gdal_creation_options__, 'transform': __gdal_transform_options__, 'warp': __gdal_warp_options__, 'open': __gdal_open_options_geotiff__}.get(type, [])
        for option in self:
            option_check = option.upper()
            if option_check in gdal_options:
                options.append(f'{option_check}={str(self[option])}')
        return options

    # Possible gdalf values: gdalwarp, gdal_translate
    def GDALoptionString(self, gdalf):
        options = self.GetGDALoptions(gdalf)
        # print(gdalf, options)
        if self.get('DataType') is not None:
            options.append(f' -ot {__gdal_data_type__[self.get("DataType", 0)]}')
        if self.get('NoDataValue') is not None:
            nodata_apx = {"gdalwarp": "dst", "gdal_translate": "a_", "gdal_rasterize": "a_"}.get(gdalf, "")
            options.append(f' -{nodata_apx}nodata {self["NoDataValue"]}')
        if gdalf == 'gdalwarp':
            srs = self.GetSpatialReference()
            if self.get('Method') is not None:
                options.append(f' -r {self["Method"]}')
            if srs is not None:
                if isinstance(srs, osr.SpatialReference):
                    options.append(f' -t_srs "{srs.ExportToProj4()}"')
            if self.get('Extent') is not None:
                options.append(f' -te {" ".join(self["Extent"])}')
                extent_srs = self.GetSpatialReference(1)
                if extent_srs is not None:
                    if isinstance(extent_srs, osr.SpatialReference):
                        options.append(f' -te_srs "{extent_srs.ExportToProj4()}"')
            rasterXsize = self.gget('RasterXSize', 'RasterSize')
            rasterYsize = self.gget('RasterYSize', 'RasterSize')
            if rasterXsize and rasterYsize:
                options.append(f' -ts {rasterXsize} {rasterYsize}')
            pixelXsize = self.gget('PixelXSize', 'PixelSize')
            pixelYsize = self.gget('PixelYSize', 'PixelSize')
            if pixelXsize and pixelYsize:
                options.append(f' -tr {pixelXsize} {pixelYsize}')
            if self.get('SourceNoDataValue') is not None:
                options.append(f' -srcnodata {self["SourceNoDataValue"]}')
        elif gdalf == 'gdal_translate':
            if self.get('Bands') is not None:
                options.append(''.join([f' -b {int(bandnum)}' for bandnum in self.get('Bands')]))
            if self.get('mask') is not None:
                options.append(' --config GDAL_TIFF_INTERNAL_MASK YES')
        elif gdalf == 'gdal_rasterize':
            if sum([self.get(key) is not None for key in ['3d', 'burn', 'a']]) != 1:
                raise Exception('One and only one of "3d", "burn" or "a" is required')
        option_type = {'gdalwarp': ['creation', 'transform', 'warp', 'open'],
                       'gdal_translate': ['creation', 'open'],
                       'gdal_rasterize': ['creation', 'open'],
                       }.get(gdalf, [])
        for type in option_type:
            options += [f' -{type[0].lower()}o {option}' for option in self.GetGDALoptionsList(type)]
        return ''.join(options)

    def GDALprocessRaster(self, rpath, tpath, gdalf):
        return os.system(f'{gdalf} {self.GDALoptionString(gdalf)} {rpath} {tpath}')

    def GDALTranslateRaster(self, rpath, tpath):
        ds1 = gdal.Open(rpath)
        option_string = self.GDALoptionString('gdal_translate')
        # print(option_string)

        # if option GDAL_TIFF_INTERNAL_MASK set to YES
        if self.get('mask') is not None:
            option_string = option_string.replace(' --config GDAL_TIFF_INTERNAL_MASK YES', '')
            gdal.SetConfigOption('GDAL_TIFF_INTERNAL_MASK', 'YES')

        translate_options = gdal.TranslateOptions(gdal.ParseCommandLine(option_string))
        gdal.Translate(tpath, ds1, options=translate_options)
        ds1 = None

    def GDALWarpRaster(self, rpath, tpath):
        option_string = self.GDALoptionString('gdalwarp')
        # print(option_string)
        warp_options = gdal.WarpOptions(gdal.ParseCommandLine(option_string))
        gdal.Warp(tpath, rpath, options=warp_options)

    def GDALRasterize(self, vpath, tpath):
        option_string = self.GDALoptionString('gdal_rasterize')
        # print(option_string)
        rasterize_options = gdal.RasterizeOptions(gdal.ParseCommandLine(option_string))
        gdal.Rasterize(gdal.Open(tpath, 1), vpath, options=rasterize_options)

class Mismatch(dict):

    def __str__(self):
        return '\n'.join([f'{key}: {self[key]}' for key in self])

    def __repr__(self, header = ''):
        result = header
        if len(self) > 0:
            for key in self:
                try:
                    result += '\n\t{}: {} -> {}'.format(key, self[key][0], self[key][1])
                except:
                    result += '\n\tERROR: Incorrect description for {}'.format(key)
        else:
            result += ' is OK'
        return result

    def Report(self, header = ''):
        print(self.__repr__(header = header))

    def Check(self, key, variable, value_exact, value_presicion, value_min, value_max):
        assert key not in self
        variable_list = AsList(variable)
        final = []
        checkNone = False

        if isinstance(value_exact, str):
            try:
                value_exact = float(value_exact)
            except ValueError:
                if value_exact.upper() == 'NONE':
                    checkNone = True
                    value_exact = None

        for var in variable_list:
            if (checkNone and var is None) or ((var is None) and (value_exact is not None)):
                final.append(1)
                self[key] = (var, value_exact)
            elif var is not None:
                result = CheckValue(var, value_exact, value_presicion, value_min, value_max)
                if result == 1:
                    self[key] = (var, value_exact)
                    value_exact = None
                elif result == 2:
                    self[key + 'Min'] = (var, value_min)
                    value_min = None
                elif result == 3:
                    self[key + 'Max'] = (var, value_max)
                    value_max = None
                else:
                    continue
                final.append(result)
        final.sort()
        return final

    def ErrorClass(self):
        unable_list = ['RasterCount', 'RasterCountMin', 'RasterCountMax', 'RasterSize', 'RasterSizeMin', 'RasterSizeMax', 'RasterXSize', 'RasterXSizeMin', 'RasterXSizeMax', 'RasterYSize', 'RasterYSizeMin', 'RasterYSizeMax']
        reproject_list = ['SpatialReference', 'ProjWkt', 'Proj4', 'GeoTransform', 'PixelSize', 'PixelSizeMin', 'PixelSizeMax', 'PixelXSize', 'PixelXSizeMin', 'PixelXSizeMax', 'PixelYSize', 'PixelYSizeMin', 'PixelYSizeMax']
        # translate_list = ['NoDataValue', 'NoDataValueMin', 'NoDataValueMax', 'DataType', 'Compression']
        cutdatapeaks_list = ['DataMinimum', 'DataMinimumMin', 'DataMinimumMax', 'DataMaximum', 'DataMaximumMin', 'DataMaximumMax']
        error_list = [0, 0]
        for key in self:
            if key in unable_list:
                return 0
            elif key in reproject_list:
                error_list[0] = 1
            # elif (key in translate_list) and error_list[0]==0:
                # error_list[0] = 1
            elif key in cutdatapeaks_list:
                error_list[1] = 1
        return error_list

#@StopFromStorage is not used by now, as it is necessary for only DataMin/Max check; all other checks may be made for files on storage directly
def CheckRasterParameters(rpath, **parameters):
    return RasterParameters(**parameters).CheckRaster(rpath)

@StopFromStorage
def ReprojectRaster(rpath, tpath, **par):
    par = RasterParameters(**par)
    rds = gdal.Open(rpath)
    if rds:
        if par.get('__preserve_original_pixel_size') and par.gget('PixelSize', 'PixelXSize', 'PixelYSize') is None:
            srs = rds.GetSpatialRef()
            if srs.IsProjected():
                trans = rds.GetGeoTransform()
                srs_unit = srs.GetLinearUnits()
                par.update(PixelXSize = abs(trans[1]*srs_unit), PixelYSize = abs(trans[5]*srs_unit))
        if par.GetSpatialReference() == 'UTM':
            par.update(SpatialReference = par.GetSpatialReference(rpath = rpath))
    rds = None
    # return par.GDALprocessRaster(rpath, tpath, 'gdalwarp')
    return par.GDALWarpRaster(rpath, tpath)

@StopFromStorage
def RasterizeVector(vpath, tpath, **par):
    par = RasterParameters(**par)
    return par.GDALRasterize(vpath, tpath)

@StopFromStorage
def Rasterize(rpath, vpath, tpath, **mask_parameters):
    mask_parameters = RasterParameters(**mask_parameters)
    SetRaster(rpath, tpath, None,
              Bands=[1], NoDataValue=0, COMPRESS='DEFLATE', NBITS=16)
    tds = gdal.Open(tpath, 1)
    tds.GetRasterBand(1).WriteArray(np.zeros((tds.RasterYSize, tds.RasterXSize), bool))
    tds = None
    if vpath is not None:
        return RasterizeVector(vpath, tpath, **mask_parameters)

@StopFromStorage
def CopyRaster(rpath, tpath, **parameters):
    # return RasterParameters(**parameters).GDALprocessRaster(rpath, tpath, 'gdal_translate')
    return RasterParameters(**parameters).GDALTranslateRaster(rpath, tpath)

@StopFromStorage
def CutDataPeaks(rpath, minimum = None, maximum = None):
    if (minimum is not None) or (maximum is not None):
        rds = gdal.Open(rpath, 1)
        if rds:
            for band_num in range(1, rds.RasterCount + 1):
                rband = rds.GetRasterBand(band_num)
                rarr = rband.ReadAsArray()
                if minimum is not None:
                    rarr[rarr < minimum] = minimum
                if maximum is not None:
                    rarr[rarr > maximum] = maximum
                rband.WriteArray(rarr)
            rds = None

@StopFromStorage
def SetRaster(rpath, tpath, miss, **parameters):
    assert rpath != tpath
    if os.path.exists(tpath) and (not parameters.get('overwrite')):
        raise Exception('File already exists: ' + tpath)
    if not isinstance(miss, Mismatch):
        miss = CheckRasterParameters(rpath, **parameters)
    error_class = miss.ErrorClass()
    # miss.Report(str(error_class))
    if error_class[0]:
        if CheckBandsMatch(rpath, parameters.get('Bands')) and \
                (not any([(key in parameters) for key in ['mask', 'scale', 'exponent', 'colorinterp']])):
            ReprojectRaster(rpath, tpath, **parameters)
        else:
            temp_path = TempName(SplitPath(rpath)[1], 'tif')
            ReprojectRaster(rpath, temp_path, **RasterParameters(**parameters).alter(DataType = None, COMPRESS = None))
            CopyRaster(temp_path, tpath, **parameters)
            DeleteFile(temp_path)
    else:
        CopyRaster(rpath, tpath, **parameters)
    if error_class[1]:
        par = RasterParameters(**parameters)
        CutDataPeaks(tpath, minimum=par.gget('DataMinimum', 'DataMimimumMin'), maximum=par.gget('DataMaximum', 'DataMaximumMax'))

@StopFromStorage
def SetRGB(rpath, tpath, method, min_max_limits, **parameters):
    par = RasterParameters(**parameters)
    rds = gdal.Open(rpath)
    bands = par.get('Bands', [1,2,3])
    assert len(bands) == 3
    if not 'Bands' in par:
        par['Bands'] = bands

    if not method in range(4):
        raise Exception(f'Incorrect method {method}')
    elif method == 0:
        min_max_limits = [(1,255),(1,255),(1,255)]
    elif len(min_max_limits) == 2:
        min_max_limits = [min_max_limits, min_max_limits, min_max_limits]
    elif len(min_max_limits) != 3:
        raise Exception(f'Incorrect min_max {min_max_limits}')

    for bandnum, min_max in zip(bands, min_max_limits):
        minimum, maximum = min_max
        total_data_limits = par.get('TotalDataLimits')
        if method == 1:
            mean, std = rds.GetRasterBand(bandnum).ComputeStatistics(0)[2:4]
            minimum = mean - minimum * std
            maximum = mean + maximum * std
        elif method == 2:
            data_min, data_max, data_range = DataTypeMinMaxRange(rds.GetRasterBand(bandnum).DataType)
            histogram = rds.GetRasterBand(bandnum).GetHistogram(data_min, data_max, data_range, approx_ok = False)
            minimum, maximum = GetPercentile(histogram, minimum, maximum)
            minimum -= data_min
            maximum -= data_min
        if total_data_limits:
            minimum = max(minimum, total_data_limits[0])
            maximum = min(maximum, total_data_limits[1])
        # print(minimum, maximum)
        par[f'scale_{bandnum}'] = f'{minimum} {maximum} 1 255'
    rds = None
    return SetRaster(rpath, tpath, None, **par.alter(DataType = 1, scale = None, mask = 4, NoDataValue = 'NONE'))

# Used in previous version
class ImageChecker:

    def __init__(self, parameters = {}):
        self.parameters = RasterParameters(parameters)
        self.results = {}

    def SetCheckParameters(self, **parameters):
        # for key in parameters:
            # if not key in self.parameters:
                # print('Unknown key:', key)
        self.parameters.update(parameters)

    def GetCheckParameters(self, **parameters):
        p = copy(self.parameters)
        p.update(parameters)
        return p

    def CheckImage(self, rpath, **parameters):
        miss = CheckImageParameters(rpath, **self.GetCheckParameters(**parameters))
        # miss.Report(rpath)
        self.results[rpath] = miss

    def Summary(self, print_summary = False):
        results_values = {}
        for val in self.results.values():
            key = ';'.join([f'{err}: {val[err][0]} -> {val[err][1]}' for err in val])
            if key in results_values:
                results_values[key] += 1
            else:
                results_values[key] = 1
        if '' in results_values:
            results_values['OK'] = results_values.pop('')
        if print_summary:
            scroll(results_values, header = 'Summary:')
        return results_values

    def ResultsToJson(self, json_path, print_summary = False):
        with open(json_path, 'w') as write_file:
            json.dump({'summary': self.Summary(print_summary),
                        'results': self.results},
                      write_file, indent=4)

    def CheckImagesInFolder(self, folder_in, target_path = None, miss_path = None, **parameters):
        for rpath in Files(folder_in, 'tif', target_path =target_path, miss_path=miss_path):
            self.CheckImage(rpath, **parameters)

    def FixImage(self, rpath, tpath, **parameters):
        p = self.GetCheckParameters()
        miss = self.results.get(rpath, {})
        if 'SpatialReference' in miss:
            p['reproject'] = True
        if any([key in miss for key in ('DataMinimumMin', 'DataMaximumMax')]):
            p['cut_data_peaks'] = True
        p.update(**parameters)
        # FixImage(rpath, tpath, **p)
