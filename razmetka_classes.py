from .check_raster_files import *
from .vector_geometry import clipRasterWithMask
from .logger import *
from .paths import openFile
from osgeo import ogr
import numpy as np
import pandas as pd
import json


__default_mask_parameters__ = {
        'empty': False,         # If True and vector field is not found an empty mask is created
        'burn': None,           # Burn mask value if the vector attribute field is not found
        'a': 'gridcode',        # Vector field name containing the attribute value
        'crop_mask': True,      # If True the mask raster is clipped with data != NoDataValue (0)
        'crop_data': False,     # If True the data raster is clipped with mask != NoDataValue (0)
        'replace': {},          # Replace mask values
        # 'delete_values'       # Abort all masks with these values -- not needed yet
        'set_satid': 'images',  # Type of satellite images on the mask
        'set_objid': 'full',    # Type of objects on the mask
        'set_appendix': '',     # Appendix to mask filenames
        'set_unmarked': None,   # Set value for all unmarked data pixels
        'overwrite': False,
        '__errors__': [],
        'Method': 'average',    # Reprojection method: set to average here
    }


@StopFromStorage
def RasterStatistics(rpath):
    rds = gdal.Open(rpath)
    report = {}
    if rds is not None:
        report['RasterCount'] = rds.RasterCount
        report['RasterXSize'] = rds.RasterXSize
        report['RasterYSize'] = rds.RasterYSize
        report['DataType'] = rds.GetRasterBand(1).DataType
        srs = rds.GetSpatialRef()
        if srs:
            report['EPSG'] = srs.GetAttrValue('AUTHORITY',1)
            report['Proj4'] = srs.ExportToProj4()
        trans = rds.GetGeoTransform(can_return_null=True)
        if trans:
            report['PixelXSize'] = abs(trans[1])
            report['PixelYSize'] = abs(trans[5])
        report['Bands'] = {}
        for bnum in range(1, rds.RasterCount+1):
            min, max, mean, std = rds.GetRasterBand(bnum).ComputeStatistics(0)
            report['Bands'][bnum] = {'Min': float(min), 'Max': float(max),
                                     'Mean': float(mean), 'StD': float(std)}
    return report


@StopFromStorage
def ReplaceMaskValues(mpath, replace):
    mds = gdal.Open(mpath, 1)
    mbnd = mds.GetRasterBand(1)
    marr = mbnd.ReadAsArray()
    report = {}
    values, counts = np.unique(marr, return_counts=True)
    for val, count in zip(values, counts):
        if val in replace:
            marr[marr==val] = replace[val]
            val = replace[val]
        if val in report:
            report[int(val)] += int(count)
        else:
            report[int(val)] = int(count)
    mbnd.WriteArray(marr)
    mds = None
    return report


@StopFromStorage
def setUnmarkedDataMaskValues(rpath, tpath, set_value):
    rds = gdal.Open(rpath)
    band = rds.GetRasterBand(1)
    mask_band = band.GetMaskBand()
    mask_arr = mask_band.ReadAsArray().astype(bool)
    tds = gdal.Open(tpath, 1)
    tband = tds.GetRasterBand(1)
    tarr = tband.ReadAsArray()
    tarr[mask_arr * (tarr == 0)] = set_value
    tband.WriteArray(tarr)
    tds = None


#@StopFromStorage
def SetVector(vpath, tpath, replace, attr_field, burn, overwrite=False):
    if vpath is None:
        return {}
    elif tpath is None:
        tds = ogr.Open(vpath)
    else:
        if tpath != vpath:
            CopySHP(vpath, tpath, overwrite=overwrite)
        tds = ogr.Open(tpath, 1)
    if tds is None:
        return {}
    else:
        tlr = tds.GetLayer()
        if burn is not None:
            return {int(replace.get(burn, burn)): len(tlr)}
        else:
            report = {}
            for feat in tlr:
                val = feat.GetField(attr_field)
                if val in replace:
                    if tpath is not None:
                        feat.SetField(attr_field, replace[val])
                        tlr.SetFeature(feat)
                    val = replace[val]
                if val in report:
                    report[int(val)] += 1
                else:
                    report[int(val)] = 1
            tds = None
            return report


@StopFromStorage
def SaveJson(json_path, dictionary):
    with open(json_path, 'w') as write_file:
        json.dump(dictionary, write_file, indent=4, ensure_ascii=False)


class MaskParameters(RasterParameters):

    def __init__(self, *args, **kwargs):
        self.update(__default_mask_parameters__)
        self.update(RasterParameters(*args, **kwargs))

    def select(self, *args):
        return MaskParameters(**{key: self[key] for key in args})

    def isValid(self, rpath, vpath, trpath, tmpath, tvpath,
                check_target=False, check_folder=False):
        rpath_check = CheckPathValidity(rpath, True, True, False, False)
        vpath_check = CheckPathValidity(vpath, not self['empty'], True,
                                        False, False)
        if check_target:
            trpath_check = CheckPathValidity(trpath, True, False,
                                        not self['overwrite'], check_folder)
            tmpath_check = CheckPathValidity(tmpath, True, False,
                                        not self['overwrite'], check_folder)
            tvpath_check = CheckPathValidity(tvpath, False, False,
                                        False, check_folder)
        else:
            trpath_check = tmpath_check = tvpath_check = None
        result = True
        for check in [rpath_check, vpath_check, trpath_check, tmpath_check,
                      tvpath_check]:
            if check is not None:
                self['__errors__'].append(check)
                result = False
        return result

    @log_me
    def analyzeData(self, rpath, vpath, trpath, tmpath, tvpath,
                    calculate_raster_stats=False):
        if self.isValid(rpath, vpath, trpath, tmpath, tvpath,
                        check_target=False, check_folder=False):
            source_data_report = {'mask_feature_count':
                        SetVector(vpath, None, {}, self['a'], self['burn'])}
            if calculate_raster_stats:
                source_data_report.update(RasterStatistics(rpath))
            return {
    'SourceData': source_data_report,
    'ResultData': {'mask_feature_count':
            SetVector(vpath, None, self['replace'], self['a'], self['burn'])},
                    }
        else:
            return {}

    def SetRaster(self, rpath, tpath):
        SetRaster(rpath, tpath, None, **self)
        return RasterStatistics(tpath)

    def SetMask(self, rpath, vpath, tpath):
        Rasterize(rpath, vpath, tpath, **self.select('a', 'burn'))
        if self.get('crop_data'):
            clipRasterWithMask(rpath, tpath)
        if self.get('crop_mask'):
            clipRasterWithMask(tpath, rpath)
        if self.get('set_unmarked') is not None:
            setUnmarkedDataMaskValues(rpath, tpath, self.get('set_unmarked'))
        # CropRasterDataMask(rpath, tpath, self.get('crop_data'), self.get('crop_mask'))
        return ReplaceMaskValues(tpath, self.get('replace'))

    def SetVector(self, vpath, tpath):
        return SetVector(vpath, tpath, replace = self['replace'], attr_field = self['a'],
                         burn = self['burn'], overwrite = self['overwrite'])

    @log_me
    @StopFromStorage
    def writeData(self, rpath, vpath, trpath, tmpath, tvpath):
        if self.isValid(rpath, vpath, trpath, tmpath, tvpath,
                        check_target=True, check_folder=True):
            report = self.SetRaster(rpath, trpath)
            report['mask_feature_count'] = self.SetVector(vpath, tvpath)
            report['mask_pixel_count'] = self.SetMask(trpath, tvpath, tmpath)
            return {'ResultData': report, 'MaskParameters': dict(self)}
        else:
            return {}


class Mask:

    def __init__(self, id, parameters, rpath=None, vpath=None):
        self.id = id
        self.par = MaskParameters(**parameters)
        self.rpath = rpath
        self.vpath = vpath
        self.trpath = None
        self.tmpath = None
        self.tvpath = None
        self.report = {}
        self.errors = []
        self.analyzed = False

    def __str__(self):
        return f'''Mask {self.id} {('Invalid', 
                'Valid')[self.par.IsValid(*self.paths())]}'''

    def paths(self):
        return self.rpath, self.vpath, self.trpath, self.tmpath, self.tvpath

    def GetTargetPaths(self, set_folder, create_subdirs=False):
        if self.rpath is None:
            print(self.id, 'Raster path not found, cannot set target paths')
            return 1
        trfolder = os.path.join(set_folder, 'images',
                                self.par.get('set_satid', ''))
        self.trpath = FullPath(trfolder, self.id +
                               self.par.get('set_appendix', ''), 'tif')
        tmfolder = os.path.join(set_folder, 'masks',
                                self.par.get('set_objid', ''),
                                self.par.get('set_satid', ''))
        self.tmpath = FullPath(tmfolder, self.id +
                               self.par.get('set_appendix', ''), 'tif')
        if create_subdirs:
            SureDir(trfolder, tmfolder)
        if self.vpath is not None:
            tvfolder = os.path.join(set_folder, 'vector')
            self.tvpath = FullPath(tvfolder, self.id +
                                   self.par.get('set_appendix', ''), 'shp')
            if create_subdirs:
                SureDir(tvfolder)

    def analyzeMask(self):
        self.report.update(self.par.analyzeData(*self.paths()),
                           calculate_raster_stats = not self.analyzed)
        self.analyzed = True
        return self.report

    def writeMask(self):
        self.report.update(self.par.writeData(*self.paths()))
        return self.report

    def reportList(self):
        result_data = self.report.get('ResultData', {})
        bands_data = result_data.get('Bands', {})
        min_vals = max_vals = []
        for bandnum in bands_data:
            band_data = bands_data[bandnum]
            if 'Min' in band_data:
                min_vals.append(band_data['Min'])
            if 'Max' in band_data:
                max_vals.append(band_data['Max'])
        if min_vals:
            minimum = min(min_vals)
        else:
            minimum = None
        if max_vals:
            maximum = max(max_vals)
        else:
            maximum = None
        keys = list(result_data.get('mask_pixel_count', {}).keys())
        keys.sort()
        return [*self.paths(), result_data.get('RasterCount'),
                result_data.get('RasterXSize'),
                result_data.get('RasterYSize'), result_data.get('DataType'),
                result_data.get('PixelXSize'), result_data.get('PixelYSize'),
                minimum, maximum, ' '.join([str(val) for val in keys])]


class Razmetka(dict):

    def __init__(self, target_folder='', **set_parameters):
        self.tfolder = target_folder
        self.spar = MaskParameters(set_parameters)
        self.legend = Legend()
        pass

    def __repr__(self):
        actual_legend = '\n\t'.join(
            [f'{id}: {self.legend.get(id, "UNKNOWN INDEX")}'
             for id in self.MaskStats()])
        source_data = '\n\n'.join(
            [f"{key}:\n\t{self[key].rpath}\n\t{self[key].vpath}"
             for key in self])
        return f'''Set of {len(self)} masks\n
            Target folder: {self.tfolder}\n  
            Mask values: \n\t{actual_legend}\n
            Masks:\n\n{source_data}'''

    def GetRaster(self, folder, parameters, target_path=None, miss_path=None):
        for file in Files(folder, 'tif',
                          target_path=target_path, miss_path=miss_path):
            id = SplitPath(file)[1]
            if id in self:
                self[id].rpath = file
            else:
                self[id] = Mask(id, self.spar.alter(**parameters), rpath=file)

    def GetVector(self, folder, parameters, target_path=None, miss_path=None):
        for file in Files(folder, ['shp'],
                          target_path=target_path, miss_path=miss_path):
            id = SplitPath(file)[1]
            if id in self:
                self[id].vpath = file
            else:
                self[id] = Mask(id, self.spar.alter(**parameters), vpath=file)

    def GetSourceData(self, folder, parameters,
                      target_path=None, miss_path=None):
        self.GetRaster(folder, parameters,
                       target_path=target_path, miss_path=miss_path)
        self.GetVector(folder, parameters,
                       target_path=target_path, miss_path=miss_path)

    def UpdateParameters(self, **set_parameters):
        self.spar.update(**set_parameters)
        for id in self:
            self[id].par.update(**set_parameters)

    def CheckTargetFolder(self, target_folder=None):
        if target_folder is not None:
            self.tfolder = target_folder
        if not self.tfolder:
            print('Target folder is empty, cannot define target paths')
        return self.tfolder

    def analyzeSet(self, target_folder=None):
        tfolder = self.CheckTargetFolder(target_folder=target_folder)
        if tfolder:
            for id in self:
                try:
                    if not self[id].GetTargetPaths(tfolder):
                        self[id].analyzeMask()
                except Exception as e:
                    self[id].errors.append(f'Estimate: {e}')
        return { id : self[id].report for id in self }

    def writeSet(self, target_folder=None):
        tfolder = self.CheckTargetFolder(target_folder=target_folder)
        if tfolder:
            for id in self:
                try:
                    if not self[id].GetTargetPaths(tfolder, create_subdirs=True):
                        if not self[id].analyzed:
                            self[id].analyzeMask()
                        self[id].writeMask()
                except Exception as e:
                    self[id].errors.append(f'WriteMask: {e}')
        return {id: self[id].report for id in self}

    def ErrorsReport(self):
        errors_report = {}
        for id in self:
            if self[id].errors:
                errors_report[id] = self[id].errors
        return errors_report

    def MaskStats(self):
        mask_stats = {}
        for id in self:
            feature_count = self[id].report.get('ResultData',{}).get('mask_feature_count',{})
            pixel_count = self[id].report.get('ResultData', {}).get('mask_pixel_count', {})
            for count_dict, count_key in zip([feature_count, pixel_count],
                                             ['feature_count', 'pixel_count']):
                for val in count_dict:
                    if val in mask_stats:
                        if count_key in mask_stats[val]:
                            mask_stats[val][count_key] += count_dict[val]
                        else:
                            mask_stats[val][count_key] = count_dict[val]
                    else:
                        mask_stats[val] = {count_key: count_dict[val]}
        for val in mask_stats:
            mask_stats[val]['legend'] = self.legend.get(val, 'UNKNOWN INDEX')
        return mask_stats

    def makeSet(self, target_folder = None):
        tfolder = self.CheckTargetFolder(target_folder=target_folder)
        summary = {'SetFolder': tfolder, 'StartTime': datetime.now().isoformat()}
        report = self.writeSet(tfolder)
        summary.update({'MasksCount': len(report),
                        'FinishTime': datetime.now().isoformat()})
        SaveJson(FullPath(tfolder, os.path.split(tfolder)[1], 'json'), {
            'Summary': summary,
            'MaskStats': self.MaskStats(),
            'Report': report,
            'Errors': self.ErrorsReport(),
        })

    def reportResults(self, report_xls_path):
        ids = list(self.keys())
        ids.sort()
        report_df = pd.DataFrame([ self[id].reportList() for id in ids ], ids,
            ['Исходный растр', 'Исходный вектор', 'Данные', 'Маска', 'Вектор',
             'Число каналов', 'Число столбцов', 'Число строк', 'Тип данных',
             'Размер пикселя X', 'Размер пикселя Y',
             'Минимум данных', 'Максимум данных', 'Значения маски'])
        report_df.to_excel(report_xls_path)

    def writeMetadata(self, target_folder=None, summary={}):

        tfolder = self.CheckTargetFolder(target_folder=target_folder)

        if tfolder:
            set_name = os.path.split(tfolder)[1]
            mask_stats = self.MaskStats()

            SaveJson(FullPath(tfolder, set_name, 'json'), {
                'Summary': summary,
                'MaskStats': mask_stats,
                'Report': {id: self[id].report for id in self},
                'Errors': self.ErrorsReport(),
            })

            if mask_stats:
                mask_stats_path = FullPath(tfolder, 'mask_values_count.xls')
                csv_legend_path = FullPath(tfolder, 'mask_values.csv')
                mask_values = list(mask_stats.keys())
                mask_values.sort()
                mask_counts_report = pd.DataFrame(index=mask_values,
                    columns=['mask_feature_count', 'mask_pixel_count',
                             'legend']).from_dict(mask_stats, orient='index')
                mask_counts_report.to_excel(mask_stats_path)
                mask_counts_report[['legend']].to_csv(csv_legend_path, sep=';',
                                                      encoding='1251')
                # Need to fix encoding for different systems

            ids = list(self.keys())
            ids.sort()
            report_df = pd.DataFrame([self[id].reportList() for id in ids], ids,
                                     ['Исходный растр', 'Исходный вектор',
                                      'Данные', 'Маска', 'Вектор',
                                      'Число каналов', 'Число столбцов',
                                      'Число строк', 'Тип данных',
                                      'Размер пикселя X', 'Размер пикселя Y',
                                      'Минимум данных', 'Максимум данных',
                                      'Значения маски'])
            report_df.to_excel(FullPath(tfolder, set_name, 'xls'))


class LegendError(Exception):

    def __int__(self, *args, **kwargs):
        pass


def legendFromXls(xlspath, sheet=0, key_column='Значение',
                  legend_column='Описание'):
    legendDataFrame = pd.read_excel(xlspath, sheet_name=sheet)
    if key_column in legendDataFrame:
        keys = legendDataFrame[key_column].values
    else:
        raise LegendError(f'Key column not found: {key_column}')
    if legend_column in legendDataFrame:
        values = legendDataFrame[legend_column].values
    else:
        raise LegendError(f'Legend column not found: {legend_column}')
    return {key: value for key, value in zip(keys, values)}


class Legend(dict):

    def updateFromXls(self, xlspath=None, sheet=0, key_column='Значение',
                      legend_column='Описание'):

        if xlspath is None:
            xlspath = self.xlspath
        else:
            self.xlspath = xlspath

        self.update(legendFromXls(xlspath, sheet=sheet, key_column=key_column,
                                  legend_column=legend_column))
        self.sheet = sheet
        self.legend_column = legend_column

    def openXls(self):
        openFile(self.xlspath)
