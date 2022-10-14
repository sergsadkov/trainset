import os
import shutil
import re
import subprocess
import platform


def SplitPath(path):
    if os.path.isdir(path):
        return path, '', ''
    folder, file = os.path.split(path)
    name, ext = os.path.splitext(file)
    return folder, name, ext[1:]


def Name(path):
    return SplitPath(path)[1]


def FullPath(folder, file, ext = ''):
    return f'{folder}\\{file}{("",".")[bool(ext)]}{ext.lstrip(".")}'


def SureDir(*folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)


def CopyFile(path_in, path_out, overwrite = False):
    if (not os.path.exists(path_out)) or overwrite:
        if os.path.exists(path_in):
            try:
                shutil.copyfile(path_in, path_out)
            except Exception as e:
                pass


def DeleteFile(file):
    if os.path.exists(file):
        try:
            os.remove(file)
        except Exception as e:
            pass


def CopySHP(file_in, file_out, overwrite = False, ext_list = ['shp', 'dbf', 'shx', 'prj', 'sbn', 'sbx', 'cpg']):
    folder_in, name_in, ext_in = SplitPath(file_in)
    folder_out, name_out, ext_out = SplitPath(file_out)
    for ext in ext_list:
        CopyFile(FullPath(folder_in, name_in, ext), FullPath(folder_out, name_out, ext), overwrite = overwrite)


def DeleteSHP(file, ext_list = ['shp', 'dbf', 'shx', 'prj', 'sbn', 'sbx', 'cpg']):
    folder, name, ext = SplitPath(file)
    for ext in ext_list:
        DeleteFile(FullPath(folder, name, ext))


def Files(folder, extensions = None, target_path = None, miss_path = None):
    files = []
    if extensions is not None:
        if isinstance(extensions, (tuple, list)):
            extensions = list(extensions)
        else:
            extensions = [extensions]
        exts = ['.' + str(ext).lower().lstrip('.') for ext in extensions]
    for corner, _folders, _files in os.walk(folder):
        if miss_path:
            if re.search(miss_path, corner):
                continue
        for file in _files:
            if miss_path:
                if re.search(miss_path, file):
                    continue
            if extensions is not None:
                if all(not file.lower().endswith(ext) for ext in exts):
                    continue
            if target_path:
                if not re.search(target_path, file):
                    continue
            files.append(corner + '\\' + file)
    return files


def CheckPathValidity(path, forbid_none, must_exist, must_not_exist, must_folder_exist):
    if path is None:
        if forbid_none:
            return 'Path not set'
    elif os.path.exists(path):
        if must_not_exist:
            return f'File exists: {path}'
    else:
        if must_exist:
            return f'File not found: {path}'
        elif must_folder_exist and not os.path.exists(os.path.split(path)[0]):
            return f'Folder does not exist: {os.path.split(path)[0]}'


def ClearFolder(folder):
    for corner, folders, names in os.walk(folder):
        for name in names:
            DeleteFile(FullPath(corner, name))
    for corner, folders, names in os.walk(folder, topdown=False):
        if len(os.listdir(corner)) == 0:
            if corner != folder:
                os.rmdir(corner)


def openFile(filepath):
    if platform.system() == 'Darwin':  # macOS
        subprocess.call(('open', filepath))
    elif platform.system() == 'Windows':  # Windows
        os.startfile(filepath)
    else:  # linux variants
        subprocess.call(('xdg-open', filepath))