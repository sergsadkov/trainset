import time
from .paths import *


__temp_dir__ = os.path.join(os.environ['TMP'], 'image_processor')


def TempName(name = 'tmp', ext = ''):
    i = 0
    temp_folder = f'{__temp_dir__}\\{i}'
    while os.path.exists(temp_folder):
        i += 1
        temp_folder = f'{__temp_dir__}\\{i}'
    os.makedirs(temp_folder)
    return FullPath(temp_folder, name, ext)


class TempFiles:

    def __init__(self):
        self.delete = []
        self.delshp = []
        self.copy = {}
        self.copyshp = {}

    def CheckArg(self, arg):
        if isinstance(arg, str):
            if re.search(r'^\\+172\.', arg):
                name = os.path.split(arg)[1]
                if name:
                    new_arg = TempName(name)
                    if os.path.exists(arg):
                        if name.lower().endswith('.shp'):
                            CopySHP(arg, new_arg)
                            self.delshp.append(new_arg)
                        else:
                            CopyFile(arg, new_arg)
                            self.delete.append(new_arg)
                    elif name.lower().endswith('.shp'):
                        self.copyshp[new_arg] = arg
                    else:
                        self.copy[new_arg] = arg
                    return new_arg
        elif isinstance(arg, (tuple, list)):
            to_list = isinstance(arg, list)
            args, kwargs = self.CheckArgsKwargs(*tuple(arg), **{})
            if to_list:
                args = list(args)
            return args
        return arg

    def CheckArgsKwargs(self, *args, **kwargs):
        args_list = list(args)
        kwargs_dict = dict(kwargs)
        for i, arg in enumerate(args_list):
            args_list[i] = self.CheckArg(arg)
        for key in kwargs_dict:
            arg = kwargs_dict[key]
            kwargs_dict[key] = self.CheckArg(arg)
        return tuple(args_list), dict(kwargs_dict)

    def __del__(self):
        for new_arg in self.copy:
            CopyFile(new_arg, self.copy[new_arg])
            DeleteFile(new_arg)
        for new_arg in self.copyshp:
            CopySHP(new_arg, self.copyshp[new_arg])
            DeleteSHP(new_arg)
        for new_arg in self.delete:
            DeleteFile(new_arg)
        for new_arg in self.delshp:
            DeleteSHP(new_arg)


def StopFromStorage(func):
    def wrapped(*args, **kwargs):
        temp_files = TempFiles()
        args, kwargs = temp_files.CheckArgsKwargs(*args, **kwargs)
        res = func(*args, **kwargs)
        del temp_files
        return res
    return wrapped


if not os.path.exists(__temp_dir__):
    os.makedirs(__temp_dir__)

# Delete all files with last change over 1 day ago from __temp_dir__
for corner, folders, names in os.walk(__temp_dir__):
    for name in names:
        file = FullPath(corner, name)
        if (time.time()-os.path.getmtime(file))/86400 > 1:
            DeleteFile(file)
    if len(os.listdir(corner))==0:
        if corner != __temp_dir__:
            os.rmdir(corner)