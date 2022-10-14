# ***************************************************************************
# Place init_logger() at the start of the program.
# Wrap function with @log_me to include it in logging.
# ***************************************************************************

import logging
import os.path
import shutil
import socket
import sys
from datetime import datetime


def init_logger(filename='razmetka.log',
                # encoding='utf-8',
                level=logging.NOTSET,
                format='%(asctime)s %(levelname)s %(message)s\n',
                datefmt='%m/%d/%Y %I:%M:%S %p'):

    logging.basicConfig(filename=filename,
                        # encoding=encoding,
                        level=level,
                        format=format,
                        datefmt=datefmt)
    logging.info('\nStarted' + '.' * 100)


def log_me(function):
    def wrapper(*args, **kwargs):
        # logging.debug(f'Function {function.__name__} args: {str(args)} kwargs: {str(kwargs)}')
        name = function.__name__
        args_kwargs = str_args_kwargs(*args, **kwargs)
        logging.debug(f'Start {name}({args_kwargs})')
        try:
            return function(*args, **kwargs)
        except Exception as e:
            # logging.error("Exception ", exc_info=True)
            logging.error(f'Exception {name}({args_kwargs}) -> {str(e)}')
    return wrapper


def upload_log(log_path=r'\\172.21.195.2\thematic\Log'):
    try:
        file = get_device_id() + '_' + str(datetime.now().strftime("%Y-%m-%d_%H-%M-%S")) + '.log'
        shutil.copyfile('razmetka.log', os.path.join(log_path, file))
    except Exception as e:
        logging.error("Exception ", exc_info=True)


def get_device_id():
    os_type = sys.platform.lower()
    if "win" in os_type:
        command = "wmic bios get serialnumber"
    elif "linux" in os_type:
        command = "hal-get-property --udi /org/freedesktop/Hal/devices/computer --key system.hardware.uuid"
    elif "darwin" in os_type:
        command = "ioreg -l | grep IOPlatformSerialNumber"
    return os.popen(command).read().replace("\n", "").replace("SerialNumber", "").replace(" ", "")


def get_ip():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
    s.connect(('<broadcast>', 0))
    return s.getsockname()[0]


def str_args_kwargs(*args, **kwargs):
    return_string = ''
    if isinstance(args, tuple):
        return_string += '", "'.join([str(arg) for arg in args])
    if isinstance(kwargs, dict):
        return_string += '", "'.join([str(key)+'='+str(kwargs[key]) for key in kwargs])
    return return_string
