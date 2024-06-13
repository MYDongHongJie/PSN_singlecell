import os


def get(key, default=None):
    return os.environ.get(key, default=default)


SC_MODULE = get('SC_MODULE', 'OESingleCell/v_3.0.0_scRNAcloud')
JOB_INTERVAL = get('SC_JOB_INTERVAL', 10)
SCC_AK = get('SCC_OBS_AK', None)
SCC_SK = get('SCC_OBS_SK', None)
SCC_SERVER = get('SCC_OBS_SERVER', None)
SCC_BN = get('SCC_OBS_BUCKET_NAME', '')


def update(kwargs):
    """
    用于将命令行的参数值更新至settings.py文件中。
    需要注意：命令行传入的参数不区分大小写，比如命令行传入--sc_module=ABC 和--SC_MODULE=ABC是一样的。
    :param kwargs: 命令行参数键值对
    :return: 访问不是命令行参数的键值对
    """
    new = kwargs.copy()
    for k, v in kwargs.items():
        attr = k.upper()
        if attr in globals():
            globals()[attr] = v
            new.pop(k)
    return new
