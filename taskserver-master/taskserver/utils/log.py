import inspect
import logging
import tempfile
import types
from collections import OrderedDict
from functools import wraps

def glob(regex, path='.', recursive=False):
    """
    遍历文件.

    Args:
        regex: 文件匹配规则，参考正则表达式
        path: 查找文件所在的文件夹，默认为当前路径，遍历当前文件夹符合regex的文件
        recursive: 遍历子文件夹，默认不遍历

    Returns:
        遍历对象生成器 即pathlib.Path类的对象集合。

    """
    pp = Path(path)
    func = pp.rglob if recursive else pp.glob
    for file in func(regex):
        yield file


class Log:
    """
    Logging decorator that allows you to log with a specific logger.

    """
    # customize these messages

    ENTRY_MESSAGE = '开始调用函数： %s'
    EXIT_MESSAGE = '%s 调用结束！！'

    def __init__(self, logger=None):
        self.logger = logger

    def __call__(self, func):
        """

        :param func:
        :return: a wrapper that wraps func.
         The wrapper will log the entry and exit points of
         the function with logging.INFO level.
        """
        if not self.logger:
            self.logger = logging.getLogger(func.__module__)

        @wraps(func)
        def wrapper(*args, **kwargs):
            self.log_befor(func, args, kwargs)
            self.logger.info(self.ENTRY_MESSAGE, func.__name__)
            f_result = func(*args, **kwargs)
            self.logger.info(self.EXIT_MESSAGE, func.__name__)
            return f_result

        return wrapper

    def log_befor(self, func, args, kwargs):
        fullargspec = inspect.getfullargspec(func)
        default_args = fullargspec.args
        default_values = fullargspec.defaults
        # 获取所有的默认参数及值
        default_map = OrderedDict(zip(reversed(default_args), reversed(default_values)))

        user_map = OrderedDict(zip(default_args, args))

        user_map.update(kwargs)
        func_name = func.__name__
        for i, arg in enumerate(default_args, 1):
            default = str(default_map.get(arg))
            user_input = user_map.get(arg)
            if len(str(user_input)) > 500:
                user = str(type(user_input)) + '对象'
            else:
                user = str(user_input)
            if default == user:
                self.logger.debug('%s 第%d个参数[%s]， 默认为`%s`，无修改', func_name, i, arg, default)
            else:
                self.logger.debug('%s 第%d个参数[%s]， 默认为`%s`，本次输入为 `%s`', func_name, i, arg, default,
                                  user)

    def __get__(self, instance, owner):
        if instance is None:
            return self
        else:
            return types.MethodType(self, instance)


def _checklogfile(log_dir='.'):
    logfiles = list(glob('oe_*.log', log_dir))
    if len(logfiles) > 0:
        return str(logfiles[0])
    else:
        _, log_tmp_fn = tempfile.mkstemp(suffix='.log', prefix='oe_', dir=log_dir)
        return log_tmp_fn


def init_log(logger, loglevel='INFO', log_dir='.'):
    """
    初始化logger对象
    Args:
        logger: 日志对象
        loglevel: console端的输出水平
        log_dir: 日志文件存储位置

    Returns:

    """
    # global log_tmp_dir, log_tmp_fn
    # log_tmp_dir = tempfile.mkdtemp(prefix='oebio', dir=log_dir)
    # log_tmp_fn = os.path.join(log_tmp_dir, 'oebio.log')

    logger.setLevel(getattr(logging, loglevel))

    console = logging.StreamHandler()
    console.setLevel(getattr(logging, loglevel))

    datefmt = '%H:%M:%S'
    if loglevel == 'DEBUG':
        debug_template = '[%(asctime)s][%(name)s][%(levelname)-1s|%(filename)s|%(lineno)d]' \
                         ': %(message)s'
        log_tmp_fn = _checklogfile(log_dir)
        file_handler = logging.FileHandler(log_tmp_fn, encoding='utf-8')
        file_handler.setLevel(getattr(logging, 'DEBUG'))
        file_handler.setFormatter(logging.Formatter(debug_template))
        logger.addHandler(file_handler)
        console.setFormatter(logging.Formatter(debug_template, datefmt=datefmt))
    else:
        info_template = '[%(asctime)s][%(levelname)-4s]%(name)s: %(message)s'
        console.setFormatter(logging.Formatter(info_template, datefmt=datefmt))

    logger.addHandler(console)


def getLogger(name=None, loglevel='DEBUG'):
    logger = logging.getLogger(name)
    # if not logger.hasHandlers():
    #     init_log(logger, loglevel=loglevel)
    return logger
