# Created by wu on 2017/8/31
from textwrap import dedent
import shutil
from pathlib import Path

def _copy(self, target):
    """
    Copy file to target file
    :param self:
    :param target:
    :return:
    """
    assert self.is_file()
    shutil.copy(str(self), str(target))

def _write_txt(self, data, encoding='utf-8', errors=None, mode=None, **kwargs):
    """
    Open the file in text mode, write to it, and close the file.
    """
    if not isinstance(data, str):
        raise TypeError('data must be str, not %s' %
                        data.__class__.__name__)
    from jinja2 import Template
    data = dedent(data)
    temp = Template(data)
    data = temp.render(**kwargs)
    data = data.format(**kwargs)
    with self.open(mode=mode, encoding=encoding, errors=errors) as f:
        return f.write(data)


def _insert_name(self, name):
    stem = self.stem
    suffix = self.suffix
    return Path('%s_%s%s' % (stem, name, suffix))


def add(self, path):
    return Path((str(self) + str(path)))


def rsplit(self, *args, **kwargs):
    return str(self).rsplit(*args, **kwargs)


def endswith(self, k):
    return str(self).endswith(k)


Path.copy = _copy
Path.rst_write = _write_txt
Path.insert = _insert_name
Path.__add__ = add
Path.rsplit = rsplit
Path.endswith = endswith
