import logging
from functools import lru_cache


from taskserver import settings

logger = logging.getLogger(__name__)


@lru_cache
def get_client(bucket_name=True):
    from obs import ObsClient
    client = ObsClient(access_key_id=settings.SCC_AK, secret_access_key=settings.SCC_SK, server=settings.SCC_SERVER)
    return client.bucketClient(settings.SCC_BN) if bucket_name else client


def format_key(key):
    name = settings.SCC_BN
    if name in key:
        index = key.find(name) + len(name) + 1
        key = key[index:]
    return key


def obs_method(method, *args, **kwargs):
    client = get_client()
    logger.info('%s开始执行，相关参数：%s, %s', method, *args, **kwargs)
    resp = getattr(client, method)(*args, **kwargs)
    if resp.status >= 300:
        logger.warning('%s执行异常,反馈原因%s,执行参数%s, %s', method, resp.reason, args, kwargs)
        return True
    logger.info('%s执行成功', method)


def list_objects(source, marker=None):
    """
    Args:
        source:OBS中要下载的文件夹
        marker:本地指定的文件夹
    Returns: 返回对象包含以下字段
        key str 对象名。
        lastModified str 对象最近一次被修改的时间。
        etag str 对象的MD5值（当对象是服务端加密的对象时，etag值不是对象的MD5值）。
        size int对象的字节数。
        owner Owner 对象的所有者。
        storageClass str 对象的存储类型。
        isAppendable bool 对象是否可被追加上传。
    """
    client = get_client()
    body = client.listObjects(source, marker=marker or source).body
    for obj in body.contents:
        if obj.get('key')[-1] != '/':
            yield obj
    if body.is_truncated:
        yield from list_objects(source, marker=body.next_marker)
