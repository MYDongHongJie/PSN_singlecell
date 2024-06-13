import requests
import os
from oebio.utils.log import getLogger

logger = getLogger('oe.cloud.sc.qsub')


# response.text
# 0, "SUCCESS" // 回调成功
# 201, "记录不存在"  // 回调了非法的任务id
# 311, "未找到output.json文件" // 任务状态回调为"success", 但未检测到对应output.json文件
# 312, "检测到同一状态重复回调" // 同一状态回调，如果正确响应过一次，云平台不会再重复接收对应任务相同状态的回调
# 313, "任务状态异常,未检测到上一个状态回调" // 状态回调需按顺序回调，不可跳跃：比如 任务还没回调"start"， 就回调"success" or "fail"
# 500, "系统错误" // 云平台内部错误
def api_upload(taskid, state, i):
    """
    :param taskid: taskid
    :param state: start,success,fail
    :param i:
    :return:
    """
    api_url = os.environ.get('SCC_CALLBACK')
    logger.info(f'upload task status to {api_url}: {taskid}:{state}')
    # api_url = f"http://sc.oebiotech.com/api/v1/tasks/callback"
    # or "http://122.112.243.140:8000/api/v1/tasks/callback"
    for num in range(1, i + 1):
        try:
            resp = requests.post(api_url, json={'task_id': taskid, 'state': state})
        except requests.exceptions.RequestException as e:
            logger.warning(f'第{num}次重试。重试原因：{e}')
            continue
        logger.info(f'code:{resp.status_code}, message: {resp.text}')
        break
