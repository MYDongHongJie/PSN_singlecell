FROM python:3.8-slim

LABEL description="qsub 投递的任务守护进程及web服务"
LABEL maintainer="OE Biotech Team(xuebiao.wu@oebiotech.com)"

RUN mkdir /app
WORKDIR /app

COPY . /app
RUN pip install -e . -i https://pypi.tuna.tsinghua.edu.cn/simple

CMD ["taskserver"]
