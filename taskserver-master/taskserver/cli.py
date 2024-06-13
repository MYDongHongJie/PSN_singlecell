"""Console script for taskserver."""
import fire
import uvicorn
LOGGING_CONFIG: dict = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "default": {
            "()": "uvicorn.logging.DefaultFormatter",
            "fmt": "%(asctime)s %(levelprefix)s %(message)s",
            "use_colors": None,
        },
        "access": {
            "()": "uvicorn.logging.AccessFormatter",
            "fmt": '%(asctime)s %(levelprefix)s %(client_addr)s - "%(request_line)s" %(status_code)s',  # noqa: E501
        },
    },
    "handlers": {
        "default": {
            "formatter": "default",
            "class": "logging.StreamHandler",
            "stream": "ext://sys.stderr",
        },
        "access": {
            "formatter": "access",
            "class": "logging.StreamHandler",
            "stream": "ext://sys.stdout",
        },
    },
    "loggers": {
        "taskserver": {'handlers': ['default'], 'level': 'DEBUG'},
        "uvicorn": {"handlers": ["default"], "level": "INFO"},
        "uvicorn.error": {"level": "INFO"},
        "uvicorn.access": {"handlers": ["access"], "level": "INFO", "propagate": False},
    },
}

# from taskserver import core
# from taskserver.oesc_task import oesc_task


def run_app(host='0.0.0.0', port=8888, **kwargs):
    from . import settings
    kwargs = settings.update(kwargs)
    kwargs['log_config'] = LOGGING_CONFIG
    uvicorn.run("taskserver.core:app", host=host, port=port, **kwargs)


def main():
    fire.Fire(run_app)
