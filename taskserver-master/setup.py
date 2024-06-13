#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages, glob

with open('README.rst', encoding='utf-8') as readme_file:
    readme = readme_file.read()

with open('CHANGELOG', encoding='utf-8') as history_file:
    history = history_file.read()

# 需要安装的依赖
requirements = ['fire', 'fastapi', 'uvicorn', 'cairosvg ', 'pandas', 'selenium', 'apscheduler', 'esdk-obs-python', ]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', 'pytest-mock' ]

# setup的更多细节配置见 https://setuptools.readthedocs.io/en/latest/
setup(
    author="OE Biotech Team",
    author_email='xuebiao.wu@oebiotech.com',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="qsub 投递的任务守护进程及web服务",
    entry_points={
        'console_scripts': [
            'taskserver=taskserver.cli:main',
        ],
    },
    install_requires=requirements,
    long_description=readme + '\n\n' + history,
    long_description_content_type='rst',
    # 该项与MANIFEST.in文件配套使用，docs和tests文件夹下的所有文件均不会打包发布，
    # 其他非明确指定的文件自动在包发布的时候包含到包中
    include_package_data=True,
    keywords='taskserver, OE Biotech',
    name='taskserver',
    packages=find_packages(include=['taskserver', 'taskserver.*']),
    # 在安装的时候自动将Scripts下的脚本注入当前Python的bin或Scripts目录作为可执行程序。
    packages_data={'taskserver': ['datasets/*']},
    # package_dir={'taskserver': 'src/mypkg'},
    scripts=glob.glob('scripts/*', recursive=True),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='http://gitlab.oebiotech.cn/oe/taskserver',
    version='0.2.0',
    zip_safe=False,
)
