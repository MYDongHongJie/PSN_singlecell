.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## 删除所有build, test, coverage命令生成的额外文件以及Python artifacts文件

clean-build: ## 删除 build命令生成的额外文件
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## 删除 Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## 删除 test and coverage 生成的额外文件
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

lint: ## 利用 flake8 检查代码风格，主要检测taskserver和tests文件夹下的python代码
	flake8 taskserver tests

test: ## 运行默认的测试
	pytest

test-all: ## 在pyhton3.7 python3.8环境下测试，利用flake8检查代码风格
	tox

coverage: ## 利用默认的python检测代码覆盖度
	coverage run --source taskserver -m pytest
	coverage report -m
	coverage html
	$(BROWSER) htmlcov/index.html

docs: ## 生成 Sphinx 文档，包含项目的API文档
	rm -f docs/source/taskserver.rst
	rm -f docs/source/modules.rst
	sphinx-apidoc -o docs/source taskserver
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes（未测试）
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

release: dist ##package and upload a release（）
	twine upload dist/*

dist: clean ## 构建源码包和 whl编译包
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean ## 将项目安装到python的 site-packages 目录下
	python setup.py install
