.. highlight:: shell

============
贡献指南
============

我们欢迎贡献，即使再小的贡献对于taskserver项目都是非常有益。

参与taskserver项目有多种方式:

贡献的类型
----------------------

报告Bug
~~~~~~~~~~~

Report bugs at http://gitlab.oebiotech.cn/oe/taskserver/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

修复Bug
~~~~~~~~

Look through the GitLab issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

实现功能
~~~~~~~~~~~~~~~~~~

Look through the GitLab issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

编写文档
~~~~~~~~~~~~~~~~~~~

taskserver could always use more documentation, whether as part of the
official taskserver docs, in docstrings, or even on the web in blog posts,
articles, and such.

提交反馈
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at http://gitlab.oebiotech.cn/oe/taskserver/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

如何开始参与？
---------------

Ready to contribute? Here's how to set up `taskserver` for local development.

0. 前置条件：
    - 安装好git，若通过 `git --version` 可以看到所安装的 git版本，若提示找不到git，请确认git的执行路径或者 安装git_
    - 安装好python， 建议python3.7以上版本。
1. Fork the `taskserver` repo on GitLab(http://gitlab.oebiotech.cn).::

    项目地址位于http://gitlab.oebiotech.cn/oe/taskserver

2. Clone your fork locally::

    $ git clone http://gitlab.oebiotech.cn/wuxuebiao/taskserver.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv taskserver
    $ cd taskserver/
    $ python setup.py develop  # 该命令将在对应的虚拟环境中以开发者模式安装taskserver

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox::

    $ flake8 taskserver tests  # 测试代码符合PEP8规范
    $ python setup.py test or pytest  # 对编写的代码进行测试
    $ tox  # 在多个环境下进行代码兼容性测试，默认测试python3.7和python3.8

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to Gitlab(http://gitlab.oebiotech.cn)::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a Merge request through the GitLab website.

    <http://gitlab.oebiotech.cn/wuxuebiao/taskserver/merge_requests/new>


.. _安装git: https://git-scm.com/downloads

Merge Request 指导
-------------------------

在提交一个合并请求前，请核对一下代码符合下述情况：
Before you submit a Merge request, check that it meets these guidelines:

1. 提交的合并请求应该包含测试，且应尽量提高代码的覆盖度
2. 如果你的合并请求增加了功能，那么应该及时更新对应的文档。把你新增加的功能放入包含文档说明的函数中，然后在README.rst中加入功能特性列表
3. 提交的合并请求应该兼容python3.7+版本

提示
-------

To run a subset of tests::

$ pytest tests.test_taskserver


部署
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run::

    $ bump2version patch # possible: major / minor / patch
    $ git push
    $ git push --tags

上述代码默认部署至PyPI，部署至服务器上的请和管理员联系。
