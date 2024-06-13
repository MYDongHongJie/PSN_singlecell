.. highlight:: shell

============
安装
============


稳定版安装
--------------

To install taskserver, run this command in your terminal:

.. code-block:: console

    $ pip install taskserver  ## 若包不在github上，则无法安装， 请参考从源码安装

This is the preferred method to install taskserver, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


源码安装
------------

taskserver 项目源码可以从 <http://gitlab.oebiotech.cn/wuxuebiao/taskserver> 下载.

PS: OE gitlab访问需要有对应的账户，以及对应的 oe组的权限。若不满足要求，请与gitlab 管理员或者主管联系。

你可以从http://gitlab.oebiotech.cn/oe clone 代码:

.. code-block:: console

    $ git clone http://gitlab.oebiotech.cn/oe/taskserver.git

从源码clone下来之后，切换到taskserver目录中进行安装：

.. code-block:: console

    $ cd taskserver
    $ python setup.py install

taskserver使用方法请参考 :ref:`taskserver_usage`。

