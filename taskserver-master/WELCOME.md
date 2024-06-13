# 欢迎来到欧易生物Python的世界

也许你觉得从头建立一个python的包很麻烦， 需要设置的地方非常多。幸运的是，你遇到了目前这个包。
[toc]

## 使用方式

在看到这个文档之前你应该生成taskserver。现在请切换工作目录到
taskserver中。

### 前言

在使用之前，若对各个文件夹及文档的说明不清楚，请参考后续的详细介绍文件夹和文件区域。以下是关于taskserver的使用方向。

### 构建开发环境

一般每个项目构建一次即可，后续可以一直使用。
```shell
# 核对可执行路径中的所有python
which -a python
# 选择合适的python版本建立虚拟环境，建议选择python3.7以上版本,
python -V  # 查看python版本，本次显示是python3.7
# -m 表示调用venv python包，第二个venv表示在当前文件夹下创建一个venv37的虚拟环境
python -m venv venv37
# windows下 激活venv37环境的方法
venv37/Scripts/activate.bat
# linux或者mac下激活venv37环境的方法
source venv37/bin/activate
# 再次核对python的执行路径，可以发现python的执行路径为venv37中的python，后续提到的python以及venv37均指本次构建的路径。
which python
# 安装开发依赖相关的包和工具 -i 表示使用对应的镜像加速安装
pip install -r requirement_dev.txt -i  https://pypi.tuna.tsinghua.edu.cn/simple
# 将当前taskserver以开发模式安装到venv37中, 注意最后的点号。
pip install -e .
# 开始编写对应的代码，利用pycharm或者VScode打开当前的项目，将对应的python设置为venv37中的Python，在taskserver文件夹中编写代码， 在tests文件夹编写对应的测试内容。
```
### 代码编写

代码编写部分，请根据需要进行代码编写。有以下几个指导说明：

1. 代码建议以函数或者类的方式进行编写，不应在文件的顶层作用域进行代码编写
	```python
	# 不建议
	import pandas as pd

	data = pd.read_table('some_file.txt')
	# 对data对象进行数据处理。
	data.to_excel('result.xlsx')

	# 建议的方式
	def format_data(infile, output='result.xlsx'): # 给对应的操作进行命名，便于理解和调用
	    """对txt数据进行格式化，保存为excel格式"""   # 对函数的功能进行介绍，若可能的话，对每个参数进行说明
	    data = pd.read_table(infile)
	    # 对数据进行处理
	    data.to_excel(output)
	```
2. 若代码功能比较多，建议按照模块的方式形成独立的文件，一方面，独立的分析在进行查看代码的时候，便于查看。另外，独立的文件在进行git 代码合并的时候，引起的冲突会少一些。
3. 若模块中的功能是比较重要的，建议在模块层面增加文档描述，这些描述的文档在后续文档生成的时候，会自动提取出来，作为项目的文档描述。


### 测试编写

测试的编写是一个非常大且重要的领域，但是测试的编写也是需要很长时间的实践。所以目前未明确规定测试的编写类型，不过taskserver项目中默认内置了一些测试框架（pytest和 unittest），其中unittest是python系统自带的测试框架，而pytest是一个新兴的框架，且pytest可以兼容unittest框架。所以建议如果入门测试编写的话，直接学习pytest即可。

测试的重要性：

1. 可以保障代码运行是正确的；
2. 可以让其他在进行代码重构的时候不用担心重构出错了；
3. 在很多自动化测试过程以及代码质量管理中，测试是一个很重要的组成部分，下一节的代码覆盖度检测就是用编写的测试对代码的覆盖度进行检测；
4. 测试实际上也是一种流程重现的方式。

测试的编写提示：

1. 每个函数均是独立的测试；
2. 若函数中存在一个if或者try语句，则说明条件不一样，一般也是要额外编写独立的测试；

### taskserver项目代码覆盖度检测

代码的覆盖度检测的目的，主要用于查看编写的代码是否满足且仅满足特定的需求。

举个栗子，假设代码编写了1000行，介绍代码可以检测A，B，C，D各种环境，但是实际上，发布出来的功能仅仅能适应A环境，起作用的区域也只有A功能的实现的区域，其他区域都是干扰项。在编写代码的过程中，有时候是不知道编写的代码是否在运行过程中有包含，通过测试调用taskserver代码，就可以衡量代码是否可以在调用过程中均产生贡献了。

代码覆盖度测试输出命令
```shell
make coverage
```

注意：若测试不通过是无法生成代码覆盖度测试报告！

### 重新生成命令行入口

命令行入口的设置在项目文件夹中的setup.py文件中，默认设置了命令taskserver 与 taskserver文件夹下的cli.py文件的main函数进行关联。
若需要增加其他命令入口，可以参考setup中的设置进行调整。
如需要增加一个`taskserver_plot 命令`与 taskserver的plot.py中的plot
函数关联，可以在setup.py文件的`setup`函数，在`console_scripts`的列表项中额外增加一项
`"taskserver_test=taskserver.plot:plot"`后，重新运行`pip install
 -e .` 即可生成可以在命令行中调用的命令`taskserver_plot`。

上述结果的修改示例：

**修改前**
```python
entry_points = {
    'console_scripts': [
        'taskserver=taskserver.cli:main',
    ]
}
```

**修改后**
```python
entry_points = {
    'console_scripts': [
        'taskserver=taskserver.cli:main',
        'taskserver_plot=taskserver.plot:plot',
    ]
}
```

修改setup.py文件之后，重新运行`pip install -e .`， 之后即可在命令行中使用命令`taskserver_plot`



需要注意的是：
>命令放的位置是与调用pip所使用的python环境一致，通常用虚拟环境进行开发的时候，记得最开始要激活对应的环境，不然允许对应的命令会提示找不到对应命令。

提示：
>为了避免对应的Python环境所依赖的包对其他人的造成影响，可以考虑安装的时候建立虚拟环境安装，然后将对应的命令行软连接到常用的可执行路径下。这样可以非常方便的隔离不同环境，又满足自己个性化的需求。

### 对项目进行GIT操作

在生成taskserver包的时候， 默认已经进行了第一次提交。所以如果你对taskserver项目中的大部分文件进行修改的时候，在命令行下可以利用git命令进行查看修改内容。

举个栗子，在taskserver中的taskserver文件夹，修改了源码文件夹中的taskserver.py中的main函数。

按以下流程：

```shell
#shell 命令新建一个分支
git checkout -b new_function
```

更改taskserver文件夹的core.py中代码，从
```python

def main():
    """该功能函数作为功能实现的入口"""
    raise NotImplementedError
```
变为

```python
def main():
    """该功能函数作为功能实现的入口"""
    print('Hello OE Biotech!')

```

在命令行下查看使用git查看
```shell
git status # 输出如下
#位于分支 master
#尚未暂存以备提交的变更：
#（使用 "git add <文件>..." 更新要提交的内容）
#（使用 "git restore <文件>..." 丢弃工作区的改动）
#修改：     taskserver/core.py

#查看改动的文件，请注意核对有没有意外的文件被修改了。
git diff
#将修改的taskserver/core.py加入暂存区
git add taskserver/core.py
#为本次提交生成一个注释说明，可以参考用标准的commit语句。
git commit -m "增加了一个打印到屏幕的提示功能"

#在http://gitlab.oebiotech.cn上建立了taskserver的项目包。
#taskserver项目在建立的时候会自动添加两个远程分支，查看远程分支信息
git remote -v   # 输出如下：
#oe	http://gitlab.oebiotech.cn/oe/taskserver.git (fetch)
#oe	http://gitlab.oebiotech.cn/oe/taskserver.git (push)
#wuxuebiao	http://gitlab.oebiotech.cn/wuxuebiao/taskserver.git (fetch)
#wuxuebiao	http://gitlab.oebiotech.cn/wuxuebiao/taskserver.git (push)

#远程仓库 oe 表示公共的仓库，一般表示各个平台的主仓库，该仓库由各个平台主管管理，可以有权限限制。当然，在gitlab
上也个人也可以建立组（group）。
#远程仓库 wuxuebiao表示gitlab
上个人的仓库，个人具有所有的修改权限，一般工作的时候，将本地的修改或者功能推送到远程的个人仓库，然后从个人仓库提出合并请求，申请代码审核。最后由项目负责人审核通过，合并到公共或者平台仓库中。

#推送本地修改的功能分支至远程个人分支
git push -u wuxuebiao new_function
#后续请在gitlab上进行操作，提出合并请求，并经代码review。

```

### 代码风格检测

统一代码风格重要性好比统一语言标准，这里的语言指的是中文，英文，俄文之类的，如果大家都说不一样的语言，沟通成本就非常高。即使是中文，可能不同地区的语言也是存在不一样的，有可能有专用的俚语，那么这些词语在沟通的时候也会额外带来理解成本。这也是咱们大中华的普通话推广的重要意义。

在编程世界中，也存在着各种不同语言，同一种语言不同人编写，可能的展现方式就不一样；同时，不同语言的统一标准不一样。幸运的是，在Python里面有一个比较统一的语言标准 [PEP8](https://www.python.org/dev/peps/pep-0008/) [中文](https://www.jianshu.com/p/26dde1315c74),我们在按照同样的标准来编写代码，沟通才会更高效。但是代码编写风格是一种习惯养成，如果刚开始，我们就要求代码规范中的几十个内容都记住，可能比较耗时。所以，在taskserver项目生成的时候，我就内置了[flake8](flake8.pycqa.org/ ) 检测。在每次利用git进行提交的时候，就自动进行代码风格检测，如果代码风格检测失败，则无法提交代码。

当然，在编写代码的过程中，也可以随时随地的检测风格是否符合规范。再举个栗子:

```shell
# 在shell运行时，请务必核对所用的python环境是否符合预期。
# 在文件taskserver.py 中的main函数前额外增加3个空行
flake8 # 输出如下：
# ./taskserver/taskserver.py:7:1: E303 too many blank lines (5)

# 错误说明，说明在taskserver.py文件中的第7行位置处有多余的空行，错误代码为E303，括号内的数字表示空行的数量。（PEP8中，顶层函数之间的距离建议空两行）。将错误修复之后，再重新运行 flake8 则不报错。

```
**注意：在git提交的时候，额外增加了git commit 的hook，在每次提交的时候均需要flake8检测没有错误才能提交。所以，作为taskserver项目的开发者，建议将这个设置保留下来，便于统一语言。**

tips：

1. 很多编辑器，如Pycharm在编写代码的时候，都可以自动进行格式化，其他工具也一样，可以通过该工具格式化代码，当然我建议是将这种代码风格形成一种习惯。
2. 待补充

###  taskserver 项目文档生成

taskserver 项目的文档是利用 Sphinx进行文档生成。
```shell
#根据项目中的makefile，直接输入下述命令即可生成对应的网页文档（linux下可行，windows下未测试）。
make docs
```

### taskserver 项目打包

打包成源码包或者whl编译的包，这些包可以在其他电脑上通过pip install  的方式进行安装。
```shell
make dist
```



## 详细介绍版(项目中文件夹)

分别介绍taskserver中每个文件夹的意义以及使用方式。

### taskserver 文件夹（源码文件）

该文件夹是python包的源码包，包内文件说明如下：

1. __init__.py 文件，文件中包含了作者、邮箱和版本信息，同时，若有需要将一些python对象设置为导入taskserver包即可使用，可以在这个文件中将对应的python对象（函数，类等）导入。
	举个栗子：源码包中的taskserver.py文件中包含了一个 `analysis` 函数，若想设置未导入taskserver即可使用，则可以在`__init__.py`文件中加入以下代码

```python
	from .taskserver import analysis


	__all__ = ['analysis', ]
```
	导入之后即可在其他地方使用` taskserver.analysis`的方式调用该方法。

2. cli.py 文件。 该文件是实现自动形成命令行借口的关键文件。 该文件中的`main`与项目文件夹中的`setup.py`文件中的`console_scripts`
是关联的，即一旦安装了taskserver包之后，即可实现从命令行调用 `taskserver`命令的方法。 命令行的编写，根据项目初始化时候的配置，可以选择`fire`、`click` 或者 `argparse`包来编写，推荐用fire，编写起来比较灵活。
3. core.py文件，该文件主要作为一个本项目包功能实现的主要入口，可以灵活运用python包的各种概念扩展该文件。

### docs文件夹 （文档文件）

该文件夹主要实现taskserver项目的文档内容。其中文档主要的书写方式是以 restructtext(.rst)格式为主，也可以兼容 Markdown格式。

- templates 是文档模板文件夹，用于存放一些文档模板，目前包含了一份 issue 的模板。
- themes文件中包含的是欧易生物自定义主题的文件夹。
- 文件夹中的authors.rst,  contributing.rst, history.rst, readme.rst 文件与外部的文件关联，这些文件若有需要调整的话， 是需要项目主路径的文件进行调整。
- conf.py 是 python `sphinx` 包生成文档的配置文件， 里面已经进行了一些OE生物的文档模板定制配置。
- make.bat和 Makefile分别是 windows下和 unix下生成项目文档的方式， 如:
	```shell
	make html # 直接在当前目录生成 _build 目录， 目录中包含对应格式的文档。
	```
- installation.rst 文件中主要介绍了如何从源码或者从代码包中安装taskserver项目的操作方法，如果该项目有特殊安装配置，也请将installation.rst文件补充完整。
- **usage.rst taskserver 包的使用方法说明**。若包的使用方法较多，可以将该文件整理成文件夹，利用 [sphinx](http://www.sphinx-doc.org/en/stable/) 的语法进行文档编写。 更复杂的方式可以参考 `sphinx.ext.autodoc`( [详细链接](http://www.sphinx-doc.org/en/stable/usage/extensions/autodoc.html )) 直接从代码中生成对应的文档。

### Scripts(非python脚本路径)

该文件夹下用于放置一些非python的脚本或者可执行程序，在安装taskserver后，这类脚本会自动注入到当前python的 Scripts(Windows)、bin(Unix)路径下。 当然，如果你对taskserver进行封装后，形成常用的python脚本，也可以放到该文件夹下，linux环境下也可以直接调用这类脚本。

需要注意是： 在执行 `pip install taskserver` 后， Scripts下的所有文件的中的第一行以 `#!`开始，且包含 `python`的字段，则会在安装后自动替换成当前python的执行路径，同时增加版本依赖限制，从而避免用默认的python执行导致的意外错误。

### tests( 测试内容文档)

该文件夹中可以包含功能测试和单元测试，默认没有区分。若存在需要区分功能测试和单元测试的情况，建议可以在当前文件夹新建 func 和unit的文件家。

- func 功能测试用于测试完整的流程。保障流程的顺利运行。
- unit 单元测试用于测试taskserver包中的函数功能实现。

文件夹中初始包含了一个 test_taskserver.py 的文件，文件中根据选择的测试框架不同，简单实现了命令行借口的一个测试函数。

更多的测试内容的编写可以参考：

- pytest <https://docs.pytest.org/en/latest/>
- unittest <https://docs.python.org/zh-cn/3.7/library/unittest.html>

由于pytest的写法比较简单和友好，格式限制没有 unitest那么严格，建议优先选择pytest，另外，pytest测试框架可以兼容 unittest和doctest的测试内容。

==测试一门艺术==

## 详细介绍（项目文件）

本区域主要介绍taskserver 项目中的相关文件说明。

### AUTHORS.rst

 该文件中主要包含taskserver项目的相关人员。比如核心开发团队，贡献者等。 一旦有其他人对本项目进行贡献，那么应该将其名字及联系方式也放在上面，便于后续贡献者记录。

 ### CONTRIBUTING.rst

 该文件包含以下几个方面的描述：

 1. 如何对taskserver项目进行贡献的方法，如通过报告Bug，修复bug，实现功能，编写文档，提交反馈等方法；
2. 包含了一个如何开始贡献的操作方法（共7步）；
3. 发起合并请求的指导说明；
4. 一些taskserver使用的小提示；
5. 部署方法（该方法适合软件部署人员，其他可以暂时不用考虑这个）。

### HISTORY.rst

taskserver包的版本更新历史。每次taskserver包更新，在这里提供更新记录。

### README.rst

关于taskserver包的功能说明。

### requirements_dev.txt

该文件作为开发者需要安装的包版本文件（通过 `pip install -r requirements_dev.txt`直接安装）， 正常使用包的依赖在 `setup.py`中已经指定。

### setup.cfg

关于taskserver包的一些设置文件。包含用于bumpversion，flake8等。

### setup.py

python 包的安装入口文件，可以使用的方式如下：

1. 直接在当前目录安装taskserver包 `pip install -e .`；
2. 利用setup.py文件将当前目录设置为安装目录，即修改代码的效果直接反映，不需要重新安装。`python setup.py develop`
3. 利用setup.py 分发可安装的whl文件，供其他人进行安装。`python setup.py bdist_wheel`,  该命令会生成taskserver-<版本号>-py2.py3-none-any.whl的安装包，若要指定使用的python版本和平台，请参考更详细的setup.py使用说明。
4. 利用setup.py 分发测试的whl文件，比如正式的预发布版本。`python setup.py egg_info -ba0 bdist_wheel`。
5. 利用setup.py生成当天时间戳的开发版本。 `python setup.py egg_info -b dev -d dbist_wheel`
6. 删除taskserver中dist文件夹中的多于的whl安装文件，仅保留2份最新的。`python setup.py rotate -m whl -k 2`。

### tox.ini

该文件是taskserver包在python3.7，python3.8下环境执行的配置文件。
命令行中输入`tox` 即可对taskserver执行python3.7，Python3.8环境下的测试代码。

换句话，意味着不同环境下的兼容性是通过 **测试代码** 来检测的。

同时，该配置文件还利用 [flake8](https://flake8.pycqa.org/en/latest/) 来执行代码格式化检测。 flake8的错误代码见<https://flake8.pycqa.org/en/latest/user/error-codes.html>

关于flake8 的几点说明：

1. flake8 是一款辅助检测 Python代码是否符合 PEP8 规范的工具， 通过pip安装之后可以直接从命令行调用；
2. taskserver包中的 flake8 配置更改了PEP8的 每行最长限制（从80更改为100）；
3. taskserver包的flake8 不检测以下文件及文件夹docs,.git,build,dist,docs/source/conf.py,venv,lib,.tox,.eggs；
4.  没用pylint是由于pylint限制太严格了。

### .editorconfig

taskserver包中的文件格式配置，具体如下：

- 全局配置：
	1. 所有文件的缩进风格为空格，每次缩进为4个空格；
	2. 每行最后不保留空格，以及文件最后插入空行；
	3. 文件编码默认为UTF-8, 换行格式为LF（即unix风格）；
- .bat文件：
	1. 缩进设置为tab；
	2. 换行设置为CRLF（windows风格）
- 所有的html，css，scss，json，yml，js文件：
	1.缩进风格为空格，每次缩进2个空格（由于这些文件容易深嵌套，4个空格会导致嵌套较深的区域缩进太深）。

### .gitignore

用于配置git的忽略文件。默认忽略python相关的无关文件，忽略的这些文件不会在git中产生修改记录。


### .makefile

该文件是将一些常用的命令打包，输入make可以查看对应的帮助信息。
