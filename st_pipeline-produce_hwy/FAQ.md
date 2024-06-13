
# FAQ
## 1、流程引擎 *snakemake* 常见报错说明
### 1.1 ：WorkflowError:
Target rules may not contain wildcards. Please specify concrete files or a rule without wildcards.

直译过来是：目标规则可能不包含通配符，请指定没有通配符的具体文件或规则。
看似是通配符的使用问题，实则和 rule all 的 input 错误有关。

Snakemake 的rule all 需要输入一个列表，这个列表包括所有的步骤会产生的文件（也就是output）。实际上，Snakemake的工作机制是先通过 shell 产生 output 文件，再检查这些 output 文件是否在指定的目录下存在，确认存在后再执行下一个 rule 。这个检查过程就是通过寻找rule all 下的input 列表完成的。

因此，在 input 列表没写清楚的情况下，如果在 rule 中还使用的通配符，就会发生相应的报错。改正方法就是修改你的 input 列表，确认每一个 rule 输出的文件和他们的目录都存在于input 列表下。

### 1.2： Nothing to be done.

和上一个问题相同，也是由于 input 列表没写好造成的，但是这里没有使用通配符。解决方法也一样，重新检查 input 列表。

### 1.3：Error:

Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory：

翻译：无法锁定目录。请确保没有其他Snakemake进程尝试在下面的目录中创建相同的文件

问题所在：之前这个目录下已经执行过其他的 snakemake 命令了。

解决办法也直接给出来了：

If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.

在这个目录下执行 ：

snakemake --unlock
### 1.4：Missing files after 5 seconds:

This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait.

翻译：这可能是由于文件系统延迟造成的。如果是这样，考虑使用延迟等待来增加等待时间。

问题所在：

1.如报错信息所说，确实是因为文件的产生没有在5秒（默认等待时间）内产生，可以使用 --latency-wait 增加等待时间

2.输出的文件并没有输出到指定的目录，修改为正确的输出目录即可

很重要的一点：只用伪执行来检查是不行的！！！伪执行只能检查相应的文件和目录，输出时的问题根本看不出来。

