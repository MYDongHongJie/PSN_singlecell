##======================================================================================================================
if [ -f ~/.bashrc ]; then
        source  ~/.bashrc
fi
hostName=`hostname`
##======================================================================================================================
if [ "$hostName" == "snakemake-master" ]; then
    namespace="snakemake-scrna"  # snakemake-master
    # >>> conda initialize >>>
    # !! Contents within this block are managed by 'conda init' !!
    __conda_setup="$('/public/scrna/software/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
        eval "$__conda_setup"
    else
        if [ -f "/public/scrna/software/miniconda3/etc/profile.d/conda.sh" ]; then
            . "/public/scrna/software/miniconda3/etc/profile.d/conda.sh"
        else
            export PATH="/public/scrna/software/miniconda3/bin:$PATH"
        fi
    fi
    unset __conda_setup
    # <<< conda initialize <<<
    conda activate  /public/scrna/software/miniconda3/envs/snakemake
    alias snakemake="/public/scrna/software/miniconda3/envs/snakemake/bin/snakemake"

elif [ "$hostName" == "server02" ]; then
    namespace="k8s-scrna" # server02
    # >>> conda initialize >>>
    # !! Contents within this block are managed by 'conda init' !!
    __conda_setup="$('/data/software/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)" 
        if [ $? -eq 0 ]; then
            eval "$__conda_setup"
        else
            if [ -f "/data/software/miniconda/etc/profile.d/conda.sh" ]; then
                . "/data/software/miniconda/etc/profile.d/conda.sh"
            else
                export PATH="/data/software/miniconda/bin:$PATH"
            fi
        fi
        unset __conda_setup
    conda activate  /data/software/IND_dev/conda_envs/snakemake_v7.19.1
    ## 查看本地混合云的空闲资源情况
    alias gss="curl -s -d type=all -d region=shanghai http://10.100.100.11:8888/v3/resource/region"
   # alias snakemake="/data/software/miniconda/envs/snakemake/bin/snakemake"
else
    echo '流程不支持在当前节点运行, 请切换到snakemake-master节点(192.168.50.10)或server02节点(10.100.100.2)!!!'
fi

#=====================当前文件夹大小评估===================================================================================
alias du_sh="time=`date +%Y-%m-%d_%H-%M-%S`;  du -h --max-depth=1 >du.sh_${time}.log ;}"
#=====================创建超链接==========================================================================================
if [ ! -d config ] ; then cp -r xenium_pipeline/config $PWD ; fi
alias mlns=" ln -sf  xenium_pipeline/{envs,scripts,Snakefile,rules} $PWD && mkdir -p logs/{00.smt,00.smq}"
#========================查看节点资源使用情况==============================================================================
#pnodes
#====================== 查看k8s任务======================================================================================
alias qstat='func() { kubectl get pods -n $namespace;}; func'
#====================== 测试k8s任务======================================================================================
alias qrun='func() { kubectl exec -ti  $1 -n $namespace  -- /bin/sh;}; func'
#====================== 查看任务日志(目前有时间限制)： geti 2023-07-04-17-14-raw-data-upload-fgbihfb-qkb4h==================
alias qlog='func() { kubectl logs -n $namespace  `kubectl get pods -n $namespace | grep $1 | cut -f 1 -d " "`;}; func'
#====================== 删除任务:qdel 2023-07-13-14-07-stereoseq-dbeigaa-2wmrb   ========================================
alias qdel="func() { jobid=`echo $1 | sed 's/-[^-]*$//'`;echo $joibid ; kubectl delete jobs $jobid -n $namespace;}; func"
#====================== 提交测试=========================================================================================
alias smt='func() { snakemake -npr -s $1 ;}; func'
 #====================== 绘制流程图:dag Snakefile dag=====================================================================
alias dag='func() {  snakemake -npr -s $1 --dag| /data/software/conda_envs/snakemake_v6.1.1/bin/dot -Tpdf > $2.pdf;}; func'
#====================== 解锁文件 ========================================================================================
alias unlock='func() { snakemake --unlock; }; func'
#====================== 提交任务，========================================================================================
#=示例1: smq_txy pre_Snakefile
#=示例2: smq_txy Snakefile
smq()  {
    snakefile=$1;
    chmod 777 scripts/scrnasub-h*
    
    # 基于config.yaml修改cluster.yaml中的信息
    project_num=`echo $(sed -n 's/Project_Num: *"\?\([^"]*\)"\?.*/\1/p' config/config.yaml) | sed -E 's/\s+//g' | sed s/_/-/g | tr '[:upper:]' '[:lower:]'`
    # 对project_num去除其中的空格, 将其中的'_'替换为'-', 将大写字母替换为小写字母, 以便在cluster.yaml中使用(投递任务的前缀)
    # task_num=`echo $(sed -n 's/Task_Num: *"\?\([^"]*\)"\?.*/\1/p' config/config.yaml) | sed -E 's/\s+//g'`
    sed -i "s/custom_prop: .*/custom_prop: $project_num/" config/cluster.yaml 

    # 如果工作目录被锁定，则先解锁再运行
    if [ -d .snakemake/locks ] && [ "$(ls -A .snakemake/locks)" ]; then unlock; fi
    # 根据当前节点调用对应的Snakemake
    if [ "$hostName" == "snakemake-master" ]; then
            export TMP=./ && nohup /public/scrna/software/miniconda3/envs/snakemake/bin/snakemake   \
                                   --cluster "./scripts/scrnasub-hwy"  \
                                   --cluster-config config/cluster.yaml  \
                                   -j 10 \
                                   --latency-wait 600 \
                                   -s ${snakefile}  > logs/00.smq/run_${snakefile}_$(date +%Y-%m-%d-%H-%M).log
    elif [ "$hostName" == "server02" ]; then
              sed -i  's:snakemake-scrna:k8s-scrna:' config/cluster.yaml
              export TMP=./ && nohup   snakemake \
                               --cluster "./scripts/scrnasub-hhy"  \
                               --cluster-config config/cluster.yaml  \
                               -j 10 \
                               --latency-wait 600 \
                               -s ${snakefile}  > logs/00.smq/run_${snakefile}_$(date +%Y-%m-%d-%H-%M).log
    else
        echo '流程不支持在当前节点运行, 请切换到snakemake-master节点(192.168.50.10)或server02节点(10.100.100.2)!!!'
    fi
}

## 报告上传云平台所需要的环境变量配置
## !!! 建议写到自己的~/.bashrc中 !!!
# export CLOUD=True && \
# export CLOUD_URL=http://192.168.20.201:9999  && \
# export CLOUDUSER=qikai.feng@oebiotech.com  && \
# export PASSWORD=xxxx
# 访问云平台: https://cloud.oebiotech.com/spa#/
for var in  CLOUDUSER PASSWORD
do
    if [ -z "${!var}" ]; then
        echo "$var 不存在，请在个人用户.bashrc中设置, 然后重新执行该脚本"
    else
        echo "$var 已设置为'${!var}''"
    fi
done

mlns
