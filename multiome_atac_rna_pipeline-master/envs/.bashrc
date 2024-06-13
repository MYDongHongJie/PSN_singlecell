#====================================modulefile setting=================================================================
export Email_Signature="
<p> -----------------------------------------------</p>
<p>  单细胞事业部</p>
<p>  上海欧易生物医学科技有限公司 </p>
<p>  服务热线：400-640-0826 </p>
<p>  服务网址：http://www.oebiotech.com/ </p>
<p>  公司地址：上海市闵行区浦江镇联航路1188号浦江智谷25号楼 </p>
"

source  /data/software/modules/modules-v4.2.1/init/bash
export MODULEPATH=$PWD/envs/:$MODULEPATH
alias mls="module list"
alias mlo="module load"
alias mun="module unload"
alias mp="module purge"
alias mwhat="module whatis"
alias minfo="module display"
alias mhelp="echo 'mls = module list\nmlo = module load\nmun = module unload\nmp = module perge\nmwhat = module whatis\nminfo = module display'"
alias mlns="cp -r st_pipeline/config . && ln -s  st_pipeline/{envs,scripts,Snakefile,rules} ."
#====================================modulefiles load ==================================================================
module load snakemake/5.30.2
export DRMAA_LIBRARY_PATH=/software/gridengine/lib/lx-amd64/libdrmaa.so

#====================================snakemake function=================================================================
smq () {
        if [[ $# -eq 0 ]]; then
                echo 'smq <rule> [more_option] [dir]'
        else
                rule="$1"
                rule=$(readlink -e $rule)
                shift

                if [[ $# -gt 1 ]]; then
                        wd="$2"
                        shift
                else
                        wd="."
                fi

                wd=$(readlink -e $wd)
                echo runing snakemake rule: $rule under: $wd @$(date +%Y%m%d%H%M%S) with other option: $more
                if [ ! -d "$wd/logs/cluster" ];then
                    mkdir -p $wd/logs/cluster
                fi
                snakemake -s $rule --use-envmodules  -d $wd --stats $wd/snakejob.$(date +%Y%m%d%H%M%S).stats \
                          --latency-wait 1200 -k --drmaa " -V \
                                   -cwd \
                                   -e $wd/logs/cluster \
                                   -o $wd/logs/cluster \
                                   -q big  \
                                   -pe smp {resources.qsub_p} \
                                   -l h_vmem={resources.qsub_mem}G" -j 1 --drmaa-log-dir  $wd/logs/cluster/  >&2 2>> $wd/snakejob.$(date +%Y%m%d%H%M%S).log
        fi
}

#======================= snakemake test
smt () {
snakemake -npr -s $1
 }

#====================== snakemake plot
smpl () {
snakemake -npr -s $1 --dag| dot -T$2 > $3
 }

#====================== snakemake unlock
smu(){
snakemake --unlock
}

#====================== report upload
reportupload() {
        python3 -m http.server 8080
        project_report=`ls -d1 *Report*`
        project_name= $1
        tar hcvf - $project_report  | pigz -6 -p 10 -k > ${project_report}.tar.gz
        module load obsutil/5.1.11
        obsutil cp ${project_report}.tar.gz   obs://oebiotech/Analysis_Report/${project_name}/ -r -f -flat -vlength -vmd5
        echo obs://oebiotech/Analysis_Report/${project_name}/${project_report}.tar.gz
}

function du_sh { time=`date +%Y-%m-%d_%H-%M-%S`;  du -h --max-depth=1 >du.sh_${time}.log ;}
function fur_dir { time=`date +%Y-%m-%d`;
                   mkdir Further_analysis_${time} ;
                   cp envs/Further_analysis_README.txt  Further_analysis_${time}/README.txt;
                   sed -i "s/开始时间:/开始时间:$time/" Further_analysis_${time}/README.txt;}
function fur_dir_update { time=`date +%Y-%m-%d`;
                   mv $1 Further_analysis_${time};
                   sed -i "s/完成时间:/完成时间:$time/" Further_analysis_${time}/README.txt;}
