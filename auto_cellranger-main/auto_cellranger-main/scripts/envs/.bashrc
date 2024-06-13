#====================================modulefile setting=================================================================
export Email_Signature="
<p> -----------------------------------------------</p>
<p>  单细胞事业部-转录组 </p>
<p>  上海欧易生物医学科技有限公司 </p>
<p>  服务热线：021-34781616 </p>
<p>  服务网址：http://www.oebiotech.com/ </p>
<p>  公司地址：上海市闵行区浦江镇联航路1188号25号楼 </p>
"

source  /data/software/modules/modules-v4.2.1/init/bash
export PATH=/data/software/obsutil/5.2.12/:$PATH
export MODULEPATH=$PWD/envs/:$MODULEPATH
alias mls="module list"
alias mlo="module load"
alias mun="module unload"
alias mp="module purge"
alias mwhat="module whatis"
alias minfo="module display"
alias mhelp="echo 'mls = module list\nmlo = module load\nmun = module unload\nmp = module perge\nmwhat = module whatis\nminfo = module display'"
alias mlns="ln -s scrna_pipeline/{envs,scripts,Snakefile,rules} . ; test -d config/ && echo config folder exists. || cp -r scrna_pipeline/config ."
# alias smpr="/home/xfyang/miniconda3/4.6.14/bin/python scripts/data_prepare.py"

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

                if [[ $# -gt 0 ]]; then
                        wd="$1"
                        shift
                        more="$@"
                else
                        wd="."
                        more="$@"
                fi
                wd=$(readlink -e $wd)
                
                echo runing snakemake rule: $rule under: $wd @$(date +%Y%m%d%H%M%S) with other option: $more
                if [ ! -d "$wd/logs/cluster" ];then
                    mkdir -p $wd/logs/cluster
                fi
                snakemake -s $rule  --use-envmodules -d $wd --stats $wd/snakejob.$(date +%Y%m%d%H%M%S).stats \
                          --latency-wait 1200 -k  --cores all \
                          --drmaa " -V \
                                   -cwd \
                                   -e $wd/logs/cluster \
                                   -o $wd/logs/cluster \
                                   -q big  \
                                   -pe smp {resources.qsub_p} \
                                   -l h_vmem={resources.qsub_mem}G" \
                          -j $more --drmaa-log-dir  $wd/logs/cluster/  >&2 2>> $wd/snakejob.$(date +%Y%m%d%H%M%S).log
                          # --cluster-config cluster.yaml \
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
if [[ $# -eq 0 ]]; then snakemake --unlock
else snakemake --unlock -s $1
fi
}



