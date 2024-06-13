#====================================modulefile setting=================================================================
export Email_Signature="
<p> -----------------------------------------------</p>
<p>  单细胞-空间转录组 </p>
<p>  上海欧易生物医学科技有限公司 </p>
<p>  服务热线：400-808-5350  </p>
<p>  服务网址：http://www.oebiotech.com/ </p>
<p>  公司地址：上海市闵行区浦江镇联航路1188号25号楼 </p>
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
chmod 755 st_pipeline/scripts/oest_data_prepare.py
if [ ! -d config ] ; then cp -r st_pipeline/config ./ ; fi
alias mlns="ln -sf  st_pipeline/{envs,scripts,Snakefile,pre_Snakefile,rules,run_pipeline.sh} . "
if [ ! -d logs ] ; then mkdir logs ; fi
function du_sh { time=`date +%Y-%m-%d_%H-%M-%S`;  du -h --max-depth=1 >du.sh_${time}.log ;}
#====================================modulefiles load ==================================================================
module load OESingleCell/v_3.0.0_visium_produce #snakemake/5.30.2
export DRMAA_LIBRARY_PATH=/software/gridengine/lib/lx-amd64/libdrmaa.so

#====================== snakemake unlock        
unlock () {
 snakemake --unlock
 }  

hostName=`hostname`
if [ "$hostName" == "server01" ]; then  queue="scrna" ; else queue="big"; fi
#====================================snakemake function=================================================================
smq () {
	if [[ $# -eq 0 ]]; then
		echo 'smq <rule> [more_option] [dir]'
	else
		# 如果工作目录被锁定，则先解锁再运行
		if [ -d .snakemake/locks ] && [ "$(ls -A .snakemake/locks)" ]; then unlock; fi

		rule="$1"
		rule=$(readlink -e $rule)
		shift
		more="$@"
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
		snakemake  -s $rule  \
			   --use-envmodules \
			   --cluster-config config/cluster.yaml \
			   -d $wd  \
			   --stats $wd/snakejob.$(date +%Y%m%d%H%M%S).stats \
			   --latency-wait 1200 \
			    -k  \
			    --drmaa " -V \
				     -cwd \
				     -e $wd/logs/cluster \
				     -o $wd/logs/cluster \
				     -q $queue  \
				     -pe smp {cluster.cpu} \
				     -l h_vmem={cluster.mem}G" \
			    -j $more \
			    --drmaa-log-dir  \
			    $wd/logs/cluster/  >&2 2>> $wd/logs/snakejob.$(date +%Y%m%d%H%M%S).log
	fi
}


#====================== snakemake test
smt () { 
snakemake -npr -s $1
 }  
#====================== snakemake plot                             
dag () {
 snakemake -npr -s $1 --dag| dot -T$2 > $3
 }                              
### 链接脚本
mlns 
