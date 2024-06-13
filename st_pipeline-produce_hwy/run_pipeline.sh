source envs/.bashrc
smt pre_Snakefile >logs/pre_smt.log
smq pre_Snakefile 
smt Snakefile >logs/smt_$(date +%Y%m%d%H%M%S).log
smq Snakefile
