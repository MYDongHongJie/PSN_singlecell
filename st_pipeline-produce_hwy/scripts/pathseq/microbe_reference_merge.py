import os
import glob

# 指定基因组序列文件所在的文件夹
data_dir = './microbe_genomes'

# 获取所有基因组序列文件的路径
files = glob.glob(os.path.join(data_dir, '*_genomic.fna.gz'))

# 合并所有基因组序列文件
with open('all_genomes.fasta.gz', 'wb') as outfile:
    for f in files:
        with gzip.open(f, 'rb') as infile:
            outfile.write(infile.read())
        print(f"Processed {f}")

print("All genomes have been merged into all_genomes.fasta.gz")
