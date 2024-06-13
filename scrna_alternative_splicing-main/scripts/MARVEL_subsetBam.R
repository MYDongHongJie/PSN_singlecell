suppressMessages(library(optparse))
suppressMessages(library(stringr))
suppressMessages(library(yaml))


option_list <- list(
    make_option("--marvelproject", type="character", default="result/Splicing/MARVEL.RData", help="MARVEL文件"),
    make_option("--diff_type", type="character", default=NULL, help="差异比较的类型，可以是样本之间，也可以是cluster之间，等。目前只能在样本之间计算"),
    make_option("--diff_control", type="character", default=NULL, help="差异比较的control"),
    make_option("--diff_case", type="character", default=NULL, help="差异比较的case"),
    make_option("--outpath", type="character", default="result/Splicing", help="输出文件所在位置" ),
    make_option("--configfile", type="character", default="config/cluster.yaml", help="config/cluster.yaml文件，从中提取 cores")
)


args <- parse_args(OptionParser(option_list=option_list))

marvelproject <- args$marvelproject
diff_type <- args$diff_type
diff_control <- args$diff_control
diff_case <- args$diff_case
outpath <- args$outpath
configfile <- args$configfile
config <- yaml.load_file(configfile)

#############################
#############################

source("scripts/MARVEL_functions.R")

#############################
#############################

message("加载MARVEL文件")
load(marvelproject)
message("加载完成")

#############################
#############################

diff_group <- paste0(diff_type, "_", black2underline(diff_case), "-vs-", black2underline(diff_control) )
create_folder(paste0(outpath,"/", diff_group ))
create_folder(paste0(outpath,"/", diff_group, "/sashimi" ))

# message("提取barcode")

barcodes <- subset_cell_according_diff_type(df=marvel$sample.metadata, 
                                diff_type=diff_type,
                                g1=diff_case,
                                g2=diff_control,
                                'rawbc',
                                'donor.id'
                                )


for (g in c('g1', 'g2')) {
  
    barcodes_dataframe <- data.frame('sampleid'=barcodes[[paste0('donor.id.', g)]], 
                                        # 在bam文件中，所有样本barcode的后缀都是-1，因此需要修改一下后缀
                                    'barcode'=paste0(str_split_fixed(barcodes[[paste0('rawbc.', g)]], pattern = '-', n=2)[,1], '-1'))
    samples <- c()
    for (sample in unique(marvel$sample.metadata$donor.id)) {
        
        barcodelist <- barcodes_dataframe[barcodes_dataframe$sampleid == sample,]$barcode

        if(length(barcodelist) > 0){
            # 记录哪些样本有对应的barcode
            samples <- c(samples, sample)
            # 输出barcode
            barcodefile <- paste0(outpath, '/', diff_group, "/sashimi/",sample, ".", g , ".barcodes.txt")
            write.table(barcodelist, file=barcodefile, quote = F, col.names = F, row.names = F)

            # 拆分bam
            tmpdir <- paste0(outpath,"/", diff_group, "/sashimi/tmp")
            create_folder(tmpdir)

            bam <- paste0(outpath,"/STARSolo/", sample, "/Bam/possorted_genome_bam.bam" )
            subbam <- paste0(outpath, "/", diff_group, "/sashimi/", sample, ".", g, ".bam" )
            com <- paste0("export TMPDIR=", tmpdir, "\n",
                            "subset-bam_linux --bam ", bam,
                            " --cell-barcodes ", barcodefile,
                            " --out-bam ", subbam,
                            " --cores ", config$reintegrateBam$cpu
            )

            message(paste0("正在拆分样本 ", sample))
            system(com)
            file.remove(barcodefile)
        }
    }

    # 合并同组的bam文件
    subbams <- ""
    for( sample in samples ){
    subbams <- paste0(subbams, outpath, "/", diff_group, "/sashimi/", sample, ".", g, ".bam ")
    }

    mergebam <- paste0(outpath, "/", diff_group, "/sashimi/", diff_group, ".", g , ".bam")
    com <- paste0("samtools merge ",
            " -o ", mergebam ,
            " ", subbams,
            " --threads ", config$reintegrateBam$cpu
            )

    message(paste0("正在合并 ", g))
    system(com)

    # 合并完成，删除中间bam文件
    for( sample in samples ){
        subbams <- paste0(outpath, "/", diff_group, "/sashimi/", sample, ".", g, ".bam")
        file.remove(subbams)
    }

}

unlink(tmpdir)

message("bam文件处理完毕")




