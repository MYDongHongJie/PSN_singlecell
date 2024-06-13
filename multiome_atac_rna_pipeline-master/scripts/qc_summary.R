#============================spaceranger QC results upload to SOR System ===============================================
spaceranger_sor_upload <- function(spaceranger,samplelist, project_id, project_type) {
    lapply(samplelist, function(sample) {
        ## 定义读入文件
        csv<-sprintf("%s/%s/outs/summary.csv", spaceranger, as.character(sample))
        ## 读入数据
        QC <- read_csv(csv)
        QC <- QC[ , ! colnames(QC) %in% c("Sample ID")]
        ## 单词首字母大写，以"_"连接字符串
        colnames(QC) <- str_replace_all(colnames(QC), "_", " ") %>%
          str_to_title() %>%
          str_replace_all(" ", "_") %>%
          str_replace_all("-", "_")
        ## 添加QC__type字段
        QC <- add_column(QC, QC__Type = "QC", .before = 1)
        ## 转换为json格式
        QC_json <- jsonlite::prettify(toJSON(QC))
        ## 定义上传参数
        infor_type <- "QC"
        information <- QC_json
        ##检查json格式并上传
        oe_sor_upload(project_id, project_type,infor_type,information,sample,operator="default")
    })
}
#=======================================================================================================================
check_and_email<- function(summary_csv){
    summary_stat<-read.table(summary_csv,
                             sep=",",
                             header = TRUE,
                             check.names=FALSE,
                             fill=TRUE,
                             na.strings = "" ,
                             row.names = 1,
                             stringsAsFactors=F)
    summary_stat<-as.data.frame(summary_stat)
    sample_list<-str_c(colnames(summary_stat),collapse=',')
    ###
    summary_stat_out <- summary_stat
    summary_stat_out$说明<-""
    # summary_stat_out<-summary_stat_out[c("Estimated number of cells",
    #                                  "ATAC Median high-quality fragments per cell",
    #                                  "ATAC Sequenced read pairs",
    #                                  "ATAC TSS enrichment score",
    #                                  "GEX Fraction of transcriptomic reads in cells",
    #                                  "GEX Median genes per cell",
    #                                  "GEX Sequenced read pairs"),]
    ####
    if (all(as.numeric(summary_stat["Estimated number of cells",]) > 500) & all(as.numeric(summary_stat["Estimated number of cells",]) < 10000)){
        print("各样本Estimated number of cells数据指标正常")
    }else{
        samples_high<-str_c(colnames(summary_stat)[which(as.numeric(summary_stat["Estimated number of cells",]) >= 10000)],collapse=',')
        samples_low<-str_c(colnames(summary_stat)[which(as.numeric(summary_stat["Estimated number of cells",]) < 500)],collapse=',')
        if(length(samples_high) & length(samples_low)){
            summary_stat_out["Estimated number of cells","说明"]<-paste0("样本",samples_high,"评估细胞数过多;","样本",samples_low,"评估细胞数过低;","过高或过低表明细胞核悬液定量不准确，核处理不当，背景RNA或DNA过高，调用算法有误.")
        }else if (length(samples_high)){
            summary_stat_out["Estimated number of cells","说明"]<-paste0("样本",samples_high,"评估细胞数过多;","过高或过低表明细胞核悬液定量不准确，核处理不当，背景RNA或DNA过高，调用算法有误.")
        }else{
            summary_stat_out["Estimated number of cells","说明"]<-paste0("样本",samples_low,"评估细胞数过低;","过高或过低表明细胞核悬液定量不准确，核处理不当，背景RNA或DNA过高，调用算法有误.")
        }
    }
    ###
    if (all(as.numeric(summary_stat["ATAC Median high-quality fragments per cell",]) > 100)){
        print("各样本ATAC Median high-quality fragments per cell数据指标正常")
    }else{
        samples<-str_c(colnames(summary_stat)[which(as.numeric(summary_stat["ATAC Median high-quality fragments per cell",]) <= 100)],collapse=',')
        summary_stat_out["ATAC Median high-quality fragments per cell","说明"]<-paste0("样本",samples,"高质量fragments数较低;","该值过低表明测序深度低，参考基因组有误，文库复杂度较低.")
    }
    ###
    if (all(as.numeric(summary_stat["ATAC Sequenced read pairs",])> 25000)){
        print("各样本ATAC Sequenced read pairs数据指标正常")
    }else{
        samples<-str_c(colnames(summary_stat)[which(as.numeric(summary_stat["ATAC Sequenced read pairs",]) <= 25000)],collapse=',')
        summary_stat_out["ATAC Sequenced read pairs","说明"]<-paste0("样本",samples,"read pairs数较低;","该值过低表明测序深度低，数据可信度存疑.")
    }
    ###
    if (all(as.numeric(summary_stat["ATAC TSS enrichment score",]) > 5)){
        print("各样本ATAC TSS enrichment score数据指标正常")
    }else{
        samples<-str_c(colnames(summary_stat)[which(as.numeric(summary_stat["ATAC TSS enrichment score",]) <= 5)],collapse=',')
        summary_stat_out["ATAC TSS enrichment score","说明"]<-paste0("样本",samples,"TSS 富集分数较低;","该值过低表明样本质量差，样本制备有误，具有高度可及性DNA群或参考基因组注释文件较差.")
    }
   ###
    if (all(as.numeric(summary_stat["GEX Fraction of transcriptomic reads in cells",]) > 0.6)){
        print("各样本GEX Fraction of transcriptomic reads in cells数据指标正常")
    }else{
        samples<-str_c(colnames(summary_stat)[which(as.numeric(summary_stat["GEX Fraction of transcriptomic reads in cells",]) <= 0.6)],collapse=',')
        summary_stat_out["GEX Fraction of transcriptomic reads in cells","说明"]<-paste0("样本",samples,"细胞中转录组reads的比例;","该值较低表明样本制备过程中背景RNA比例较高，也可能是算法的问题.")
    }
    ###
    if (all(as.numeric(summary_stat["GEX Median genes per cell",]) > 800)){
        print("各样本GEX Median genes per cell数据指标正常")
    }else{
        samples<-str_c(colnames(summary_stat)[which(as.numeric(summary_stat["GEX Median genes per cell",]) <= 800)],collapse=',')
        summary_stat_out["GEX Median genes per cell","说明"]<-paste0("样本",samples,"细胞中检测到的基因数过低;","该值过低可能是测序深度不足，样品质量较差，文库构建有误，参考基因组有误或基因注释较差.")
    }
    ###
    if (all(as.numeric(summary_stat["GEX Sequenced read pairs",]) > 20000)){
        print("各样本GEX Sequenced read pairs数据指标正常")
    }else{
        samples<-str_c(colnames(summary_stat)[which(as.numeric(summary_stat["GEX Sequenced read pairs",]) <= 20000)],collapse=',')
        summary_stat_out["GEX Sequenced read pairs","说明"]<-paste0("样本",samples,"read pairs数过低;","该值过低会影响聚类，差异分析以及特征关联.")
    }
    #####
    if ( all(as.numeric(summary_stat["Estimated number of cells",])> 500) &&
         all(as.numeric(summary_stat["Estimated number of cells",])< 10000) &&
         all(as.numeric(summary_stat["ATAC Median high-quality fragments per cell",]) > 100) &&
         all(as.numeric(summary_stat["ATAC Sequenced read pairs",]) > 25000) &&
         all(as.numeric(summary_stat["ATAC TSS enrichment score",]) > 5) &&
         all(as.numeric(summary_stat["GEX Fraction of transcriptomic reads in cells",]) > 0.6) &&
         all(as.numeric(summary_stat["GEX Median genes per cell",]) > 800) &&
         all(as.numeric(summary_stat["GEX Sequenced read pairs",])> 20000)){
        stat<-"正常"
        subject <- as.character(glue::glue("【质控反馈】{project_id} {customer}老师-单细胞多模态项目质控反馈"))  #主题
        body <- as.character(glue::glue('<p> 朱卉，</p>
        <p>您好，{customer}老师的10x多模态项目{project_id} {sample_number}个{species}样本（{sample_list})质控已完成。</p>
        <p>各样本cellranger运行质控结果文件详见附件，请查收。</p> ')) # 邮件正文
        merge_body <- paste0(body,c(Sys.getenv("Email_Signature")))
    }else{
        stat<-"异常"
        subject <- as.character(glue::glue("【质控异常】{project_id} {customer}老师-单细胞多模态项目质控反馈"))  #主题
        body <- as.character(glue::glue('<p> 朱卉，</p>
        <p>您好，{customer}老师的10x多模态项目{project_id} {sample_number}个{species}样本（{sample_list})质控已完成。</p>
        <p>质控结果显示部分样本的指标存在异常，详细情况见下表说明，请知悉。</p>
        <p>各样本cellranger运行质控结果文件详见附件，请查收。</p> ')) # 邮件正文
        merge_body <- paste0(body,htmlTable(summary_stat_out,rowlabel="参数"),c(Sys.getenv("Email_Signature")))
    }

    to <- email_to  # 收件人地址
    attach.files <- c(summary_csv,list.files(path = outdir, pattern ="web_summary.html",full.names =TRUE))  # 附件路径
    send.mail(from="singlecell@oebiotech.com",
                to=c(to),
                subject=subject,
                body=merge_body,
                encoding = "utf-8",
                html = TRUE,
                inline = FALSE,
                smtp = list(host.name = "smtp.oebiotech.com",
                            port = 465,
                            user.name = "singlecell@oebiotech.com" ,
                            passwd = "9LqYhXGFfkwgVhh5" ,
                            ssl = TRUE),
                authenticate = TRUE,
                send = TRUE,
                debug = FALSE,
                attach.files=attach.files)
}
#=======================================================================================================================
statisc_summary<- function(spaceranger,samplelist, outdir){
    lapply(samplelist, function(sample) {
        html<-sprintf("%s/%s/outs/web_summary.html", spaceranger, as.character(sample))
        csv<-sprintf("%s/%s/outs/summary.csv", spaceranger, as.character(sample))
        file.copy(from = html, to = paste0(outdir, "/", as.character(sample), "_", "web_summary.html"))
        file.copy(from = csv, to = paste0(outdir, "/",as.character(sample), "_", "summary.csv"))
    })
    summarylist<-list.files(path = outdir, pattern = "_summary.csv",full.names =TRUE)
    summary<- do.call(cbind, lapply(summarylist,function(i) { t(read_csv(i, col_names = TRUE)%>% mutate_if(is.numeric, round, 4))}))
    write.table(summary, file=paste0(outdir, "/", "summary.csv"), sep=",",quote = FALSE, col.names=F,na="")
    write.table(as.data.frame(t(summary)), file=paste0(spaceranger, "/", "summary.csv"), sep=",",quote = FALSE, row.names=F,col.names=T,na="")
}
#=======================================================================================================================
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("oeRtools"))
suppressPackageStartupMessages(library("mailR"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("rjson"))
suppressPackageStartupMessages(library("httr"))
suppressPackageStartupMessages(library("oeRtools"))
#suppressPackageStartupMessages(library("gdata"))
suppressPackageStartupMessages(library("htmlTable"))

# create parser object
parser <- ArgumentParser()
# default ArgumentParser will add an help option
parser$add_argument("-s", "--metadata", default="./config/samples.csv", help="metadata information csv[default \"%(default)s\"]")
parser$add_argument("-i", "--cellrnager", default="./result/cellranger", help="spaceranger output directory[default \"%(default)s\"]")
parser$add_argument("-o", "--outdir", default="./result/report", help = "output directory [default \"%(default)s\"]")
parser$add_argument("-c", "--config", default="./config/config.yaml", help="config.yaml files[default \"%(default)s\"]")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()
#=======================================================================================================================
##config.yaml parse
config<-yaml.load_file(args$config)
project_id<-config$report$Project_Num
project_type<-config$report$Project_Type
customer<-config$report$Customer
species<-config$report$Species
raw_reads<-config$params$raw_reads
email_to<-config$email$emailto

###outdir check and make
outdir <- paste0(args$outdir, "/", project_id, "-QC_multimodal_report")
if (file.exists(outdir) == "FALSE") { dir.create(outdir, recursive =TRUE)}
outdir<-suppressWarnings(normalizePath(outdir))

## sample list get
metadata<- read.table(args$metadata, sep=",",header = TRUE, fill=TRUE, na.strings = "" ,row.names=1)
samplelist<- as.vector(rownames(metadata))
sample_number<-dim(metadata)[1]
## copy and merge spaceranger statistic results
summary_anno <- statisc_summary(args$cellrnager,samplelist, outdir)
## check and send email
check_and_email(paste0(outdir, "/", "summary.csv"))

## sor uploading
#spaceranger_sor_upload(args$spacernager,samplelist, project_id, project_type)