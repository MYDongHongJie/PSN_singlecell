
#============================spaceranger QC results upload to SOR System ===============================================
# spaceranger_sor_upload <- function(spaceranger,samplelist, project_id, project_type) {
#     lapply(samplelist, function(sample) {
#         ## 定义读入文件
#         csv<-sprintf("%s/%s/outs/metrics_summary.csv", spaceranger, as.character(sample))
#         ## 读入数据
#         QC <- read_csv(csv)
#         QC <- QC[ , ! colnames(QC) %in% "Sample ID"]
#         ## 单词首字母大写，以"_"连接字符串
#         colnames(QC) <- str_replace_all(colnames(QC), "_", " ") %>%
#           str_to_title() %>%
#           str_replace_all(" ", "_") %>%
#           str_replace_all("-", "_")
#         ## 添加QC__type字段
#         QC <- add_column(QC, QC__Type = "QC", .before = 1)
#         ## 转换为json格式
#         QC_json <- jsonlite::prettify(toJSON(QC))
#         ## 定义上传参数
#         infor_type <- "QC"
#         information <- QC_json
#         ##检查json格式并上传
#         oe_sor_upload(project_id, project_type,infor_type,information,sample,operator="default")
#     })
# }

#=======================================================================================================================
# check_and_email<- function(summary_csv){
#     summary_stat<-read.table(summary_csv,
#                                            sep=",",
#                                            header = TRUE,
#                                            check.names=FALSE,
#                                            fill=TRUE,
#                                            na.strings = "" ,
#                                            row.names=1)
#     summary_stat<-as.data.frame(summary_stat)
#     if(library_type != "fresh"){Reads_mappeds <- "Reads Mapped to Probe Set"}else{ Reads_mappeds <- "Reads Mapped to Genome"}
#     sample_list<-str_c(colnames(summary_stat),collapse=',')
#     ###判断存在bug，需要调整修正
#     if ( all(summary_stat["Number of Reads",] >= as.integer(raw_reads)*1e6*0.9) &&
#          #all(summary_stat["Median Genes per Spot",] >= 3000) &&
#          #all(summary_stat["Median UMI Counts per Spot",] >= 10000) &&
#          all(summary_stat["Valid Barcodes",] >= 0.75) &&
#          all(summary_stat["Valid UMIs",] >= 0.75) &&
#          all(summary_stat[Reads_mappeds,] >= ifelse(library_type != "fresh",0.50,0.85)) &&
#          all(summary_stat["Fraction Reads in Spots Under Tissue",] >= 0.50)){
#         stat<-"正常"
#         body <-as.character(glue::glue('<p> 朱卉，</p>
#         <p>您好，{customer}老师的空间转录组项目{project_id} {sample_number}个{species}样本({sample_list})质控已完成。</p>
#         <p>质控结果显示{sample_list}等{sample_number}个样本的基因，UMI等各项指标处于正常范围内，请知悉。</p>
#         <p>各样本质控概况如下，spaceranger运行结果文件见附件，请查收。</p> '))
#     }else{
#         stat<-"异常"
#         body <- as.character(glue::glue('<p> 朱卉，</p>
#         <p>您好，{customer}老师的空间转录组项目{project_id} {sample_number}个{species}样本（{sample_list})质控已完成。</p>
#         <p>质控结果显示{sample_list}等{sample_number}个样本的部分指标存在异常，详细情况见下表说明，请知悉。</p>
#         <p>各样本spaceranger运行质控结果文件详见附件，请查收。</p> ')) # 邮件正文
#     }
#    # write.table(body,file=paste0(summary_csv,".email"),col.names =FALSE,row.names = FALSE,quote = FALSE)
#     #summary_stat_out<-as.data.frame(read.table(summary_csv, sep=",",header = TRUE,fill=TRUE, na.strings = ""))
#     summary_stat_out<-summary_stat
#     summary_stat_out$说明<-""
#    ###
#     if ( all(summary_stat["Number of Reads",] >= as.integer(raw_reads)*1e6*0.9)){
#         summary_stat_out["Number of Reads","说明"]<-as.character(glue::glue('数据量足够(要求数据量:{raw_reads}M)'))
#         max_min<-round(max(summary_stat["Number of Reads",])/min(summary_stat["Number of Reads",]),2)
#         if( max_min > 1.2){
#         summary_stat_out["Number of Reads","说明"]<-as.character(
#           glue::glue('数据量足够(要求数据量:{raw_reads}M),
#           数据量最高与最低样本相差{max_min}倍'))
#         }
#     }else{
#         samples<-str_c(colnames(summary_stat)[which(summary_stat["Number of Reads",]<as.integer(raw_reads)*1e6*0.9)],collapse=',')
#         summary_stat_out["Number of Reads","说明"]<-paste0("样本",samples,"数据量不足")
#     }
#    ####
#     # if ( all(summary_stat["Median Genes per Spot",] >= 3000)){
#     #     summary_stat_out["Median Genes per Spot","说明"]<-"正常"
#     # }else{
#     #     samples<-str_c(colnames(summary_stat)[which(summary_stat["Median Genes per Spot",]<3000)],collapse=',')
#     #     summary_stat_out["Median Genes per Spot","说明"]<-paste0("样本",samples,"低于现有质控标准3000")
#     # }
#     ####
#     # if ( all(summary_stat["Median UMI Counts per Spot",] >= 10000 )){
#     #     summary_stat_out["Median UMI Counts per Spot","说明"]<-"正常"
#     # }else{
#     #     samples<-str_c(colnames(summary_stat)[which(summary_stat["Median UMI Counts per Spot",]< 10000)],collapse=',')
#     #     summary_stat_out["Median UMI Counts per Spot","说明"]<-paste0("样本",samples,"低于现有质控标准10000")
#     # }
#     if ( all(summary_stat["Valid Barcodes",] >= 0.75)){
#         summary_stat_out["Valid Barcodes","说明"]<-"正常"
#     }else{
#         samples<-str_c(colnames(summary_stat)[which(summary_stat["Valid Barcodes",]<0.75)],collapse=',')
#         summary_stat_out["Valid Barcodes","说明"]<-paste0("样本",samples,"低于现有质控标准0.75,提示可能与测序质量有关")
#     }
#     ###
#     if ( all(summary_stat["Valid UMIs",] >= 0.75)){
#         summary_stat_out["Valid UMIs","说明"]<-"正常"
#     }else{
#         samples<-str_c(colnames(summary_stat)[which(summary_stat["Valid UMIs",]<0.75)],collapse=',')
#         summary_stat_out["Valid UMIs","说明"]<-paste0("样本",samples,"低于现有质控标准0.75，提示可能与文库或者测序质量有关")
#     }
#    ####
#     if ( all(summary_stat[Reads_mappeds,] >= ifelse(library_type != "fresh",0.50,0.85))){
#         summary_stat_out[Reads_mappeds,"说明"]<-"正常"
#     }else{
#         samples<-str_c(colnames(summary_stat)[which(summary_stat[Reads_mappeds,]<ifelse(library_type != "fresh",0.50,0.85))],collapse=',')
#        summary_stat_out[Reads_mappeds,"说明"]<-glue::glue("样本{samples}低于现有质控标准{ifelse(library_type != 'fresh',0.50,0.85)}")
#     }
#    ####
#
#     if ( all(summary_stat["Fraction Reads in Spots Under Tissue",] >= 0.50)){
#         summary_stat_out["Fraction Reads in Spots Under Tissue","说明"]<-"正常"
#     }else{
#         samples<-str_c(colnames(summary_stat)[which(summary_stat["Fraction Reads in Spots Under Tissue",]<0.50)],collapse=',')
#         summary_stat_out["Fraction Reads in Spots Under Tissue","说明"]<-paste0("样本",samples,"低于现有质控标准0.50")
#     }
#     merge_body <- paste0(body, htmlTable(summary_stat_out,rowlabel="SampleID"), Sys.getenv("Email_Signature"))
#     # write.fwf(summary_stat_out,file=paste0(summary_csv,".email") ,append =TRUE,rownames=FALSE,colnames = FALSE,scientific=FALSE)
#     #write.table(formattable(summary_stat),file=paste0(summary_csv,".email"),append =TRUE)
#     write.table(merge_body,file=paste0(summary_csv,".email.html"))
#     #from <- Sys.getenv("Email_From")  # 发件人邮箱
#     to <- email_to  # 收件人地址
#     #password<-Sys.getenv("Email_PassWord")  # 邮箱授权码
#     subject <- as.character(glue::glue("【质控反馈】{project_id} {customer}老师-Visium项目质控反馈"))  #主题
#     attach.files <- c(summary_csv,list.files(path = outdir, pattern ="web_summary.html",full.names =TRUE))  # 附件路径
#     tryCatch(
#         error = function(cnd) {
#             cat("Failed to send email. Possibly because 1. a connection issue occur  OR 2.the email attachments are too large.
# Please check 'summary.csv.email.html' and send the email manually.
# ")
#         },
#         send.mail(from="singlecell@oebiotech.com",
#                 to= to,
#                 subject=subject,
#                 body=merge_body,
#                 encoding = "utf-8",
#                 html = TRUE,
#                 inline = FALSE,
#                 smtp = list(host.name = "smtp.oebiotech.com",  # 此处以QQ邮箱为例
#                             port = 465,  # QQ邮箱发送服务器端口
#                             user.name = "singlecell@oebiotech.com" ,
#                             passwd = "9LqYhXGFfkwgVhh5" ,
#                             ssl = TRUE),
#                 authenticate = TRUE,
#                 send = TRUE,
#                 debug = FALSE,
#                 attach.files=attach.files)
#     )
# }
#=======================================================================================================================
statisc_summary<- function(starsolo, samplelist, outdir){
    # lapply(samplelist, function(sample) {
    #     html<-sprintf("%s/%s/outs/web_summary.html", spaceranger, as.character(sample))
    #     csv<-sprintf("%s/%s/outs/metrics_summary.csv", spaceranger, as.character(sample))
    #     file.copy(from = html, to = paste0(outdir, "/", as.character(sample), "_", "web_summary.html"))
    #     file.copy(from = csv, to = paste0(outdir, "/",as.character(sample), "_", "metrics_summary.csv"))sctool
    # })
    summarylist<-list.files(path = starsolo, pattern = "Summary.csv",all.files=T,full.names =TRUE,recursive =TRUE)
    summarylist<-summarylist[str_detect(summarylist,pattern="GeneFull")]
    summary<- do.call(cbind, lapply(summarylist,function(i) { t(read_csv(i, col_names = FALSE)%>% mutate_if(is.numeric, round, 4))}))
    rownames(summary)<-c("sampleid",samplelist)
    print(samplelist)
    write.table(summary, file=paste0(outdir, "/", "STARsolo.statistic.tsv"), sep="\t",quote = FALSE, col.names=F)
}
#=======================================================================================================================
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("argparse"))
# suppressPackageStartupMessages(library("oeRtools"))
# suppressPackageStartupMessages(library("mailR"))
# suppressPackageStartupMessages(library("yaml"))
# suppressPackageStartupMessages(library("rjson"))
# suppressPackageStartupMessages(library("httr"))
# suppressPackageStartupMessages(library("oeRtools"))
# suppressPackageStartupMessages(library("gdata"))
# suppressPackageStartupMessages(library("htmlTable"))

docstring <- " example1:\\n\\n\\
sctool    -o results   summary  -s   config/samples.csv  -i result/1.STARsolo/ "
sub_email <- subparsers$add_parser(
    "summary",
    description = docstring,
    formatter_class = "argparse.RawTextHelpFormatter",
    # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
    argument_default = "True",
    help = "Using RCTD to refer Spaital transcriptome cell type"
)
sub_email$add_argument(
"-s", "--metadata", default="config/samples.csv", help="metadata information csv[default \"%(default)s\"]")
sub_email$add_argument(
"-i", "--starsolo", default="result/1.STARsolo/", help="spaceranger output directory[default \"%(default)s\"]")
# sub_email$add_argument(
# "-c", "--config", default="config/config.yaml", help="config.yaml files[default \"%(default)s\"]")

# ================ Subcmd: st_deconv, deconvolute the cell type composition for spatial transcriptomics or bulk transcriptomics using scRNA-seq data ========
args <- commandArgs(TRUE)
if ("summary" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "summary") {
        #=======================================================================================================================
        ###step1 outdir check and make
        if (file.exists(output_dir) == "FALSE") { dir.create(output_dir, recursive =TRUE)}
        outdir<-suppressWarnings(normalizePath(output_dir))
        ## sample list get
        metadata<- read.table(opt$metadata, sep=",",header = TRUE, fill=TRUE, na.strings = "" ,row.names=1)
        samplelist<- as.vector(rownames(metadata))
        print(samplelist)
        print(outdir)
        ## copy and merge spaceranger statistic results
        statisc_summary(opt$starsolo, samplelist, output_dir)
        #=======================================================================================================================
        ###step3  check and send email
        #check_and_email(paste0(outdir, "/", "summary.csv"))
        ## sor uploading
        #spaceranger_sor_upload(opt$spacernager,samplelist, project_id, project_type)
        # ==============================================================================================================
        ## save session information
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
 }
}