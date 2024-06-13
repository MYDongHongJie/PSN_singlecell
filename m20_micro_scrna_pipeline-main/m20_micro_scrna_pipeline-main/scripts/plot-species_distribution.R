sub_read_distribution<- subparsers$add_parser("species_distribution", help = "获取kraken2分布饼图")
sub_read_distribution$add_argument("--sc_taxonomy", type = "character", help = "sc_taxonomy")
sub_read_distribution$add_argument("--prefix", type = "character", help = "sample name id ")
sub_read_distribution$add_argument("--palette",
     type = "character",
     default = "customecol2",
     choices=c("blindless","customecol2","cold","glasbey","ditto","alphabet","alphabet2","colx22","cyclic","tableau20", "paired", "purplegray12"),
     help = "the discrete color schema for each cell cluster, blindless as default. Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,tableau20:20, paired:12, purplegray12:12.[default: %(default)s]"
)
args <- commandArgs(TRUE)
if ("species_distribution" %in% args) {
  opt <- intial_setting()
  if(opt$sub_name == "species_distribution" ){
    futile.logger::flog.info("step1:导入sc_taxonomy结果")#========================================================
    sc_taxonomy <- data.table::fread(opt$sc_taxonomy,header = T)%>% dplyr::count(species)
    readr::write_tsv(sc_taxonomy,file=glue::glue("{output_dir}/all_species_distribution_for_{opt$prefix}.tsv"))
    top10_species <-sc_taxonomy %>%
                    dplyr::arrange(desc(`n`))%>%
                    dplyr::top_n(10,`n`)
    sc_taxonomy <- sc_taxonomy %>%
                   dplyr::arrange(desc(`n`))%>%
                   dplyr::rename(nCount=`n`)%>%
                   dplyr::mutate(species=ifelse(species %in% top10_species$species,species,"others")) %>%
                   dplyr::group_by(species) %>%
                   dplyr::summarise(`nCount` = sum(`nCount`)) %>%
                   dplyr::arrange(desc(`nCount`))%>%
                   dplyr::mutate(prop = round(nCount*100 / sum(nCount),2))%>%
                   dplyr::mutate(label =ifelse(prop<2,"",paste0(as.character(round(prop,2)),"%")))%>%
                   dplyr::mutate(ypos = cumsum(prop)- 0.5*prop )
    if("others" %in% sc_taxonomy$species){
      fc_level<- c(top10_species$species,"others")}else{
      fc_level<- c(top10_species$species)
    }
    sc_taxonomy$species <-factor(sc_taxonomy$species,levels=fc_level)
    sc_taxonomy <- sc_taxonomy %>% dplyr::arrange(`species`)
    readr::write_tsv(sc_taxonomy,file=glue::glue("{output_dir}/top10_species_distribution_for_{opt$prefix}.tsv"))
    futile.logger::flog.info("step2:绘制species分布饼图")#=============================================
    sum_colors <- OESingleCell::SelectColors(sc_taxonomy$species,opt$palette)
    ## png格式
    png(glue::glue("{output_dir}/piechart.top10_species_for_{opt$prefix}.png"),height = 500,width=800+max(stringr::str_length(sc_taxonomy$species))*0.1*100)
    par(mar = c(0,0,2,0))
    pie(sc_taxonomy$prop, radius=0.8,main=opt$prefix,col=sum_colors,labels=sc_taxonomy$label,cex=1.2,cex.main=2)
    legend("right", names(sum_colors), cex=1.5, bty = "n",fill = sum_colors)
    dev.off()
    pdf(glue::glue("{output_dir}/piechart.top10_species_for_{opt$prefix}.pdf"),height = 5,width=10+max(stringr::str_length(sc_taxonomy$species))*0.1)
    par(mar = c(0,0,2,0))
    pie(sc_taxonomy$prop, radius=0.8,main=opt$prefix,col=sum_colors,labels=sc_taxonomy$label,cex=1.2,cex.main=2)
    legend("right", names(sum_colors), cex=1.5, bty = "n",fill = sum_colors)
    dev.off()

    futile.logger::flog.info("step3：saving session information ") #=============================================
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
}}

