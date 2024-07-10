#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
#import sys
#reload(sys)
#sys.setdefaultencoding('utf-8')
import os
from sys import argv
from jinja2 import Template
import jinja2
import re
import glob
def get_enrichment():
    go_bar = ""
    go_dot = ""
    kegg_bar = ""
    kegg_dot = ""
    f = r"./pictures/seurat/enrichment"
    before = '<option value="./pictures/seurat/enrichment/'
    mid_go_bar = '/GO_enrichment_pvalue_barplot.png">'
    mid_go_dot = '/GO.richfactor.png">'
    mid_kegg_bar = '/KEGG_enrichment_pvalue_barplot.png">'
    mid_kegg_dot = '/KEGG.richfactor.png">'
    after = '</option>'
    num1 = 0
    heading=""
    for i, j, k in os.walk(f):
        for n in j:
            num1 += 1
            if (num1 == 1):
                heading=n
            go_bar = go_bar + before + n + mid_go_bar + n + after + '\n'
            go_dot = go_dot + before + n + mid_go_dot + n + after + '\n'
            kegg_bar = kegg_bar + before + n +mid_kegg_bar + n + after + '\n'
            kegg_dot = kegg_dot + before + n +mid_kegg_dot + n + after + '\n'
        break

    return go_bar,go_dot,kegg_bar,kegg_dot,heading

#生成top10基因UMAP图列表
def get_top10_gene_umap():
    top10_gene_umap_array = glob.glob(
        "./pictures/marker/Each_celltype_marker/*_top10_umap.png")
    top10_gene_umap = ""
    top10_gene_umap_head = ""
    print(top10_gene_umap_array)
    i = 0
    for s in top10_gene_umap_array:
        i += 1
        if (i == 1):
            top10_gene_umap_head = s
        top10_gene_umap = top10_gene_umap + '<option value=' + '"' + s + '"' + '>' + s[(s.index("/c")+1):(s.index("top")-1)] + '</option>' + "\n"
    return top10_gene_umap,top10_gene_umap_head

    # top10_gene_umap = ""
    # f = r"./05_Marker/Top10_marker_each_cluster"
    # top10_gene_umap_before = '<option value="./05_Marker/Top10_marker_each_cluster/'
    # top10_gene_umap_mid = '/top10_umap.png">'
    # top10_gene_umap_after = '</option>'
    # num1 = 0
    # top10_gene_umap_heading = ""
    # link1 = "./05_Marker/Top10_marker_each_cluster/"
    # for i, j, k in os.walk(f):
    #     for n in j:
    #         num1 += 1
    #         if (num1 == 1):
    #             top10_gene_umap_heading = n
    #             link1 = link1 + top10_gene_umap_heading + "/top10_umap.png"
    #
    #         top10_gene_umap = top10_gene_umap + top10_gene_umap_before + n + top10_gene_umap_mid + n + top10_gene_umap_after + '\n'
    #     break
    # print(top10_gene_umap)
    # print(top10_gene_umap_heading)
    # print(link1)
    # return top10_gene_umap,top10_gene_umap_heading

def get_top10_gene_violin():
    top10_gene_violin_array = glob.glob(
        "./pictures/marker/Each_celltype_marker/*_top10_vlnplot.png")
    top10_gene_violin = ""
    top10_gene_violin_head = ""
    print(top10_gene_violin_array)
    i = 0
    for s in top10_gene_violin_array:
        i += 1
        if (i == 1):
            top10_gene_violin_head = s
        top10_gene_violin = top10_gene_violin + '<option value=' + '"' + s + '"' + '>' + s[(s.index("/c")+1):(s.index("top")-1)] + '</option>' + "\n"
    return top10_gene_violin,top10_gene_violin_head

    # top10_gene_violin = ""
    # f = r"./05_Marker/Top10_marker_each_cluster"
    # top10_gene_violin_before = '<option value="./05_Marker/Top10_marker_each_cluster/'
    # top10_gene_violin_mid = '/top10_vilion.png">'
    # top10_gene_violin_after = '</option>'
    # num1 = 0
    # top10_gene_violin_heading = ""
    # link2 = "./05_Marker/Top10_marker_each_cluster/"
    # for i, j, k in os.walk(f):
    #     for n in j:
    #         num1 += 1
    #         if (num1 == 1):
    #             top10_gene_violin_heading = n
    #             link2 = link2 + top10_gene_violin_heading + "/top10_vilion.png"
    #
    #         top10_gene_violin = top10_gene_violin + top10_gene_violin_before + n + top10_gene_violin_mid + n + top10_gene_violin_after + '\n'
    #     break
    # print(top10_gene_violin)
    # print(top10_gene_violin_heading)
    # print(link2)
    # return top10_gene_violin,top10_gene_violin_heading


example = """
# 宏基因组报告

## 项目整体流程概况

### 1.1 项目信息

&[](./images/general.txt)

### 1.2 分析项目及基本要求

| 分析项目                             | 分析要求      |
| :----------------------------------- | ------------- |
| 2.1 生物信息学分析                   |               |
| 2.1.1 测序原始数据的预处理和质控     | —             |
| 2.1.2 高质量序列的筛查过滤           | —             |
| 2.1.3 基于高质量序列的物种注释       | —             |
| 2.1.4 序列的组装拼接                 | —             |
| 2.1.5 非冗余序列集的构建             | —             |
| 2.1.6 基因预测                       | —             |
| 2.1.7 蛋白功能注释                   | —             |
| 2.2 功能组成分析                     |               |
| 2.2.1 各等级功能注释相对丰度分布分析 | —             |
| 2.2.2 相对丰度差异分析               | 样本（组）≥ 2 |
| 2.2.3 共有/独有KEGG代谢通路分析      | 样本（组）≥ 2 |
| 2.2.4 KEGG代谢通路富集分析           | 样本（组）≥ 2 |
| 2.2.5 NOG/CAZy功能比较分析           | 样本（组）≥ 2 |
| 2.2.6 功能注释丰度聚类分析           | 样本（组）≥ 2 |

...

## 项目分析结果

### 2.1 生物信息学分析

#### 2.1.1 测序原始数据的预处理和质控

本项目基于**Illumina NovaSeq**高通量测序平台...

![测序读长、插入片段和接头（Adapter）的构造图](./images/icon/adapter.png)

...

#### 2.1.3 基于高质量序列的物种注释

##### 2.1.3.1 kraken2

...
"""

html = """
<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="renderer" content="webkit">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
        <title>{{ title }}</title>
        <link type="text/css" rel="stylesheet" href="static/css/bootstrap.min.css">
        <link type="text/css" rel="stylesheet" href="static/css/report.css">
        <link type="text/css" rel="stylesheet" href="static/css/bootstrap-table.css">
        <link type="text/css" rel="stylesheet" href="static/css/table.css">
        <script src="static/js/jquery-3.1.0.min.js" type="text/javascript"></script>
        <script src="static/js/bootstrap.min.js" type="text/javascript"></script>
        <script src="static/js/tables.js" type="text/javascript"></script>
        <script src="static/js/base.js" type="text/javascript"></script>
        <script src="static/js/bootstrap-table.js"></script>
        <script src="static/js/tableExport.js"></script>
        <script src="static/js/bootstrap-table-export.js"></script>
        <script src="static/js/tableList.js"></script>
        <style>.div-inline{ display:inline} </style>
        <style>
        @media print {  
            .pagination {  
             display: none;  
         }  
        }        
            select {
                width: 500px; /* 设置下拉元素的宽度为200像素 */
            }
        </style>
        <script type="text/javascript">
            function show_umap() {
                document.getElementById("heading_umap").src = document.getElementById("selection_umap").value;
            }
        </script>
        <script type="text/javascript">
            function show_violin() {
                document.getElementById("heading_violin").src = document.getElementById("selection_violin").value;
            }
        </script>

      <style>

.collapsible {
  background-color: white;
  color: #5484c7;
  cursor: pointer;
  width: 51%;
  border: none;
  text-align: left;
  outline: none;
  font-size: 15px;
}


.collapsible:after {
  color: black;
  font-weight: bold;
  float: right;
  margin-left: 5px;
  ont-size: 20px
}


button:active, .collapsible:hover {
  background-color: #E0E0E0;
}
.mycontent {
  padding: 0 18px;
  display: none;
  overflow: hidden;
  background-color: white;
}
    #pageTitle {  
        margin-top: 150px; /* 将pageTitle向下移动50像素 */  
    }  
</style>

    </head>
    <body data-spy="scroll" data-target=".scrollspy">
        <div id="slide-nav" class="d-print-none">
            <a href="http://www.personalbio.cn/" class="slide-nav-logo">
                <img src="./static/icon/logo.png" alt="" align="left">
            </a>
            <h2 class="slide-nav-title">{{ title }}</h2>
        </div>
        <div id="sidebar" class="fixed-column fixed-column-top bg-blue d-print-none bg-color">

        </div>
        <div id="frontpage" class="d-none d-print-block">
            <img src="static/icon/page_bg.jpg">
            <div id="pageTitle">
                <p>项目名称：{{ title }}</p>
                <p>委托单位：{{ partner }}</p>
                <p>制定日期：{{ date }}</p>
            </div>
        </div>
        
        <div id="floatbar" class="d-print-none">
            <div id="floatbar_goTop">
                <span class="glyphicon glyphicon-chevron-up"></span>
                <p>返回顶部</p>
            </div>
            <div class="floatbar-wrap">

                <div class="float_btn float-mouse">
                    <span class="glyphicon glyphicon-envelope"></span>
                    <p>联系方式</p>
                </div>
                <!-- <div class="float_btn"><span class="glyphicon glyphicon-qrcode"></span><p>微信订阅号</p></div> -->
                <div class="floatbar-info float-mouse">
                    <a id="fankui" href="http://c.eqxiu.com/s/U1PFE5Pk" target="_blank">
                        <span class="glyphicon glyphicon-info-sign"></span>
                        <p>信息反馈</p>
                    </a>
                </div>
                <!-- <div id="floatbar_remove"><span class="glyphicon glyphicon-remove"></span></div> -->
                <!-- <div id="floatbar_goTop"><span class="glyphicon glyphicon-chevron-up"></span><p>返回顶部</p></div> -->
            </div>
            <div class="packup packup-active">
                <i class="glyphicon glyphicon-plus"></i>
            </div>
        </div>
        <div class="float_tips tel-box">
            <div class="tel-left">
                <p class="tel-left-p">
                    <b>版权所有：</b>
                    <br>派森诺生物单细胞空转产品部
                </p>
                <p class="tel-left-p">
                    <b>邮箱：</b>
                    <br>scsupport@personalbio.cn
                </p>
            </div>
            <div class="tel-left-border">
                <div></div>
            </div>
            <img src="static/icon/2Dplot_1.jpg" alt="" class="tel-right">
        </div>
        <div class="float_tips">
            <img src="static/icon/2Dplot_1.jpg" alt="">
        </div>
        <div class="float_tips">
            <img src="static/icon/2Dplot_2.png" alt="">
        </div>
        <div class="tab-content">
            {% for page in pages %}
            <div id="part{{ page.id }}_tab_pane" class="tab-pane {% if loop.first %}active {% endif %}bg-gray d-print-block">
                <div class="catalog fixed-column bg-white d-print-none catalog-id1">
                    <!-- <a href="http://www.personalbio.cn" target="_blank" class="tab-content-a"><img src="static/icon/logo.png" link="www.personalbio.cn"></a> -->
                    <nav class="nav-pills scrollspy nav-pills-top nav-pills-one">
                    {% for part in page.parts %}
                        <a class="nav-link {% if loop.first %}{% endif %}" href="#part{{ page.id }}_ch{{ part.id }}">{{ part.id }} &nbsp;{{ part.title }}</a>
                        {% if part.subparts|length > 0 %}
                        <nav class="nav-pills nav-pills-two">
                            {% for subpart in part.subparts %}
                            <a class="nav-link" href="#part{{ page.id }}_ch{{ part.id }}_{{ subpart.id }}">{{ part.id }}.{{ subpart.id }} &nbsp;{{ subpart.title }}</a>
                            {% endfor %}
                        </nav>
                        {% endif %}
                    {% endfor %}
                    </nav>
                </div>
                <div class="content bg-white content-id{{ page.id }}">
                   
                    {% set outer_loop = loop %}
                    {% for part in page.parts %}
                    <div id="part{{ page.id }}_ch{{ part.id }}" class="mar-bot-80">
                        {% if not outer_loop.first or loop.index != 2 %}
                        <div class="col-12 border-bottom mb-5 d-none d-print-block" style="page-break-before: always;">
                            <img src="static/icon/logo.png" class="col-2 offset-10 mb-2"/>
                        </div>
                        {% endif %}
                        <h2 class="chapter-title base-item1">{{ part.id }} &nbsp;{{ part.title }}</h2>
                        {% if part.subparts|length > 0 %}
                        {% for subpart in part.subparts %}
                        <div id="part{{ page.id }}_ch{{ part.id }}_{{ subpart.id }}" class="mar-bot-80">
                            <h3 class="section-title base-item10">{{ part.id }}.{{ subpart.id }} &nbsp;{{ subpart.title }}</h3>
                            {{ subpart.text }}
                        </div>
                        {% endfor %}
                        {% else %}
                        {{ part.text }}
                        {% endif %}
                    </div>
                    {% endfor %}
                </div>
            </div>
            {% endfor %}
        </div>
        <div class="carousel-img-zoom">
            <img src="./result/img/1.png" alt="">
        </div>
        <script src="static/js/main.js"></script>

<script>
var coll = document.getElementsByClassName("collapsible");
var i;

for (i = 0; i < coll.length; i++) {
  coll[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    if (content.style.display === "block") {
      content.style.display = "none";
    } else {
      content.style.display = "block";
    }
  });
}
</script>
<script>
      $('.selected_img img').click(function (){
      $(this).attr("width", $(this).attr("width") == "100%" ? "100%" : "100%")
      })
      $('.img-select').change(function (){
      $(this).next(".selected_img").children("img").attr("src", $(this).val())
      })
      $('.iframe-select').change(function(){
      $(this).next("iframe").attr("src", $(this).val())
      })




</script>
    </body>

</html>
"""

para = """
<p class="paragraph">{text}</p>
"""

bold = """<strong>{text}</strong>"""

link = """<a href="{url}" target="_blank">{name}</a>"""

img = """
                            <div component_type="img">
                                <div class="row">
                                    <div class="col-xl-12 col-lg-12 col-12 offset-xl-0 offset-lg-0 offset-0">
                                        <img class="img-fluid" alt="" src="{url}" width="60%">
                                        <b class="imgname col-12">{title}</b>
                                        <p class="tablenote col-12">
                                        {text}
                                        </p>
                                    </div>
                                </div>
                            </div>
"""

table = """
                            <div class="row" component_type="table">
                                <div class="col-xl-10 col-lg-10 col-10 offset-xl-1 offset-lg-1 offset-1">
                                    <table class="onehead_table three-line-table table-sm" id="">
                                        <tbody>
                                        {% for item in table %}
                                        <tr>
                                            {% for it in item %}
                                            <th>{{ it }}</th>
                                            {% endfor %}
                                        </tr>
                                        {% endfor %}
                                        </tbody>
                                    </table>
                                    </div>
                            </div>
"""

table_link = """            
                            <div class="row" component_type="table">
                                <div class="col-xl-10 col-lg-10 col-12 offset-xl-1 offset-lg-1">
                                    <b class="tablename col-12">{{ title }}</b>
                                    <div id="reportTableDiv">
                                    <table id="reportTable{{ num }}"></table>
                                    </div>
                                    <script type="text/javascript">
                                        $('#reportTable{{ num }}').bootstrapTable({
                                            height: 300,
                                            pagination: true,
                                            pageSize: 10,
                                            pageNumber: 1,
                                           
                                            exportDataType: 'basic',
                                            exportTypes: ['csv', 'xls', 'doc', 'txt'],

                                            tableTitle: '{{ table.title }}',
                                            optionalParams: {
                                                fileName: '{{ table.title }}',
                                            },
                                            titleContent: {{ table.content }},
                                            columns: [
                                            {% for col in table.cols %}
                                            {field:"{{ col }}",title:"{{ col }}",align:"center",valign:"middle",sortable:"true"},
                                            {% endfor %}
                                            ],
                                            data: [
                                            {% for data in table.datas %}
                                            {{ data.__str__() }},
                                            {% endfor %}
                                            ],
                                            ellipsis: false,
                                            formatLoadingMessage: function () {
                                                return "请稍等，正在加载中...";
                                            },
                                            onSearch: function (text) {
                                            },
                                            onPageChange: function (size, number) {

                                            },
                                            formatNoMatches: function() {
                                                return '暂无相关内容！';
                                            }
                                        });
                                        $(window).resize(function () {
                                            $('#reportTable{{ num }}').bootstrapTable('resetView');
                                        });
                                    </script>
                                </div>
                            </div>
"""

def parse_table(item):
    with open(item) as t:
        content = t.read().strip().split("\n")

    content = [ line.split("\t") for line in content ]
    data = [ dict(zip(content[0], line)) for line in content[1:] ]
    table = {
        "title" : "",
        "cols" : content[0],
        "content" : content[0].__str__(),
        "datas" : data
    }
    return table

def insert_table(tbl):
    global table
    t = Template(table)
    return t.render(table=tbl)

def cons_page(page):
    pages = []
    for p in page:
        if p["level"] == "h2":
            pages.append(p)
        elif p["level"] == "h3":
            pages[-1]["parts"].append(p)
        elif p["level"] == "h4":
            pages[-1]["parts"][-1]["subparts"].append(p)

    return pages
    
    
def parse(md_f):
    html_text = ""
    global html
    page = []
    table = []
    table_cnt = 0
    id = 0
    num = 0
    subnum = 0
    text = ""
    with open(md_f) as md:
        for line in md:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith("|"):
                if not "----" in line:
                    row = line.split("|")
                    row = [ item.strip() for item in row ]
                    table.append(row[1:-1])
            else:
                if table:
                    text += insert_table(table)
                    table = []
                
            row = line.split(" ")
            title = " ".join(row[1:])

            if row[0] == "#":
                main_title = title
            elif row[0] == "##":
                if page:
                    page[-1]["text"] = text
                    text = ""
                id += 1
                num = 0
                subnum = 0
                page.append({
                    "level" : "h2",
                    "id" : id,
                    "title" : title,
                    "parts" : []
                })
                
            elif row[0] == "###":
                if page:
                    page[-1]["text"] = text
                    text = ""
                num += 1
                subnum = 0
                page.append({
                    "level" : "h3",
                    "id" : num,
                    "title" : title,
                    "text" : "",
                    "subparts" : []
                })
            elif row[0] == "####":
                if page:
                    page[-1]["text"] = text
                    text = ""
                subnum += 1
                page.append({
                    "level" : "h4",
                    "id" : subnum,
                    "title" : title,
                    "text" : "",
                })
            elif row[0] == "#####":
                text += '<h4 class="chapter-title base-item8">{}</h4>'.format(title.replace(' ', ' &nbsp;'))
            elif line.startswith("!"):
                name = line[2:line.find("]")]
                url = line[line.find("(")+1:line.find(")")]
                global img
                text += img.format(url=url, title=name, text="")
            elif line.startswith("&"):
                name = line[2:line.find("]")]
                url = line[line.find("(")+1:line.find(")")]
                table_cnt += 1
                global table_link
                t = Template(table_link)
                text += t.render(table=parse_table(url), title=name, num=table_cnt)
            elif line.startswith("|"):
                pass
            elif line.startswith("<"):
                text += line
            elif line.endswith(">"):
                text += line
            else:
                # bold
                strong = re.compile('\*\*([^\*]*)\*\*')
                image = re.compile('\[([^\[]*)\]\(([^\(]*)\)')
                italic = re.compile('\*([^\*]*)\*')
                line = re.sub(strong, lambda x: bold.format(text=x.groups()[0]), line)
                line = re.sub(image, lambda x: link.format(name=x.groups()[0], url=x.groups()[1]), line)
                line = re.sub(italic, lambda x: '<i>{}</i>'.format(x.groups()[0]), line)
                
                html_ordered_list = re.compile(r'(<ol>(?:\s*<li>.*?</li>)+\s*</ol>)', re.DOTALL)
                line = html_ordered_list.sub(lambda match: "<ol>\n" + 
                                            "\n".join("<li>" + x.strip()[len(re.match(r'<li>\s*', x).group(0)):-len('</li>')] + "</li>" 
                                            if re.match(r'<li>\s*', x)
                                            else "</ol>\n<ol>\n<li>" + x.strip() + "</li>"
                                            for x in match.group(1).split("\n")) + 
                                            "</ol>", line)
                html_ordered_list = re.compile(r'(<ol>(?:\s*<li>.*?</li>)+\s*</ol>)', re.DOTALL)

# Convert HTML ordered list to HTML
                line = html_ordered_list.sub(lambda match: "<ol>\n" +  
                                            "\n".join("<li>" + x.strip()[len(re.match(r'<li>\s*', x).group(0)):-len('</li>')] + "</li>" 
                                            if re.match(r'<li>\s*', x)
                                            else "</ol>\n<ol>\n<li>" + x.strip() + "</li>"
                                            for x in match.group(1).split("\n")) + 
                                            "</ol>", line)

                text += para.format(text=line)
                
    page[-1]["text"] = text
    pages = cons_page(page)
    global html
    t = Template(html)
    top10_gene_umap,top10_gene_umap_head = get_top10_gene_umap()
    top10_gene_violin,top10_gene_violin_heading=get_top10_gene_violin()
    go_bar, go_dot, kegg_bar, kegg_dot, heading=get_enrichment()
    template = Template(t.render(pages=pages,
                    title=main_title,partner = argv[2],date = datetime.date.today()))
    return template.render(
        top10_gene_umap=top10_gene_umap,
        top10_gene_umap_head=top10_gene_umap_head,
        top10_gene_violin=top10_gene_violin,
        top10_gene_violin_heading=top10_gene_violin_heading,
        go_bar=go_bar,
        go_dot=go_dot,
        kegg_bar=kegg_bar,
        kegg_dot=kegg_dot,
        heading=heading
                           )


def main():
    #print(parse(argv[1]))
    with open("report.html", 'w+', encoding='utf-8') as f:
        f.write(parse(argv[1]))
        f.close()

if __name__ == "__main__":
    main()
