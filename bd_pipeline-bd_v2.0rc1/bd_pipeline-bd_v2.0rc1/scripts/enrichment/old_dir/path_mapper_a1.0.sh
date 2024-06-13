#! /bin/bash
# by The Coder, 20160902
[ -d ~/.config/ ] || mkdir ~/.config/
[ -d ~/.config/ImageMagick ] || cp -a /home/fanyucai/.config/ImageMagick ~/.config/
export PATH=/home/fanyucai/lib/usr/bin/:$PATH
export LD_LIBRARY_PATH=/home/fanyucai/lib/usr/lib64/:$LD_LIBRARY_PATH
export MAGICK_CONFIGURE_PATH=$HOME/.config/ImageMagick
export MAGICK_CODER_MODULE_PATH=/home/fanyucai/lib/usr/lib64/ImageMagick-6.7.2/modules-Q16/coders

if [ ! -s anno-kegg.backgroud.xls ] || [ ! -d html_png ] || \
[ $(ls *-vs-*-diff-*.xls | wc -l) -eq 0 ]; then exit; fi


for i in $(ls *-vs-*-diff-*.xls); do
	dirn=${i%-diff-*.xls}
	mkdir $dirn

	awk 'BEGIN{FS=OFS="\t"} {p=match($5, "."); if(p>0) $5=substr($5,1,p+4)}
	NR>1{print $1, "FC: "$5, $9}' $i |
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3; next}
	(a[$1]){print $0,a[$1]}' anno-kegg.backgroud.xls - |
	sed 's/path://g; s/PATH://g' |
	awk -v dirn=$dirn 'BEGIN{FS="\t"; OFS="    "} {split($5, a, ",")}
	{for(i in a) print "&#13;",$1,$2,$3"\t"$4 > dirn"/"a[i]".xls"}'

for xls in $(ls $dirn/*.xls); do
	path=$(basename $xls | sed 's/.xls$//')
	html=html_png/$path.html
	png=html_png/$path.png
	HTML=$dirn/$path.html

	if [ -s $html ] && [ -s $png ]; then
		cp $png $dirn/
	else
		echo "WARNING: missing $path!"
		rm $xls
		continue
	fi

	awk 'NR==FNR{if($0~"^<head>") p1=FNR; if($0~"^<img src=\"") p2=NR; next}
	(FNR<=p1 || FNR>=p2){print $0}' $html $html |
	sed '/^<img id=/d' |
	sed 's#"/kegg/pathway/.*\.png#"'$path'.png#' |
	sed 's#"/tmp/mark_pathway.*\.png#"'$path'.png#' |
	sed 's#"/dbget-bin/#"http://www.genome.jp/dbget-bin/#g' |
	sed 's#"/kegg-bin/#"http://www.kegg.jp/kegg-bin/#' |
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=1; next}
	(($1~"rect" || $1~"poly") && $4!~":"){ \
	for(i in a){ p=match($4, "[^[:alnum:]]"i"[^[:alnum:]]"); if(p>0) $0=$0"\t"i}}
	{print $0}' $xls - | awk 'BEGIN{FS=OFS="\t"}
	NR==FNR{a[$2]=(a[$2] ? a[$2]""$1 : $2""$1); next}
	(NF>4){$4="title=\""a[$5]; for(i=6; i<=NF; i++){ $4=$4"&#13;"a[$i] }; \
	$4=$4"\" />"; $0=$1"\t"$2"\t"$3"\t"$4} {print $0}' $xls - > $HTML

	awk 'BEGIN{FS=OFS="\t"} ($1~"rect" || $1~"poly"){
	print $2, match($4, "Up"), match($4, "Down")}' $HTML |
	awk 'BEGIN{FS=OFS="\t"} ($2+$3)>0{sub("coords=", "", $1); \
	if($2*$3>0) $4="yellow"; if($3==0) $4="red"; \
	if($2==0) $4="green"; print $1,$4}' |
	awk 'BEGIN{FS=OFS="\t"} {split($1, a, ",")}
	length(a)==4{print "region", $2 , "46x18+"a[1]"+"a[2] }
	length(a)==8{print "draw", $2 , "line;"a[1]","a[2]";"a[3]","a[4] }' |
	while read a b c; do
	if [ "$a" == "region" ]; then
		convert $dirn/$path.png -fuzz 80% -fill $b -region $c -opaque white $dirn/$path.png
	fi
	if [ "$a" == "draw" ]; then
		draw=$(echo $c | sed 's/;/ /g')
		convert $dirn/$path.png -stroke $b -strokewidth 2 -draw "$draw"  $dirn/$path.png
	fi
	done

	rm $xls
done
done
rm *-vs-*.xls html_png anno-kegg.backgroud.xls
