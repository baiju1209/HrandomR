date
	db=/E/BYSJ/EasyMicrobiome-main
    ea=/E/BYSJ/EasyAmplicon-master/EasyAmplicon-master
    wd=/E/BYSJ/sra-row/SRR121
    cd ${wd}/result


mkdir tree
cd tree 
#查看现有OTU数量
tail -n+2 ../otutab_rare.txt | wc -l

#筛选高丰度otu
usearch -otutab_trim ../otutab_rare.txt \
        -min_otu_freq 0.002 \
        -output otutab1.txt
#查看高丰度otu数量
tail -n+2 otutab1.txt | wc -l

#修改特征ID列名
sed -i '1 s/#OTU ID/OTUID/' otutab1.txt
#提取ID用于提取序列
cut -f 1 otutab1.txt > otutab_high.id

# 筛选高丰度菌/指定差异菌对应OTU序列
usearch -fastx_getseqs ../otus.fa -labels otutab_high.id \
    -fastaout otus1.fa
#筛选otu的物种注释
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../taxonomy.txt \
        otutab_high.id > otutab_high.tax

 #获得OTU对应组均值，用于样本热图
 #依赖之前otu_mean.R计算过按Group分组的均值
awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean.txt otutab_high.id \
        | sed 's/#OTU ID/OTUID/' > otutab_high.mean


#合并物种注释和丰度为注释文件
cut -f 2- otutab_high.mean > temp
paste otutab_high.tax temp > annotation.txt
head -n 3 annotation.txt
