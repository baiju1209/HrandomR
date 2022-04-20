  date
	db=/E/BYSJ/EasyMicrobiome-main
    ea=/E/BYSJ/EasyAmplicon-master/EasyAmplicon-master
    wd=/E/BYSJ/sra-row/SRR552
    cd ${wd}

    mkdir -p temp
    mkdir -p result 
	mkdir -p result/alpha
    cp metadata.txt result
	
#	for i in `tail -n+2 metadata.txt | cut -f 1`;do
#	vsearch -fastq_mergepairs ${i}_1.fastq.gz \-reverse ${i}_2.fastq.gz \-fastqout temp/${i}.merged.fq -relabel ${i}.
#	done 

#单端的这个
#	gunzip *.gz
#	time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
 #      usearch -fastx_relabel ${i}.fastq -fastqout temp/${i}.merged.fq -prefix ${i}.
  #  done &




	 #合并所有样品至同一文件
	cat temp/*.merged.fq > temp/all.fq
	 #查看文件大小223M，软件不同版本结果略有差异
	ls -lsh temp/all.fq
	 #查看序列名，“.”之前是否为样本名，样本名绝不允许有.
	head -n 6 temp/all.fq|cut -c1-60
	 ## 3. 切除引物与质控 Cut primers and quality filter
	 
	     # 左端10bp标签+19bp上游引物V5共为29，右端V7为18bp下游引物
	     # Cut barcode 10bp + V5 19bp in left and V7 18bp in right
	     # 务必清楚实验设计和引物长度，引物已经去除可填0，27万条序列14s
	time vsearch --fastx_filter temp/all.fq \
	  --fastq_stripleft 58 --fastq_stripright 45 \
	  --fastq_maxee_rate 0.01 \
	  --fastaout temp/filtered.fa
	     # 查看文件了解fa文件格式
	head temp/filtered.fa
		 
	 #4.1 序列去冗余 Dereplicate
	 # 并添加miniuniqusize最小为8或1/1M，去除低丰度噪音并增加计算速度
	 # -sizeout输出丰度, --relabel必须加序列前缀更规范, 1s
	vsearch --derep_fulllength temp/filtered.fa \
	  --minuniquesize 1 --sizeout --relabel Uni_ \
	  --output temp/uniques.fa 
	 #高丰度非冗余序列非常小(此处852KB)，名称后有size和频率
	ls -lsh temp/uniques.fa
	head -n 2 temp/uniques.fa
	 
	#  4.2 聚类OTU/去噪ASV Cluster or denoise
	 #方法2. ASV去噪 Denoise: predict biological sequences and filter chimeras
	 #6s, 1530 good, 41 chimeras, 序列百万条可能需要几天/几周
	usearch -unoise3 temp/uniques.fa -minsize 10 \
	  -zotus temp/zotus.fa
	 #修改序列名：Zotu为改为ASV方便识别
	sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
	head -n 2 temp/otus.fa
	 
	 
	mkdir -p result/raw
	 
	 # 方法1. vsearch+rdp去嵌合(快但容易假阴性)
	 # 可自行下载silva并解压(替换rdp_16s_v18.fa为silva_16s_v123.fa)，极慢但理论上更好
	vsearch --uchime_ref temp/otus.fa \
	  -db ${db}/usearch/rdp_16s_v18.fa \
	  --nonchimeras result/raw/otus.fa
	 # RDP: 7s, 143 (9.3%) chimeras; SILVA：9m, 151 (8.4%) chimeras
	 # Win vsearch结果添加了windows换行符^M需删除，mac不要执行此命令
	sed -i 's/\r//g' result/raw/otus.fa
	 
	
	  
	  # 方法2. vsearch生成特征表
	  # id(1)：100%相似度比对49.45%序列，1m50s
	  # id(0.97)：97%相似度比对83.66%序列，1m10s(更高数据使用率，更快)
	  
	  #我的小破电脑用了十来分钟才搞定这一个步骤
	time vsearch --usearch_global temp/filtered.fa \
	  --db result/raw/otus.fa \
	  --id 0.97 --threads 4 \
	  --otutabout result/raw/otutab.txt 
		
		
	  #224236 of 268019 (83.66%)可比对
	  # vsearch结果windows用户删除换行符^M校正为标准Linux格式
	sed -i 's/\r//' result/raw/otutab.txt
	head -n6 result/raw/otutab.txt | cut -f 1-6 |cat -A
	  # csvtk统计表行列
	 # csvtk -t stat result/raw/otutab.txt
	  
	  
	  
	  # 物种注释-去除质体和非细菌/古菌并统计比例(可选)
	  # RDP物种注释(rdp_16s_v18)更快，但缺少完整真核来源数据,可能不完整，耗时15s;
	  # SILVA数据库(silva_16s_v123.fa)更好注释真核、质体序列，极慢耗时3h起
	  # 置信阈值通常0.6/0.8，vserch最低0.1/usearch可选0输出最相似物种注释用于观察潜在分类
	vsearch --sintax result/raw/otus.fa \
	  --db ${db}/usearch/rdp_16s_v18.fa \
	  --sintax_cutoff 0.1 \
	  --tabbedout result/raw/otus.sintax 
	head result/raw/otus.sintax | cat -A
	sed -i 's/\r//' result/raw/otus.sintax
	  
	  
	  
	 
	  # 方法1. 原始特征表行数
	wc -l result/raw/otutab.txt
	    #R脚本选择细菌古菌(真核)、去除叶绿体、线粒体并统计比例；输出筛选并排序的OTU表
	    #输入为OTU表result/raw/otutab.txt和物种注释result/raw/otus.sintax
	    #输出筛选并排序的特征表result/otutab.txt和
	    #统计污染比例文件result/raw/otutab_nonBac.txt和过滤细节otus.sintax.discard
	    #真菌ITS数据，请改用otutab_filter_nonFungi.R脚本，只筛选真菌
	    # Rscript ${db}/script/otutab_filter_nonBac.R -h # 显示参数说明
	Rscript ${db}/script/otutab_filter_nonBac.R \
	  --input result/raw/otutab.txt \
	  --taxonomy result/raw/otus.sintax \
	  --output result/otutab.txt\
	  --stat result/raw/otutab_nonBac.stat \
	  --discard result/raw/otus.sintax.discard
		  # 筛选后特征表行数
	wc -l result/otutab.txt
		  #过滤特征表对应序列
	cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
	usearch -fastx_getseqs result/raw/otus.fa \
	  -labels result/otutab.id -fastaout result/otus.fa
		  #过滤特征表对应序列注释
	awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
		result/raw/otus.sintax result/otutab.id \
		> result/otus.sintax
			  
			  #不统计了，说是32位不支持
			  
		### 5.3 等量抽样标准化
		    # Normlize by subsample
		
		    #使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
		    #可指定输入文件、抽样量和随机数，输出抽平表result/otutab_rare.txt和多样性alpha/vegan.txt
		   
	Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
	  --depth 10000 --seed 1 \
	  --normalize result/otutab_rare.txt \
	  --output result/alpha/vegan.txt
	usearch -otutab_stats result/otutab_rare.txt \
	  -output result/otutab_rare.stat
			  
	cat result/otutab_rare.stat
		## 6. α多样性 alpha diversity
		
		### 6.1. 计算α多样性 calculate alpha diversity
		
		    # 使用USEARCH计算14种alpha多样性指数(Chao1有错勿用)
		    #details in http://www.drive5.com/usearch/manual/alpha_metrics.html
	usearch -alpha_div result/otutab_rare.txt \
	  -output result/alpha/alpha.txt
			  
			  
		### 6.2. 计算稀释丰富度 calculate rarefaction richness
		
		    #稀释曲线：取1%-100%的序列中OTUs数量，每次无放回抽样
		    #Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
		  #  usearch -alpha_div_rare result/otutab_rare.txt \
		  #    -output result/alpha/alpha_rare.txt \
		  #    -method without_replacement
		    #预览结果
		 #   head -n2 result/alpha/alpha_rare.txt
		    #样本测序量低出现非数值"-"的处理，详见常见问题8
		 #   sed -i "s/-/\t0.0/g" result/alpha/alpha_rare.txt
		 
		 
		 #在进行此步骤的时候发现-alpha_div_rare不支持，因此暂且取消
		
		
		
		
		### 6.3. 筛选高丰度菌 Filter by abundance
		
		    #计算各特征的均值，有组再求分组均值，需根据实验设计metadata.txt修改组列名
		    #输入文件为feautre表result/otutab.txt，实验设计metadata.txt
		    #输出为特征表按组的均值-一个实验可能有多种分组方式
		    #-h显示脚本帮助(参数说明)
	Rscript ${db}/script/otu_mean.R -h
		    #scale是否标准化，zoom标准化总和，all输出全部样本均值，type计算类型mean或sum
	Rscript ${db}/script/otu_mean.R --input result/otutab.txt \
	  --metadata result/metadata.txt \
	  --group Group --thre 0 \
	  --scale TRUE --zoom 100 --all TRUE --type mean \
	  --output result/otutab_mean.txt
		    # 结果为全部和各组均值
	head -n3 result/otutab_mean.txt
		
	
	    #如以平均丰度>0.1%筛选，可选0.5或0.05，得到每个组的OTU组合
	awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
	    else {for(i=2;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
	    result/otutab_mean.txt > result/alpha/otu_group_exist.txt
	head result/alpha/otu_group_exist.txt
	wc -l result/alpha/otu_group_exist.txt
	    # 试一试：不同丰度下各组有多少OTU/ASV
	    # 可在http://ehbio.com/test/venn/中绘图并显示各组共有和特有维恩或网络图
	    # 也可在http://www.ehbio.com/ImageGP绘制Venn、upSetView和Sanky
	
	
	#   β多样性 Beta diversity(不支持)
	
	    #结果有多个文件，需要目录
	 #   mkdir -p result/beta/
	 #   基于OTU构建进化树 Make OTU tree, 4s
	 #   usearch -cluster_aggd result/otus.fa -treeout result/otus.tree
	    #生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
	  ##   -filename_prefix result/beta/
	
	
	 8. 物种注释分类汇总
	
	    #OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
	cut -f 1,4 result/otus.sintax \
	  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
	  > result/taxonomy2.txt
	head -n3 result/taxonomy2.txt
	
	    #OTU对应物种8列格式：注意注释是非整齐
	    #生成物种表格OTU/ASV中空白补齐为Unassigned
	awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
	 split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
	 print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
	 result/taxonomy2.txt > temp/otus.tax
	sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
	sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
	 > result/taxonomy.txt
    head -n3 result/taxonomy.txt
	 #统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
	mkdir -p result/tax
	for i in p c o f g;do
	  usearch -sintax_summary result/otus.sintax \
	  -otutabin result/otutab_rare.txt -rank ${i} \
	  -output result/tax/sum_${i}.txt
	done
	sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
	    # 列出所有文件
	wc -l result/tax/sum_*.txt
	head -n3 result/tax/sum_g.txt
	
	 9. 有参定量特征表
	
	    # 比对Greengenes97% OTUs比对，用于PICRUSt/Bugbase功能预测
	mkdir -p result/gg/
	
	    #方法1. usearch比对更快，但文件超限报错选方法2
	    # 默认10核以下使用1核，10核以上使用10核
	usearch -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fasta \
	  -otutabout result/gg/otutab.txt -threads 4
	    # 比对率80.0%, 1核11m，4核3m，10核2m，内存使用743Mb
	head -n3 result/gg/otutab.txt
	
	
	    # #方法2. vsearch比对，更准更慢，但并行24-96线程更强
	  #   vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fasta \
	  #     --otutabout result/gg/otutab.txt --id 0.97 --threads 12
	    # 比对率81.04%, 1核30m, 12核7m
	
	
	    #统计
	 #   usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat
	 #   cat result/gg/otutab.stat
	
	
	  #删除中间大文件
	rm -rf temp/*.fq
	
	    # 分双端统计md5值，用于数据提交
	>md5sum1.txt 
	>md5sum2.txt 
	md5sum *_1.fastq.gz > md5sum1.txt
	md5sum *_2.fastq.gz > md5sum2.txt
	paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | sed 's/*//g' > result/md5sum.txt
	rm md5sum*
	cat result/md5sum.txt
	    
	
