#' @title HDR
#' @description BALABALA,BALABA
#' @details Input AND output
#' @param x A a data.
#' @param W The size.
#' @return A list
#' @export
#' @import
#' @importFrom
#' @examples
hello <- function(a = 3, b = 6) {
  result <- a * b
  print(result)
}


#记得设置工作目录
#内置的数据，需要处理自己的数据的话记得换路径
#数据预处理
#读取 OTUs 丰度表
setwd('data')

otu <- read.table('otutab.txt',row.names = 1,header=T)
otu1 <- read.table('otutab_1.txt',row.names = 1,header=T)

#过滤低丰度 OTUs 类群，它们对分类贡献度低，且影响计算效率
#例如剔除总丰度低于100 的值
otu <- otu[which(rowSums(otu,na.rm = TRUE,)>=100),]
dim(otu)

otu1 <- otu1[which(rowSums(otu1,na.rm = TRUE,)>=100),]
dim(otu1)
name=intersect(row.names(otu),row.names(otu1))
otu=otu[name,]
otu1=otu1[name,]
#合并有关于植物的信息
plant <- read.table("group.txt", row.names = 1,header=T)


otu <- data.frame(t(otu))
#otu <- otu[rownames(plant), ]
otu <- cbind(otu, plant)
otu_train <- otu
otu_train$group=factor(otu_train$group)
otu1 <- data.frame(t(otu1))


##randomForest 包的随机森林
install.packages('randomForest', repos = "https://mirrors.pku.edu.cn/CRAN/")
library(randomForest)

#随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
set.seed(123)
otu_train.forest <- randomForest(group~., data = otu_train, importance = TRUE)
otu_train.forest

importance_otu <- otu_train.forest$importance
head(importance_otu)



#作图展示 top30 重要的 OTUs
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'Top 30 - variable importance')

#可以根据某种重要性的高低排个序，例如根据“IncNodePurity”指标
importance_otu <- importance_otu[order(importance_otu[,4], decreasing = TRUE), ]
head(importance_otu)

##交叉验证辅助评估选择特定数量的 OTU
#5 次重复十折交叉验证
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu_train[-ncol(otu_train)], otu_train$group, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv

#提取验证结果绘图
install.packages('reshape2', repos = "https://mirrors.pku.edu.cn/CRAN/")
library(reshape2)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus),FUN = mean)
head(otu_train.cv.mean, 10)

#拟合线图
install.packages('ggplot2', repos = "https://mirrors.pku.edu.cn/CRAN/")
library(ggplot2)

p=ggplot(otu_train.cv.mean, aes(Group.1, x)) +
  geom_line() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')

p=p+geom_vline(aes(xintercept=23), colour="#BB0000", linetype="dashed")
p


#然后取出排名靠前的 OTU
importance_otu.select <- importance_otu[1:23, ]
importance_otu.select

#输出表格
write.table(importance_otu.select, 'importance_otu.select.txt', sep = '\t', col.names = NA, quote = FALSE)

#不妨简单查看下这些重要的 OTU 丰度与植物的关系
#可以看到趋势非常明显
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu[ ,c(otu_id.select, 'group')]

otu_train.select <- otu.select
#随机森林计算（默认生成 500 棵决策树）
set.seed(123)
otu_train.select.forest <- randomForest(factor(group)~., data = otu_train.select, importance = TRUE)
otu_train.select.forest

install.packages('pROC', repos = "https://mirrors.pku.edu.cn/CRAN/")
library(pROC)

fc<-as.numeric()
mod_pre<-as.numeric()
model <- randomForest(factor(group)~., data = otu_train.select, importance = TRUE)

model_pre<-predict(model,type="prob")
fc<-append(fc,as.numeric(ifelse(otu_train.select$group=="HHM",0,1)))
mod_pre<-append(mod_pre,model_pre[,1])
df<-cbind(fc,as.numeric(mod_pre))

# ROC curve
x<-plot.roc(df[,1],df[,2],
            smooth=F,
            lwd=2,
            ylim=c(0,1),
            xlim=c(1,0),
            legacy.axes=T,
            main="",
            col="seagreen3")

# the AUC value
x[["auc"]]
legend.name <- c("Area under the curve: 0.808")
legend("bottomright",
       legend=legend.name,
       lwd = 2,
       col = "seagreen3",
       bty="n")



t=predict(model,newdata=otu1)
write.table(t,"predict.txt")
