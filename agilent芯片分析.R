library(GEOquery)
library(limma)
setwd("G:\\垂体瘤GEO数据库分析\\Gene expression profiling of non-invasion and invasion NFPAs")
#获取数据
gse51618<-getGEO(filename = "GSE51618_series_matrix.txt.gz")
gpl6480<-getGEO('GPL6480',destdir = ".")
#制作探针注释文件
colnames(Table(gpl6480))
Table(gpl6480)[1:10,1:8]
write.csv(Table(gpl6480)[,c(1,7)],"GPL6480.csv",row.names=F)
genename<-read.csv("GPL6480.csv")
#构建表达矩阵
exprSet<-as.data.frame(exprs(gse51618))#得到矩阵，以下代码为将行名转换为ID
exprSet$ID<-rownames(exprSet)
express<-merge(x=exprSet,y=genename,by="ID",all.x = T)
express$ID<-NULL
#删除重复的基因，保留每个基因最大表达量的结果
rowMeans<-apply(express,1,function(x)mean(as.numeric(x),na.rm=T))
express<-express[order(rowMeans,decreasing = T),]
express<-express[!duplicated(express[,11]),]#express第11列为基因名
rownames(express)<-express[,11]
express<-express[,-11]
#构建分组矩阵
pdata<-pData(gse51618)
group_list<-subset(pdata,select=title)
group_list$condition<-c("nin","nin","nin","nin","i","i","i","n","n","n")
design<-model.matrix(~0+factor(group_list$condition))
colnames(design)<-levels(factor(group_list$condition))
rownames(design)<-colnames(express)
            ##以下代码为构建tumor-normal分组矩阵
pdata<-pData(gse51618)
group_list<-subset(pdata,select=title)
group_list$condition<-c("t","t","t","t","t","t","t","n","n","n")
design<-model.matrix(~0+factor(group_list$condition))
colnames(design)<-levels(factor(group_list$condition))
rownames(design)<-colnames(express)
#构建差异比较矩阵
contrast.matrix<-makeContrasts(contrasts = c("i-nin","i-n","nin-n"),levels=design)
fit<-lmFit(express,design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
            ##以下代码为构建tumor-normal比较矩阵
contrast.matrix<-makeContrasts(contrasts = c("t-n"),levels=design)
fit<-lmFit(express,design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
#得到差异表达结果
x<-topTable(fit2,coef = "i-nin",n=Inf,adjust.method = "BH",sort.by = "P")
sum(x$adj.P.Val<0.05)
write.csv(x,"invasion vs non-invasion_total_DEG.csv",quote = F)
re<-x[x$adj.P.Val<0.05 & (x$logFC>1 | x$logFC< -1),]
write.csv(re,"invasion vs non-invasion_DEG.csv",quote = F)

x<-topTable(fit2,coef = "i-n",n=Inf,adjust.method = "BH",sort.by = "P")
sum(x$adj.P.Val<0.05)
re<-x[x$adj.P.Val<0.05 & (x$logFC>1 | x$logFC< -1),]
write.csv(re,"invasion vs normal_DEG.csv",quote = F)

x<-topTable(fit2,coef = "nin-n",n=Inf,adjust.method = "BH",sort.by = "P")
sum(x$adj.P.Val<0.05)
re<-x[x$adj.P.Val<0.05 & (x$logFC>1 | x$logFC< -1),]
write.csv(re,"non-invasion vs normal_DEG.csv",quote = F)

x<-topTable(fit2,coef = "t-n",n=Inf,adjust.method = "BH",sort.by = "P")
sum(x$adj.P.Val<0.05)
re<-x[x$adj.P.Val<0.05 & (x$logFC>1 | x$logFC< -1),]
write.csv(re,"tumor vs normal_DEG.csv",quote = F)
