setwd("G:\\垂体瘤GEO数据库分析\\Gene expression profiling of non-invasion and invasion NFPAs\\up regulation_invasion vs noninvasion")
pathway<-read.csv("KEGG pathway.csv")
par(mar=c(5,25,4,2))
barplot(rev(pathway$X..LOG10.Pvalue.),col="blue",
        horiz = TRUE,xlab = "-log10(PValue)")
text(seq(from=0.7,length.out=10,by=1.2),x=-1.75,label=rev(pathway$Term),xpd=TRUE)
title(main="up-regulated",cex.main=2.0)

#以下是绘制下调通路
setwd("G:\\垂体瘤GEO数据库分析\\Gene expression profiling of non-invasion and invasion NFPAs\\down regulation_invasion vs noninvasion")
pathway<-read.csv("KEGG pathway.csv")
par(mar=c(5,25,4,2))
barplot(rev(pathway$X.log10.Pvalue.),col="blue",
        horiz = TRUE,xlab = "-log10(PValue)")
text(seq(from=0.7,length.out=4,by=1.2),x=-1.2,label=rev(pathway$Term),xpd=TRUE)
title(main="down-regulated",cex.main=2.0)
