library(ggplot2)
data<-read.csv("invasion vs non-invasion_total_DEG.csv")
#设置数据分组
data$label<-ifelse(
  data$logFC< -1 & data$P.Value<0.05,0,
  ifelse(
    data$logFC>1 & data$P.Value<0.05,2,1))
#设置横轴及纵轴
r03<-ggplot(data,aes(logFC,-1*log10(P.Value),color=as.factor(label)))
#显示火山图
r03+geom_point()+scale_color_manual(values =c("green","black", "red"))
#设置坐标轴范围
r04<-r03+geom_point()+scale_color_manual(values =c("green","black", "red"))+
  xlim(-10,10)+ylim(0,7.5)+labs(title = "Volcanoplot")+theme(plot.title = element_text(hjust = 0.5))
#添加阈值线
volcano<-r04+geom_hline(yintercept = -1*log10(0.05))+geom_vline(xintercept = c(-1,1))
#保存
ggsave("volcano.png",volcano)
