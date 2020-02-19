setwd("G:\\垂体瘤GEO数据库分析\\Gene expression profiling of non-invasion and invasion NFPAs\\up regulation_invasion vs noninvasion")
data<-read.csv("KEGG pathway.csv")
p<-ggplot(data = data,aes(-log10(PValue),Term))
pp<-p+geom_point(aes(size=Count,color=-log10(PValue)))

pp+scale_color_gradient(low = "green",high = "red")+
  labs(color="-log10(P)",size="gene number",x="-log10(P)",title = "TOP10 pathway")+
  theme(plot.title = element_text(hjust = 0.5))