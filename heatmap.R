library(gplots)
data<-express[,1:7]
heatmap.2(as.matrix(data),col = redgreen(75),scale = "row",
          key = TRUE,symkey = FALSE,density.info = "none",
          trace = "none",cexRow=0.5)
