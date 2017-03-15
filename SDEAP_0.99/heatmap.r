library(made4)


#par(ps = 20, cex = 1, cex.main = 1)
args <- commandArgs(trailingOnly = TRUE)
#png(filename="heatmap.png")
pdf()

mydata <- as.matrix(read.table(args[1], header=T, sep = "\t",row.names = 1,as.is=TRUE))

#mydata <- mydata[,-1]
sTree <- heatplot(mydata, lowcol="blue", highcol="red",scale="none" , dend= "both", returnSampleTree=TRUE  , cexCol=0.2 , cexRow=0.2)

dev.off()