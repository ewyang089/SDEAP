library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
png(filename="reg.png")

in_data <- read.table(args[1], header=T,row.names=1 , sep = "\t",as.is=TRUE,check.names=F)
in_mtx <- as.matrix(in_data)
rname<- rownames(in_mtx)
M_mean <-c()
cv2 <- c()
name<-c()

#print(rname)
for(i in 1:nrow(in_data)){
	#print(mean(in_mtx[i,]))
        tmp_mean <- mean(in_mtx[i,])
	#print(i)
        #print(in_mtx[i,])
        if(tmp_mean>0){
	        if(log2(tmp_mean) > -2){

		 M_mean<-append(M_mean, tmp_mean)
		 cv2<-append(cv2, var(in_mtx[i,])/(tmp_mean*tmp_mean))
	         name <- append(name ,rname[i])
		}
	}
}


df <- data.frame(name,M_mean,cv2)
df<- df[order(df$M_mean),]
sim_fit <- nls(df$cv2 ~  a*df$M_mean^(-1) + b, start=list(a=1, b=0),data = df)
summary(sim_fit)
plot( log2(df$M_mean) , df$cv2 )
lines(log2(df$M_mean) ,predict(sim_fit), col="red" )

pred <- predict(sim_fit)

outdf <- data.frame(df$name ,  df$M_mean ,  df$cv2 , pred  )
colnames(outdf) <- c("name","M_mean", "cv2","pred" )
#outdf <- outdf[which(outdf$cv2/outdf$pred > 2 & outdf$M_mean > 0.5  ),]
outdf <- outdf[which(outdf$cv2/outdf$pred >as.numeric( args[3]) & outdf$M_mean > as.numeric(args[4])  ),]

#print(outdf$name) 

in_data<- in_data[rownames(in_data) %in% outdf$name, ]

write.table(in_data, args[2]  ,sep="\t",quote=F, col.names = T, row.names = T)
dev.off()
