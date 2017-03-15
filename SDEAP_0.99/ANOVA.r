args <- commandArgs(trailingOnly = TRUE)
dat <- as.matrix(read.table(args[1], header=T,row.names = 1, check.names = F  , sep = "\t"))
c <- as.matrix(read.table(args[2], header=F, sep = "\t"))
outMTX <- matrix(, nrow = 0, ncol = ncol(dat))
pvals <- c()
for (i in 1:nrow(c)){
	
	data <- dat[i,]
	cluster <- c[i,]
	#cluster<-cluster[-length(cluster)]	
	cluster<-cluster +1
        tmpvec <- rep(1,1)
	pvalue <- 1	
	
        if(length(unique(cluster)) > 1){
				
		if(min(table(cluster))> 0  ){
			
			Q = data.frame(y = data, group = cluster)
		
			fit = lm(y ~ group, Q)
			res = anova(fit)
			
			pvalue <- res$"Pr(>F)"[-length(res$"Pr(>F)")]
			#tmpvec <- c( data,cluster, pvalue)
		}
		else{
			#tmpvec <- c( data,cluster, 0.9999)
		
		}
		
		
	}else{
		#tmpvec <- c( data,cluster, 0.9999)
		
	}
	#length(tmpvec)
	#names(tmpvec) <- NULL
	pvals<- append(pvals, pvalue)
        
}

FDR = p.adjust(pvals, "BH")


out_df <- data.frame(dat, pvals , FDR ,check.names = F )

outfname <-paste (args[1],".pvalues", sep = "")
write.table(out_df, outfname , sep="\t", quote= FALSE, col.names = T, row.names = T)



rname <- rownames(dat)

out_rname <- c()

for (i in 1:nrow(dat)){
	if(FDR[i]<as.numeric(args[3])){
		outMTX<-rbind(outMTX , dat[i,])
		out_rname <- append(out_rname , rname[i]) 
	}
}


rownames(outMTX)<- out_rname
outfname <-paste (args[1],".heatmap", sep = "")
write.table(outMTX,outfname , sep="\t", quote= FALSE, col.names = T, row.names = T)


