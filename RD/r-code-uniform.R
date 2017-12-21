iso <- read.table("liver_uc_1",head=F,sep="\t",stringsAsFactors=F)

nnn <- 0
enm <- 0
SG <- 0
Rsquare <- 0
slope <- 0
TL <- 0
pdf("isoform_regress_uc.pdf",width=5,height=5)
par(mfrow=c(1,1))
   ID <- names(table(iso[,1]))  
   for(iii in ID)  #  23787
   {
    sub <- iso[iso[,1]== iii, ]
    infor <- sub[1,]
    tot <- sub[-1,]
    if(sum(tot[,6])>200) 
    { if(infor[1,3]=="+")
      {tot <- tot[nrow(tot):1,]} 
      name <- strsplit(infor[1,7],",")[[1]]
      INDEX <- matrix(0,nrow=nrow(tot),ncol=length(name))
      DIST <- matrix(0,nrow=nrow(tot),ncol=length(name))
      lexon <- rep(0,nrow(tot))
      exon_c <- rep(0,nrow(tot))
      for(k in 1:nrow(tot))
      {
          INDEX[k,] <- as.numeric(strsplit(tot[k,7],",")[[1]])
          exon_c[k] <- tot[k,6]
          lexon[k] <- (tot[k,5]-tot[k,4] + 1)
          if(k==1)
          {DIST[k,] <- lexon[k] * INDEX[k,]/2  }
          if(k>1)    
          {DIST[k,] <- INDEX[k,] * (lexon[1:k]%*% INDEX[1:k,] - lexon[k] * INDEX[k,]/2)}
      }
      # ltranscript <- rep(0,length(name))
      ### EM  
        filter <- ( (exon_c >0 ) &  (lexon>150) ) 
        if(sum(filter)>2)
        {
          nnn <- nnn+1 
          tname <- paste(iii," (",name[1]," ",infor[1,2]," ",infor[1,3],")",sep="")
          dist <- DIST[filter,1]
          dens <- log(exon_c[filter]/lexon[filter])
          plot(dist,dens,xlim=c(0,sum(lexon)),xlab=expression(d[gj]),ylab=expression(log(N[gj]/l[gj])),main=tname,type="p")
          ndist <- DIST[!filter,1]
          if(length(ndist)>0)
          { 
              ypos <- min(dens) - (max(dens)-min(dens)) *0.015 
              points(ndist,rep(ypos,length(ndist)),pch=4) 
          }
          e.lm <- lm(dens~dist)
          abline(e.lm,col="blue",lwd=2) 
          sss <-summary(e.lm)          
          infor <-round(c(coefficients(e.lm), sss$r.squared ),4 )
            
          enm[nnn] <- sum(filter)
          Rsquare[nnn] <- infor[3]
          slope[nnn] <-  infor[2]       
          TL[nnn] <- sum(lexon)
          SG[nnn] <-  sss$sigma
          xpos <- 0.2*sum(lexon)
          ypos <- min(dens) + (max(dens)-min(dens)) *0.5
          ccc <- paste("R-square:", infor[3],sep=" ")
          text(xpos,ypos,ccc,offset=1)  
          ypos <- min(dens) + (max(dens)-min(dens)) *0.4
          ccc <- paste("Slope:", infor[2],sep=" ")
          text(xpos,ypos,ccc,offset=1)  
          ypos <- min(dens) + (max(dens)-min(dens)) *0.3
          ccc <- paste("Intercept:", infor[1],sep=" ")
          text(xpos,ypos,ccc,offset=1) 
          ypos <- min(dens) + (max(dens)-min(dens)) *0.2 

      }    
     } ## tot >? 
   } ## 
par(mfrow=c(1,1)) 
dev.off()

reg <- data.frame(enm= enm,Rsquare=Rsquare,slope=slope,TL=TL)

write.table(reg,"uniformexon.txt",quote=F,sep="\t",col.names=T,row.names=F)

pdf("isoform_regress_summary_uc.pdf",width=5,height=5)
hist(reg[,2],main="",breaks=seq(0,1,0.1),xlab="R-Square")
abline(v=median(reg[,2]),col="blue",lwd=2,lty=2)
ym <- -1*reg[,3]
xm <- reg[,4]
y<- ym[(ym>0) & (reg[,2]>0.7) ]
x<- xm[(ym>0) & (reg[,2]>0.7) ]

plot(x,y,xlab="Transcript Length",ylab=expression(alpha[g]),type="p",cex=0.2,col="grey")

y.loess <- loess(y ~ x, span=0.75, data.frame(x=x, y=y))
y.predict <- predict(y.loess, data.frame(x=x[order(x)]))
lines(x[order(x)],y.predict,col="blue",lwd=2)
hist(y,main="",xlab= expression(alpha[g]))
abline(v=median(y),col="blue",lwd=2,lty=2)
abline()
dev.off()



print( summary(SG))

print( summary(SG[(ym>0) & (reg[,2]>0.7)]))



