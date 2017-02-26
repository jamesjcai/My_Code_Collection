#Load required R.matlab package
#library(R.matlab)
 
#Load data saved from Matlab
source('Data.R')
tiff(filename="output.tif")

a<-density(x1)
b<-density(y1)
c<-density(x2)
d<-density(y2)

plot(a,type="l",col="blue",lwd=2,ylim=c(0, 0.4),xlab='MAF',main='')
lines(b,type="l",col="red",lwd=2)
lines(c,type="l",col="blue",lwd=2,lty=3)
lines(d,type="l",col="red",lwd=2,lty=3)

legend("topright", c("a","b","c","d"),inset=.05,
    lty=c(1, 1, 3, 3),col=c('blue','red','blue','red'))


#require('lattice')
#dat <- data.frame(dens = c(x, y)
#                   , lines = rep(c("a", "b"), each = 100))
#densityplot(~dens,data=dat,groups = lines,
#            plot.points = FALSE, ref = TRUE, 
#            auto.key = list(space = "right"))
