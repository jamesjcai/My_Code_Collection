#Load required R.matlab package
#library(R.matlab)
 
#Load data saved from Matlab
source('Data.R')
tiff(filename="output.tif")


# Add a Kernel Density Plot
myhist <- hist(x, breaks=20, density=20, col=NULL, xlab="Accuracy", main="Overall")
multiplier <- myhist$counts / myhist$density
mydensity <- density(x)
mydensity$y <- mydensity$y * multiplier[1]

#plot(myhist)
lines(mydensity,type = "l",col="red",lwd=2)



# Add a Normal Curve
#h<-hist(x, breaks=10, col="red", xlab="Miles Per Gallon", 
#  	 main="Histogram with Normal Curve") 
xfit<-seq(min(x),max(x),length=40) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(myhist$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=1)


dev.off()




#    m<-mean(g)
#    std<-sqrt(var(g))
#    hist(g, density=20, breaks=20, prob=TRUE, 
#         xlab="x-variable", ylim=c(0, 2), 
#         main="normal curve over histogram")
#    curve(dnorm(x, mean=m, sd=std), 
#          col="darkblue", lwd=2, add=TRUE, yaxt="n")