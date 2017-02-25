source('Data.R')

tiff(filename="output.tif")
#pdf(file="output.pdf")
plot(x, y, type = "p",  xlim = NULL, ylim = NULL, log = "", main = NULL, sub = NULL, xlab = "Ordered MD^2", ylab = "Chi-square", cex=1.8, cex.lab=1.8)
#plot(x, y, type = "p",  xlim = NULL, ylim = NULL, log = "", main = NULL, sub = NULL, xlab = "Gene 1", ylab = "Gene 2", cex=1.8, cex.lab=1.8)


points(tail(x, n=3),tail(y, n=3),type="p",col="red",pch=8,cex=2.2)
dev.off()


#library(ggplot2)
#tiff(filename="output.tif")
#qplot(x, y, xlab="X (units)", ylab="Y (units)")+ theme_bw()
#dev.off()

