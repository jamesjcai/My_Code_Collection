source('Data.R')

library(hexbin)

tiff(filename="output.tif")
bin<-hexbin(x, y, xbins=50) 
plot(bin, main="Hexagonal Binning")
dev.off()


#pdf("output.pdf")
#plot(x,y, main="PDF Scatterplot Example", 
#col=rgb(0,100,0,50,maxColorValue=255), pch=16)
#dev.off()
