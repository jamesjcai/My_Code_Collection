source('Data.R')

tiff(filename="output.tif")
boxplot(X~G, main="Title", xlab="XLabel", ylab="YLabel")
dev.off()
