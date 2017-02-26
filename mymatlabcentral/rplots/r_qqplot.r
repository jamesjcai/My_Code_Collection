source('Data.R')
tiff(filename="output.tif")
qqplot(x, y)
dev.off()


