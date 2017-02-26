source('Data.R')


tiff(filename="output.tif")

## colour brewing
library(RColorBrewer)
g = 11
my.cols <- rev(brewer.pal(g, "RdYlBu"))
#compute 2D kernel density

# kernel density using MASS 
library(MASS)
z <- kde2d(X, Y, n=50)
plot(X, Y, xlab="X", ylab="Y", pch=19, cex=.3, col = "gray60")
contour(z, drawlabels=FALSE, nlevels=g, col=my.cols, add=TRUE, lwd = 2)
abline(h=mean(Y), v=mean(X), lwd=2, col = "black")
legend("topleft", paste("r=", round(cor(X, Y),2)), bty="n")

dev.off()
