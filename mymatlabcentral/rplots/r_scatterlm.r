source('Data.R')

library(ggplot2)

tiff(filename="output.tif")
qplot(x, y, xlab="X (units)", ylab="Y (units)", geom=c("point", "smooth"), method="lm", formula= y~x) + theme_bw()
dev.off()


# http://www.cookbook-r.com/Graphs/Scatterplots_(ggplot2)/

# http://scs.math.yorku.ca/index.php/R_Graphs_Gallery

 library(lattice)
 data(Davis, package="car")
 # Package -- lattice
 xyplot(weight ~ height, data=Davis, subset=-12,
 			 main="Davis Data - Weight by Height",
 			 sub="Data collected by Davis, data set available in package car",
 			 panel = function(x, y) {
 			 	panel.grid(h = -1, v = 2) # -1 align horiz grid to ticks, 2 indicates 2 vertical lines
 			 	panel.xyplot(x, y, col='purple')
 			 	panel.loess(x, y, span=0.5, col='red')  # Span is the loess smoothness parameter, default 2/3 
 			 	panel.lmline(x, y)  # Adds a least squares line 
 			 } )