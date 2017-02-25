source('Data.R')

require(beeswarm)

tiff(filename="output.tif")

method<-2;


if (method==1){
	boxplot(X ~ G, outline = FALSE, main = 'boxplot + beeswarm', labels = c('bimodal', 'uniform'))
	beeswarm(X ~ G, col = 1, pch = 16, add = TRUE)
}

if (method==2){
	beeswarm(X ~ G, method = 'swarm', pch = 16, pwcol = as.numeric(c),
		 xlab = '', ylab = 'xv', labels = c('bimodal', 'uniform'))
	#legend('topright', legend = LETTERS[1:2], title = 'ER', pch = 16, col = 1:2)
	legend('topright', legend = c("XXX","YYY"), title = 'ER', pch = 16, col = 1:2)
	#legend('topright', inset=.05, legend = c(paste("P =", 0.00001)), bty='n', cex=1.5)

	bxplot(X ~ G, add = TRUE)
}

dev.off()
# http://denishaine.wordpress.com/2011/09/16/beeswarm-plot-with-ggplot2/
# http://www.cbs.dtu.dk/~eklund/beeswarm/
