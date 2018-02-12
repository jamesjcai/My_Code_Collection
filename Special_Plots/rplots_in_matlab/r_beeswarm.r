source('Data.R')

require(beeswarm)

tiff(filename="output.tif")
beeswarm(X ~ G, method = 'swarm', pch = 16, pwcol = as.numeric(c),
         xlab = '', ylab = 'xv', labels = c('bimodal', 'uniform'))
#legend('topright', legend = LETTERS[1:2], title = 'ER', pch = 16, col = 1:2)
legend('topright', legend = c("XXX","YYY"), title = 'ER', pch = 16, col = 1:2)

dev.off()
# http://denishaine.wordpress.com/2011/09/16/beeswarm-plot-with-ggplot2/