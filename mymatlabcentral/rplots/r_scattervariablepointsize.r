source('Data.R')

library(ggplot2)

tiff(filename="output.tif")

mydf <- data.frame(x = x,
                   y = y,
                   count = c)

ggplot(mydf, aes(x = x, y = y)) + geom_point(aes(size = count),shape=21) + geom_point(aes(colour = count)) + scale_shape(solid = FALSE)+ theme_bw()

#ggplot(mydf, aes(x = x, y = y)) + 
#  geom_point(aes(size = count)) +
#  scale_size_continuous(range = c(3, 7))
# ggsave(file = '2013-11-25.png', height = 5, width = 5)

dev.off()


