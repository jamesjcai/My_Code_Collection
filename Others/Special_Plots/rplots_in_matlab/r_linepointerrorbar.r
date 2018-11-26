source('Data.R')

require(ggplot2)
require(grid)


tiff(filename="output.tif")



#pdf("output.pdf")
#plot(x,y, main="PDF Scatterplot Example", 
#col=rgb(0,100,0,50,maxColorValue=255), pch=16)
#dev.off()

myd <- data.frame (X = xx,
                   Y = yy,
                   group = rep (c("Kidney", "Gill"), each = 12),
                   error = rep (c(2.5, 3.0), each = 12))

# group = rep (c("Kidney", "Gill"), each = 12)
# error = rep (c(2.5, 3.0), each = 12)

# line and point plot
f1 = ggplot(data = myd, aes(x = X, y = Y, group = group) )  # lesion becomes a classifying factor
f2 <- f1 + geom_errorbar(aes(ymin = Y - error, ymax = Y + error), width=0.3) +
	geom_line() + geom_point(aes(shape=group, fill=group), size=5)

 f3 <- f2 +  scale_x_continuous("X (units)", breaks=1:12) +
     scale_y_continuous("Y (units)", limits = c(0, 40), breaks=seq(0, 40, by = 5)) +
     scale_shape_manual(values=c(24,21)) +
     scale_fill_manual(values=c("white","black")) +
     stat_abline(intercept=0, slope=0, linetype="dotted") +
     annotate("text", x=11, y=10, label="X") +
     theme_bw()

   optns <- theme (
          plot.title = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=12),
          axis.title.y = element_text(face="bold", size=12, angle=90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.2,0.8),
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          legend.key.size = unit(1.5, "lines"),
          legend.key = element_blank()
     )
f3 +  ggtitle ( "MY awsome plot for publication") + optns

dev.off()