source('Data.R')


tiff(filename="output.tif")

Method<-2;

# http://www.inside-r.org/packages/cran/car/docs/Ellipses
# Method 1
if (Method==1){
    library(car)
    dataEllipse(x[,1],y[,1], levels=0.95, lty=2)
}

# Method 2
if (Method==2){
	library(ggplot2)
	library(devtools)
	library(digest)
	source('stat-ellipse.r')
	df <- data.frame(x=x, y=y, group="A");
	qplot(data=df, x=x, y=y)+stat_ellipse() + theme_bw()
}


# Method 3
if (Method==3){

	df <- data.frame(x=x, y=y, group="A");

	library(ellipse)
	df_ell <- data.frame()
	for(g in levels(df$group)){
	df_ell <- rbind(df_ell, cbind(as.data.frame(with(df[df$group==g,], ellipse(cor(x, y), 
						 scale=c(sd(x),sd(y)), 
						 centre=c(mean(x),mean(y))))),group=g))
	}

	library(ggplot2)
	ggplot(data=df, aes(x=x, y=y)) + geom_point(size=1.5, alpha=.6) + geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=1, linetype=2) + theme_bw()

}

dev.off()


# http://stackoverflow.com/questions/2397097/how-can-a-data-ellipse-be-superimposed-on-a-ggplot2-scatterplot


