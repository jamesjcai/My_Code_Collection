##data <- read.table("data.csv", header = FALSE, sep = ",")
##results <- shapiro.test(data$V1)
##results2 <- c( results$statistic[["W"]], results$p.value )
##write.table(results2, file="testResults.csv", sep = ",", col.names = FALSE, row.names = FALSE, qmethod = "double")

#setwd("C:/Users/JCai/Desktop/geo/model3popgeno2")

dglm.Pvalues <- function(dglm.fit){
	P.disp = anova.dglm(dglm.fit)$Adj.P[2]
	P.mean = summary(dglm.fit)$coef[2,4]
	list(P.mean=P.mean, P.disp=P.disp)
}

##An example using dglm.Pvalues(.) with 200 simulated observations
## set.seed(123)

set.seed(12323)

require(dglm)

pop <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2);


files <- list.files(path="/scratch/jcai/input",pattern="*",full.names=FALSE)


for (fileidx in 1:length(files)){
        print(fileidx);

	output_x <- sprintf("output_%d.txt",fileidx);
	input_x <- sprintf("input_%d.txt",fileidx);
	index_x <- sprintf("index/%d.txt",fileidx);

	if (file.exists(index_x)) next;
	file.copy("index/seed",index_x,overwrite=TRUE);

	#unlink(output_x,recursive=FALSE);
	#unlink(input_x,recursive=FALSE);


	file.copy(paste("/scratch/jcai/input/",files[fileidx],sep=""),input_x,overwrite=TRUE);
	x <- read.table(input_x,header=TRUE);
	y <- x[,1];
	#pop <- x[,2];
	#sex <- x[,3];
	a<-dim(x);
	a<-a[2];

	for (i in c(4:a)){
		gen <- x[,i];

		#result <- try(d2.fit <- dglm( formula = y ~ gen, dformula = ~ gen));
		#if(class(result)[1] == "try-error") next;

		d2.fit <- try(dglm( formula = y ~ gen + pop, dformula = ~ gen), TRUE);
		if(class(d2.fit)[1] == "try-error") next;

		P.values <- try(dglm.Pvalues( d2.fit ),TRUE);
		if(class(P.values)[1] == "try-error") next;

		Pres <- list(i-3, P.values);
		write.table(Pres, file=output_x, sep = "\t", col.names = FALSE, row.names = FALSE, qmethod = "double",append=TRUE)
		#if (P.values[2] < 0.05) {print( P.values )}
	}
	file.copy(output_x,paste("/scratch/jcai/output_x/",files[fileidx],sep=""),overwrite=TRUE);
	unlink(output_x,recursive=FALSE);
	unlink(input_x,recursive=FALSE);
}
