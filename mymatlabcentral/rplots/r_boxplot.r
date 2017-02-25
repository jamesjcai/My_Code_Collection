source('Data.R')

tiff(filename="output.tif")

# boxplot(X~G, main="Title", xlab="XLabel", ylab="YLabel")

 df=data.frame(X,G);
 df$type <- with(df, reorder(G, X, median))
 
 boxplot(X ~ df$type, data = df,
          xlab = "XLabel", ylab = "YLabel",
          main = "Title", varwidth = FALSE,
          col = "lightgray")
dev.off()
