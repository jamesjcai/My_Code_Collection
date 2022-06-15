#setwd('G:\\My Drive\\Collab\\Li_SCA\\matlab_processed_data\\sctenifold_net_analysis\\knk_research')
library("rhdf5")
setwd('C:\\Users\\jcai\\OneDrive\\h5test')

d <- h5read(file = "my_example_filex.h5", name = "/c")
a <- h5read(file = "my_example_filex.h5", name = "/s")

#c <- h5read(file = "test.mat", name = "/c")
#g <- h5read(file = "teststring.mat", name = "/g")

d2 <- h5read(file = "aaa.mat", name = "/testdat1")
# a2 <- h5read(file = "aaa.mat", name = "/testdat2")

a2<-read.table('genelist.txt')
