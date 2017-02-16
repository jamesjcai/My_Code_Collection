source('Data.r')
options(warn=-1)
require(dglm,quiet=TRUE)

d2.fit=dglm(formula = y~x, dformula = ~x)
p.disp = anova.dglm(d2.fit)$Adj.P[2]
#p.mean = summary(d2.fit)$coef[2,4]
cat(p.disp)
