
# weights only for rare variants
# 	Z1: genotypes of rare variants
# 	Z2: genotypes of common variants
# 		To use flat weight set weights.beta=c(1,1)

G_Inverse<-function(S, tol = sqrt(.Machine$double.eps))
{
  dnx <- dim(S)
  if(dnx[1] == 1){
  	S_inv<-1/S
  	rank<-1
  	
  } else {
  
  	s <- eigen(S,symmetric=TRUE)
  	nz <- which(s$values > tol * s$values[1])
  	S_inv<- s$vectors[,nz] %*% (t(s$vectors[,nz]) / s$values[nz])
  	rank<-length(nz)
  }
  
  re<-list(S_inv = S_inv, rank=rank)
  return(re)
}


Hotellings_T<-function(X, Y){


	nx<-dim(X)[1]
	ny<-dim(Y)[1]
	
	X1<-t(X) -  colMeans(X)
	Y1<-t(Y) - colMeans(Y)
	xy<- colMeans(X) - colMeans(Y)
	
	S = (X1 %*% t(X1)  + Y1 %*% t(Y1))/(nx + ny -2)
	S_inv<-G_Inverse(S)
	v<-S_inv$rank
	
	T = nx * ny / (nx + ny) * ( t(xy) %*% S_inv$S_inv %*% xy)
	F<-(nx+ny -v-1)/v / (nx + ny -2) * T
	p.val<-1- pf(F, df1=v, df2=nx+ny-v-1)
	return(list(p.val=p.val, df1=v, df2=nx+ny-v-1, T2=T, F=F))
}


#
#	burden 1: weighted burden
#	burden 2: CAST (dichotomized)
#
CMC_Method<-function(Z1, Z2, y, weights.beta=c(1,25), weights = NULL){


	n<-length(y)
	if(is.null(weights)){
		MAF<-colMeans(Z1)/2
		weights<-SKAT:::Beta.Weights(MAF,weights.beta)
	}
	
	Z1 = t(t(Z1) * (weights))
	burden1<-rowSums(Z1)
	burden2<-rep(0,n)
	IDX<-which(burden1 > 0)
	
	if(length(IDX) > 0){
		burden2[IDX]<-1
	}
	
	Z2.all1<-cbind(Z2, burden1)
	Z2.all2<-cbind(Z2, burden2)
	
	# Case/Control
	ID_CASE<-which(y==1)
	ID_Control<-which(y==0)
	
	A1.1<-Z2.all1[ID_CASE,]
	A1.2<-Z2.all1[ID_Control,]

	A2.1<-Z2.all2[ID_CASE,]
	A2.2<-Z2.all2[ID_Control,]
	
	out1<-Hotellings_T(A1.1, A1.2)
	out2<-Hotellings_T(A2.1, A2.2)

	return(list(p.val1=out1$p.val, p.val2=out2$p.val))
} 

Burden_Method<-function(Z, y, weights.beta=c(1,25), weights = NULL){


	n<-length(y)
	if(is.null(weights)){
		MAF<-colMeans(Z)/2
		weights<-SKAT:::Beta.Weights(MAF,weights.beta)
	}
	
	Z = t(t(Z) * (weights))
	burden1<-rowSums(Z)
	
	l<-glm(y~ burden1, family=binomial)
    
    p.val<-summary(l)$coefficients[2,4]
    
    return(list(p.value=p.val))
    
}


