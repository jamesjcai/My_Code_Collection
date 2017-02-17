#gsSKAT Package functions

#####################################################
##P-value based analyses with linear kernel functions
#####################################################

gmt.to.design <- function(gmt.file, res.genes, min.gene, max.gene){
  require(Matrix)
  gmt.lst <- list()
  for(i in 1:length(gmt.file)){
    gmt.lst[[i]] <- readLines(gmt.file[i])
  }
  gmt.dat <- do.call(c,gmt.lst)
  gs.names <- c()
  gs.genes <- list()
  for(i in 1:length(gmt.dat)){
    dat.i <- unlist(strsplit(gmt.dat[i],"\t"))
    gs.names[i] <- dat.i[1]
    gs.genes[[i]] <- dat.i[-c(1:2)][which(dat.i[-c(1:2)]%in%res.genes[,1])]
  }
  #Which sets meet min/max gene criterion?
  gs.len <- unlist(lapply(gs.genes,length))
  which.out <- which(gs.len>=min.gene&gs.len<=max.gene)
  genes.used <- unique(unlist(gs.genes[which.out]))
  sparse.i <- c()
  sparse.j <- c()
  sparse.x <- c()
  for(i in 1:length(which.out)){
    col.i <- which(genes.used%in%gs.genes[[which.out[i]]])
    n.i <- length(col.i)
    row.i <- rep(i,n.i)
    sparse.i <- append(sparse.i,row.i)
    sparse.j <- append(sparse.j,col.i)
    sparse.x <- append(sparse.x,rep(1/sqrt(n.i),n.i))
  }
  D.mat <- sparseMatrix(i = sparse.i,j = sparse.j, x = sparse.x)
  gs.names.out <- gs.names[which.out]
  gene.ord <- rep(NA,length(genes.used))
  for(i in 1:length(gene.ord)){
    gene.ord[i] <- which(res.genes[,1]==genes.used[i])
  }
  colnames(D.mat) <- genes.used
  rownames(D.mat) <- gs.names.out
  return(list(D.mat = D.mat,genes.used = genes.used, gs.used = gs.names[which.out],gene.ord = gene.ord))
}


n.eff.test.gal <- function(D.mat){
  D.mat.lam <- eigen(tcrossprod(D.mat),symmetric = TRUE, only.values = TRUE)$values
  D.mat.lam.nz <- D.mat.lam[which(D.mat.lam>0)]
  return(sum(sqrt(D.mat.lam.nz))^2 / sum(D.mat.lam.nz))
}

n.eff.test.gao <- function(D.mat){
  #Based off Gao et al. 2009
  D.mat.lam <- abs(eigen(tcrossprod(D.mat), symmetric = TRUE, only.values = TRUE)$values)
  out <- min(which(cumsum(D.mat.lam)/sum(D.mat.lam)>0.995))
  return(out)
}

gsSKAT <- function(p.gene, path.gmt, min.gene = 10, max.gene = 500){
  ##Arguments
  #p.gene: Gx2 data.frame where first column corresponds to the gene names, second column gene-level p-values
  #path.gmt: character string or vector of character strings pointing to .gmt file(s) to load
  #min.gene: minimum number of genes in a gene set to conduct test
  #max.gene  maximum number of genes in a gene set to conduct test
  
  #SKAT/SKAT-O will sometimes report p-values of 1
  #To avoid -Inf Z scores, we set these to 0.95
  if(sum(p.gene[,2]==1)>0)
    p.gene[which(p.gene[,2]==1),2] <- 0.95
  
  gs.info <- gmt.to.design(path.gmt, p.gene, min.gene = min.gene, max.gene = max.gene)
  n.gs <- length(gs.info$gs.used)
  n.eff <- n.eff.test.gao(gs.info$D.mat)
  n.gene <- rep(NA,n.gs)
  z.gs <- rep(NA,n.gs)
  p.gs <- rep(NA,n.gs)
  
  for(i in 1:n.gs){
    genes.i <- gs.info$gene.ord[which(gs.info$D.mat[i,]!=0)]
    pval.i <- p.gene[genes.i,2]
    gs.z.i <- sum(qnorm(1-pval.i,0,1))/sqrt(length(pval.i))
    n.gene[i] <- length(pval.i)
    p.gs[i] <- pnorm(gs.z.i,0,1,lower.tail = FALSE)
    z.gs[i] <- gs.z.i
  }
  gs.df <- data.frame("Path.Name" = gs.info$gs.used,"N.genes" = n.gene,"Z.stat" = z.gs, "Pval" = p.gs)
  p.ord <- order(gs.df$Pval)
  gs.df <- gs.df[p.ord,]
  return(list(res.df = gs.df, n.eff = n.eff,gs.Dmat = gs.info$D.mat[p.ord,]))
}


gene.res.extract <- function(p.gene,gs.out,which.path){
  #Get row
  if(is.numeric(which.path)){
    which.row <- which.path}else
      if(is.character(which.path)){
        which.row <- which(rownames(gs.out$gs.Dmat)==which.path)
      }
  which.col <- which(gs.out$gs.Dmat[which.row,]!=0)
  which.genes <- colnames(gs.out$gs.Dmat)[which.col]
  p.gene.sub <- p.gene[which(p.gene[,1]%in%which.genes),]
  return(p.gene.sub)
}


