setwd('G:\\My Drive\\James_files\\talkr')
library(talklr)
library(dplyr)

library(talklr)

expressed_thresh=4
pseudo_count=1

expression_data<-read.table(system.file("extdata", "glom_normal_data.txt", package = "talklr"),header=T,sep="\t",stringsAsFactors =F)


  expression_data$genes<-toupper(expression_data$genes)
  n_cell<-ncol(expression_data)-1
  good<- rowSums(expression_data[,2:ncol(expression_data)]>expressed_thresh)>=1 # expressed above threshold in at least one cell type
  expressed_genes<-expression_data$genes[good]
  expression_data<-expression_data[good,]
  expression_data[,2:ncol(expression_data)]<- expression_data[,2:ncol(expression_data)] + pseudo_count

  expressed_net<- receptor_ligand[receptor_ligand$Ligand.ApprovedSymbol %in% expressed_genes & receptor_ligand$Receptor.ApprovedSymbol %in% expressed_genes,]
  ind<-match(expressed_net$Ligand.ApprovedSymbol,expression_data$genes)
  expressed_net<-cbind(expressed_net,expression_data[ind,2:ncol(expression_data)])
  colnames(expressed_net)[17:(17+n_cell-1)]<-paste("ligand_",colnames(expression_data)[2:ncol(expression_data)],sep="")

  ind<-match(expressed_net$Receptor.ApprovedSymbol,expression_data$genes)
  expressed_net<-cbind(expressed_net,expression_data[ind,2:ncol(expression_data)])
  colnames(expressed_net)[(17+n_cell):(17+2*n_cell-1)]<-paste("receptor_",colnames(expression_data)[2:ncol(expression_data)],sep="")

#  expressed_net$KL<-pairwise.interaction.KL(as.matrix(expressed_net[,17:(17+n_cell-1)]),as.matrix(expressed_net[,(17+n_cell):(17+2*n_cell-1)]),method=KL_method)

ligand_mat<-as.matrix(expressed_net[,17:(17+n_cell-1)])
receptor_mat<-as.matrix(expressed_net[,(17+n_cell):(17+2*n_cell-1)])

  n<-ncol(ligand_mat)
  p<-nrow(ligand_mat)
  ligand_mat<-ligand_mat
  receptor_mat<-receptor_mat
  interaction_mat<-matrix(0,nrow=nrow(ligand_mat),ncol=1)

#      for (i in 1:n){
#      for (j in 1:n){
#        interaction_mat<-cbind(interaction_mat,apply(cbind(ligand_mat[,i],receptor_mat[,j]),1,min))
#      }
#    }

    for (i in 1:n){
      temp_mat<-matrix(rep(ligand_mat[,i],n),nrow=p,ncol=n,byrow=F)
      interaction_mat<-cbind(interaction_mat,temp_mat*receptor_mat)
    }    

  interaction_mat<-interaction_mat[,-1]
  interaction_mat<- interaction_mat/matrix(rep(rowSums(interaction_mat),n^2),nrow=p,ncol=n^2,byrow=F)
  interaction_mat<-interaction_mat*log2(interaction_mat*(n^2))
  interaction_mat[is.nan(interaction_mat)]<-0


ligand_exprs<-as.numeric(expressed_net[1,17:19])
receptor_exprs<-as.numeric(expressed_net[1,20:22])

cell_labels<-c("podo","mesa","endo")
thresh=1

  n_cell<-length(cell_labels)
  norm_mat<-matrix(0,n_cell,n_cell)
  ligand_exprs[ligand_exprs<thresh]<- 0 #if not expressed, set to 0
  receptor_exprs[receptor_exprs<thresh]<- 0

  for (i in 1:n_cell){
    norm_mat[i,]<-ligand_exprs[i]*receptor_exprs
  }
  # min_val<-min(norm_mat)
  # max_val<-max(norm_mat)
  # final_mat<-(norm_mat - min_val)/(max_val-min_val)
  total<-sum(norm_mat)
  final_mat<-norm_mat/total
  