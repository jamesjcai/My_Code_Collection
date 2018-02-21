

######################################################
######################################################

########### Convert random indices to indices from the big ped file

SmallIndToBig<-function(ro, indzJ, fpv ){
    File_PValue <- fpv
    small_curr_set_logfile<-paste(File_PValue, '_', indzJ, '.ro', ro, '.curr_set.ro', ro, sep='')
    rand_indx_logfile<-paste(file=File_PValue, '_', indzJ,'.random_idx', sep='')
    
    
    if (file.exists(small_curr_set_logfile) & file.exists(rand_indx_logfile)){
        cat(paste(small_curr_set_logfile, "\n"))
        cat(paste(rand_indx_logfile, "\n"))
        small_curr_set<-scan(file=small_curr_set_logfile)
        rand_indx<-scan(file=rand_indx_logfile)
        indx_from_big<-rand_indx[small_curr_set]
        big_curr_set_logfile<-paste(File_PValue,'.ro',ro,'.big_curr_set.ro',ro, sep='')
        write(indx_from_big, file=big_curr_set_logfile, ncolumns=length(indx_from_big),append=TRUE)
        big_set20_logfile<-paste(File_PValue, '.ro',ro,'.big_set20.ro',ro, sep='')
        write(rand_indx, file=big_set20_logfile, ncolumns=length(rand_indx),append=TRUE)
        
        #unlink(x=small_curr_set_logfile, recursive = FALSE)
        #unlink(x=rand_indx_logfile, recursive = FALSE)
        
    }
    
    else{
        cat(paste('No file XX ', small_curr_set_logfile, "\n"))
        cat(paste('No File XX ', rand_indx_logfile, "\n"))
        
    }
    
}

######################################################
######################################################

########### SKAT backwards elimination with resampling

skat.be <- function( Z, y, w, File.Out, ro=0, N.SIMR=300, basedir=getwd() ) {
    
    basename=File.Out
    curr_home_dir=paste(basedir, '/', basename, sep='')
    
    if (file.exists( curr_home_dir )){
        cat (paste("Directory", curr_home_dir, "exists\n"))
    } else {
        dir.create(file.path(basedir, basename))
        cat (paste("Created", curr_home_dir, '\n'))
    }
    
    File_PValue=paste(curr_home_dir, '/', File.Out, sep='')
    
    ## Load ped
    
    ## Get genotypes
    wfunc=diag(w);
    Z.all=Z%*%wfunc;
    cat ("Z.all", dim(Z.all), "\n" )
    
    scatit(Z.all,  y, File.Out=paste(File_PValue, '_', 'BIG', ".ro",ro, sep=''), ro=ro)
    
    
    for(j in 1:N.SIMR){
        
        rand.ind<-sort(sample(1:dim(Z.all)[2], 10, replace=F))
        write(rand.ind, file=paste(File_PValue,'_', j, ".random_idx", sep=''), ncolumns=length(rand.ind))
        Z.rand<-Z.all[,rand.ind]
        ###
        cat(paste("Processing random ", j, 'of ', 'Big', '\n'))
        scatit(Z=Z.rand, y, File.Out=paste(File_PValue,  '_', j, ".ro",ro, sep=''), ro=ro)
        
        
        SmallIndToBig(ro=ro, indzJ=j, File_PValue)
        
        unlink(x=paste(File_PValue,  '_', j, ".ro",ro, sep=''), recursive = FALSE)
        
        
        
        rand_indx_logfile<-paste(file=File_PValue, '_', j,'.random_idx', sep='')
        unlink(x=rand_indx_logfile, recursive = FALSE)
        
        small_curr_set_logfile<-paste(File_PValue, '_', j, '.ro', ro, '.curr_set.ro', ro, sep='')
        unlink(x=small_curr_set_logfile, recursive = FALSE)
        
    }
    
}

skat.cl <- function( File.Ped, File.W, File.Out, ro=0, N.SIMR=300, basedir=getwd() ) {
    
    basename=File.Out
    curr_home_dir=paste(basedir, '/', basename, sep='')
    
    if (file.exists( curr_home_dir )){
        cat (paste("Directory", curr_home_dir, "exists\n"))
    } else {
        dir.create(file.path(basedir, basename))
        cat (paste("Created", curr_home_dir, '\n'))
    }
    
    File.Ped=paste(basedir, '/', File.Ped, sep='')
    File.W=paste(basedir, '/', File.W, sep='')
    File_PValue=paste(curr_home_dir, '/', File.Out, sep='')
    cat ( File_PValue, "\n" )
    
    ## Load ped
    BigData<-read.table(File.Ped, header=FALSE)
    cat ("BigData", dim(BigData), "\n" )

    ## Get genotypes
    Z.org<-as.matrix(BigData[,-c(1:6)])
    wf=scan(File.W);
    wfunc=diag(wf);
    Z.all=Z.org%*%wfunc;
    cat ("Z.all", dim(Z.all), "\n" )
    ## Case-Ctrl information
    y<-BigData[,6]-1
    
    
    scatit(Z=Z.all,  y, File.Out=paste(File_PValue, '_', 'BIG', ".ro",ro, sep=''), ro=ro)
    
    
    for(j in 1:N.SIMR){
        
        rand.ind<-sort(sample(1:dim(Z.all)[2], 10, replace=F))
        write(rand.ind, file=paste(File_PValue,'_', j, ".random_idx", sep=''), ncolumns=length(rand.ind))
        Z.rand<-Z.all[,rand.ind]
        ###
        cat(paste("Processing random ", j, 'of ', 'Big', '\n'))
        scatit(Z=Z.rand, y, File.Out=paste(File_PValue,  '_', j, ".ro",ro, sep=''), ro=ro)
        
        
        SmallIndToBig(ro=ro, indzJ=j, File_PValue)
        
        unlink(x=paste(File_PValue,  '_', j, ".ro",ro, sep=''), recursive = FALSE)
        
        
        
        rand_indx_logfile<-paste(file=File_PValue, '_', j,'.random_idx', sep='')
        unlink(x=rand_indx_logfile, recursive = FALSE)
        
        small_curr_set_logfile<-paste(File_PValue, '_', j, '.ro', ro, '.curr_set.ro', ro, sep='')
        unlink(x=small_curr_set_logfile, recursive = FALSE)
        
    }

}