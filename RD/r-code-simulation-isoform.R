      name <- c("A","B","C")
      lexon <- c(400,200,150,400)
      exon_c <- rep(0,length(lexon))
      INDEX  <- matrix(c( 1,1,1,  1,0,1,  1,1,0,  1,1,1),nrow=length(lexon),ncol=length(name),byrow=T)

      DIST   <- matrix(0,nrow=length(lexon),ncol=length(name))
      NM     <- matrix(0,nrow=length(lexon),ncol=length(name))
      CCC <- 1  
      AAA <- 0.001
      TTT <- c(0.5,0.2,0.3)    
      for(k in 1:length(lexon))
      {
          if(k==1)
          {DIST[k,] <- lexon[k] * INDEX[k,]/2  }
          if(k>1)    
          {DIST[k,] <- INDEX[k,] * (lexon[1:k]%*% INDEX[1:k,] - lexon[k] * INDEX[k,]/2)}
      }
      EEE   <- INDEX *  exp(-AAA*DIST)
      for(k in 1:length(name))       
      {
        noise <- exp(rnorm(length(lexon),0,0.1))
        NM[,k] <- ( CCC * noise  * TTT[k] * EEE[,k])* lexon
      }
       exon_c <- apply(NM,1,sum)    
 

      NM     <- matrix(0,nrow=length(lexon),ncol=length(name)) 
      # ltranscript <- rep(0,length(name))
      ### EM  
      alpha  <- 0
      theta  <- rep(1/length(name),length(name))
      for(iter in 1:1000)  
      {     ## E-Step
            EEE  <-  exp(-alpha*DIST) 
            for(k in 1:length(lexon))
            {  eeek <- sum(INDEX[k,]*theta*EEE[k,])
              # if(eeek> 1e-12)
               {NM[k,] <- exon_c[k]*(INDEX[k,]*theta*EEE[k,])/eeek}
              # else
              # {NM[k,] <- rep(0,length(name))}
            }
            ## M-Step
                a <- alpha
                for(n in 1:10)
                {
                    eee  <-  exp(-a*DIST)
                    D0 <- lexon %*% (INDEX*eee) 
                    D1 <- lexon %*% (INDEX*eee*DIST)
                    D2 <- lexon %*% (INDEX*eee*DIST^2)
                    TI <- apply(INDEX*NM,2,sum)
                    First <-  -sum(INDEX*NM*DIST)  
                    Second <- 0   
                    for(ttt in 1:ncol(NM))
                    {
                       First  <- First+TI[ttt]*(D1[ttt]/D0[ttt]) 
                       Second <- Second +TI[ttt]*((D1[ttt]^2-D0[ttt]*D2[ttt])/D0[ttt]^2)            
                    }
                    a <- a - First/Second  
                }   
                alpha <- a 
                EEE <- exp(-alpha*DIST)
                WWW <- lexon %*% (INDEX*EEE)    
                BBB <- (apply(INDEX*NM,2,sum)/sum(exon_c))/WWW
                theta <- BBB/sum(BBB)
      } ## EM Iteration
      print(round(c(1,alpha,theta),4))
      
      NM     <- matrix(0,nrow=length(lexon),ncol=length(name))
      # ltranscript <- rep(0,length(name))
      ### EM
      alpha  <- 0
      theta  <- rep(1/length(name),length(name))
      for(iter in 1:1000)
      {     ## E-Step
            EEE  <-  exp(-alpha*DIST)
            for(k in 1:length(lexon))
            {  eeek <- sum(INDEX[k,]*theta*EEE[k,])
              # if(eeek> 1e-12)
               {NM[k,] <- exon_c[k]*(INDEX[k,]*theta*EEE[k,])/eeek}
              # else
              # {NM[k,] <- rep(0,length(name))}
            }
            ## M-Step
                EEE <- exp(-alpha*DIST)
                WWW <- lexon %*% (INDEX*EEE)
                BBB <- (apply(INDEX*NM,2,sum)/sum(exon_c))/WWW
                theta <- BBB/sum(BBB)
      } ## EM Iteration
      print(round(c(2,alpha,theta),4))
      


