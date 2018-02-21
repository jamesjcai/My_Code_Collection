
##########################################################################################
# Adaptive Combination of P-values (ADA) Algorithm for Case-Control Sequence Data
# To Pinpoint Rare Causal Variants from A Large Number of Variants in A Gene
##########################################################################################
# If you use this code to analyze data, please cite the following three papers:   
# Lin W-Y (2016). Beyond rare-variant association testing: Pinpointing rare causal variants in case-control sequencing study. Scientific Reports, 6: 21824.
# Lin W-Y, Lou XY, Gao G, Liu N (2014). Rare Variant Association Testing by Adaptive Combination of P-values. PLoS ONE, 9: e85728. [PMID: 24454922]. 
# Cheung YH, Wang G, Leal SM, Wang S (2012). A fast and noise-resilient approach to detect rare-variant associations with deep sequencing data for complex disorders. Genetic epidemiology 2012, 36(7):675-685.
# Any question, please contact: Wan-Yu Lin, linwy@ntu.edu.tw, Institute of Epidemiology and Preventive Medicine, National Taiwan University 
# Thank you.
# Acknowledgements: This program was developed based on Cheung et al.'s code. Thank them for sharing their R code.
##########################################################################################


suppressMessages(library('mgcv'));   

ADATest <- function(file, mafThr = 0.05, num_perm = 1000,
        midp = TRUE, mode = 'additive', twoSided = TRUE) {
    data = as.matrix(read.table(file)); phen = data[, 1];
    cases = data[phen == 1, -1];
    ctrls = data[phen == 0, -1];

    n_alleles = colSums(rbind(cases * (cases >= 0), ctrls * (ctrls >= 0)));
    n = colSums(rbind(cases >= 0, ctrls >= 0));
    idx1 = which(n_alleles > 0 & n_alleles <= (n * 2 * mafThr));
    cases = cases[, idx1]; ctrls = ctrls[, idx1];
    
    result = getPvalADA(cases, ctrls, num_perm, midp, mode, twoSided);
    p1 = result[[1]];
    optimal.t = result$optimal_threshold
    posit = idx1[result$posit]
    method = (if (midp) 'Sigma-MidP' else 'Sigma-P');
    return(list(pval=p1, optimal.t=optimal.t, posit=posit));
}


getPvalADA <- function(cases, ctrls, num_perm = 1000,
        midp = TRUE, mode = 'additive', twoSided = TRUE) {
    n_cases = nrow(cases);
    n_ctrls = nrow(ctrls);
    n = n_cases + n_ctrls;
    casesori <- cases
    ctrlsori <- ctrls

    s1 = c();
    
    for (i in 1: (num_perm + 1)) {
  res <- c(getScoreADA(cases, ctrls, midp, mode, twoSided));
	s1 = rbind(s1, res[[1]]);
	data = rbind(cases, ctrls);
	data = data[sample(1: n, n), ];
	cases = data[1: n_cases, ];
	ctrls = data[(n_cases + 1): n, ];
    }


    can.threshold <- seq(0.10,0.20,0.01)
    p1 <- c()
    p1.sub <- c()
    for(kkk in 1:length(can.threshold)){
      num_perm1 <- num_perm-sum(is.na(s1[,kkk]))
      ns1 <- s1[is.na(s1[,kkk])==0,kkk]
      p1.sub[1] <- (sum(ns1[2: (num_perm1 + 1)] >= ns1[1])+1) / (num_perm1 + 1)
      for(lll in 2:nrow(s1)){
        p1.sub[lll] <- (sum(ns1[2: (num_perm1 + 1)] >= ns1[lll])) / (num_perm1)
      }
      p1 = cbind(p1,p1.sub);
    }

    p1.min <- c()
    for(lll in 1:nrow(s1)){
      p1.min[lll] <- min(p1[lll,])
    }
    num_perm.min <- num_perm-sum(is.na(p1.min))
    nsp1 <- p1.min[is.na(p1.min)==0]
    pval1 = (sum(nsp1[2: (num_perm.min + 1)] <= p1.min[1]) + 1) / (num_perm.min + 1); 
    optimal_threshold <- mean(can.threshold[which(p1[1,]<=p1.min[1])]) 
    posit <- getScoreADAoriginal(casesori, ctrlsori, midp = FALSE, mode = 'additive', twoSided = TRUE, optimal_threshold)
        
    return(list(pval1 = pval1, optimal_threshold = optimal_threshold, posit=posit));  
}


getScoreADA <- function(cases, ctrls, midp = FALSE, mode = 'additive', twoSided = TRUE) {

    hasMiss = FALSE;
    if (any(ctrls < 0) || any(cases < 0)) {
        hasMiss = TRUE;
    }


    casesValidGen = (cases >= 0);
    ctrlsValidGen = (ctrls >= 0);

    c = 2;
    switch(tolower(mode),
        dominant = {
            ctrls[ctrlsValidGen] = (ctrls[ctrlsValidGen] >= 1);
            cases[casesValidGen] = (cases[casesValidGen] >= 1);
            c = 1;
        },                                                                  
        recessive = {
            ctrls[ctrlsValidGen] = (ctrls[ctrlsValidGen] == 2);
            cases[casesValidGen] = (cases[casesValidGen] == 2);
            c = 1;
        }
    )

    ku = t(colSums(ctrls * ctrlsValidGen));
    ka = t(colSums(cases * casesValidGen));

    n_s = 1;
    if (twoSided) {
        n_s = 2;
    }

    score = c();
    ww = c();
    can.threshold <- seq(0.10,0.20,0.01)
    score2v <- matrix(NA,2,length(can.threshold))
    for (j in 1: (if (twoSided) 2 else 1)) {
        if (j == 1) {
            ptr = ((ka/nrow(cases)) > (ku/nrow(ctrls)));
            k1 = ku[ptr]; k2 = ka[ptr];
        } else {
            ptr = ((ku/nrow(ctrls)) > (ka/nrow(cases)));
            k1 = ka[ptr]; k2 = ku[ptr];
        }

        if (all(ptr == 0)) {
            score[j] = 0;
            next;
        }

        if (hasMiss) {
            if (j == 1) {
                m = colSums(ctrlsValidGen) * c;
                k = colSums(casesValidGen) * c;
            } else {
                m = colSums(casesValidGen) * c;
                k = colSums(ctrlsValidGen) * c;
            }
            m = t(m[ptr]); k = t(k[ptr]);
        } else {
            if (j == 1) {
                m = nrow(ctrls) * c; k = nrow(cases) * c;
            } else {
                m = nrow(cases) * c; k = nrow(ctrls) * c;
            }
        }


        if (hasMiss) {
            gp = uniquecombs(cbind(k1, k2, m, k));
            m = gp[, 3]; k = gp[, 4];
        } else {
            gp = uniquecombs(cbind(k1, k2));
        }
        idx = attr(gp, 'index');
        
        q = (gp[, 1] + gp[, 2] + 1) / (m + k + 2);  
        
        

        w1 = dbeta(q, 1, 25)
               
        
                
        n_alleles = gp[, 1] + gp[, 2];
        if (midp) {
            w2 = phyper(gp[, 1] - 1, m, k, n_alleles) + 
                    dhyper(gp[, 1], m, k, n_alleles) * 0.5;
        } else {
            w2 = phyper(gp[, 1], m, k, n_alleles);
        }
        w2log = -log(w2);
        
        n = 0;
        if (length(gp) > 0) {
            for (i in 1: nrow(gp)) {
                n[i] = sum(idx == i);
            }
        }
        
        wwv <- matrix(NA,length(w1),length(can.threshold))
        for(threshold.i in 1:length(can.threshold)){
          wwv[,threshold.i] = w1 * w2log * (w2<can.threshold[threshold.i]);
          score2v[j,threshold.i] <- sum(wwv[,threshold.i] * n)
        }
    }
    max.score2v <- apply(score2v,2,max)
    score = list(max.score2v);
}


getScoreADAoriginal <- function(cases, ctrls, midp = FALSE, mode = 'additive', twoSided = TRUE, optimal_threshold) {

    hasMiss = FALSE;
    if (any(ctrls < 0) || any(cases < 0)) {
        hasMiss = TRUE;
    }


    casesValidGen = (cases >= 0);
    ctrlsValidGen = (ctrls >= 0);

    c = 2;
    switch(tolower(mode),
        dominant = {
            ctrls[ctrlsValidGen] = (ctrls[ctrlsValidGen] >= 1);
            cases[casesValidGen] = (cases[casesValidGen] >= 1);
            c = 1;
        },                                                                  
        recessive = {
            ctrls[ctrlsValidGen] = (ctrls[ctrlsValidGen] == 2);
            cases[casesValidGen] = (cases[casesValidGen] == 2);
            c = 1;
        }
    )

    ku = t(colSums(ctrls * ctrlsValidGen));
    ka = t(colSums(cases * casesValidGen));

    n_s = 1;
    if (twoSided) {
        n_s = 2;
    }

    score = c();
    ww = c();
    can.threshold <- seq(0.10,0.20,0.01)
    score2v <- matrix(NA,2,length(can.threshold))
    posi <- c()
    for (j in 1: (if (twoSided) 2 else 1)) {
        if (j == 1) {
            ptr = ((ka/nrow(cases)) > (ku/nrow(ctrls)));
            k1 = ku[ptr]; k2 = ka[ptr];
        } else {
            ptr = ((ku/nrow(ctrls)) > (ka/nrow(cases)));
            k1 = ka[ptr]; k2 = ku[ptr];
        }

        if (all(ptr == 0)) {
            score[j] = 0;
            next;
        }

        if (hasMiss) {
            if (j == 1) {
                m = colSums(ctrlsValidGen) * c;
                k = colSums(casesValidGen) * c;
            } else {
                m = colSums(casesValidGen) * c;
                k = colSums(ctrlsValidGen) * c;
            }
            m = t(m[ptr]); k = t(k[ptr]);
        } else {
            if (j == 1) {
                m = nrow(ctrls) * c; k = nrow(cases) * c;
            } else {
                m = nrow(cases) * c; k = nrow(ctrls) * c;
            }
        }


        if (hasMiss) {
            gp = uniquecombs(cbind(k1, k2, m, k));
            m = gp[, 3]; k = gp[, 4];
        } else {
            gp = uniquecombs(cbind(k1, k2));
        }
        idx = attr(gp, 'index');
        
        q = (gp[, 1] + gp[, 2] + 1) / (m + k + 2);  
        
        w1 = dbeta(q, 1, 25)
               
                        
        n_alleles = gp[, 1] + gp[, 2];
        if (midp) {
            w2 = phyper(gp[, 1] - 1, m, k, n_alleles) + 
                    dhyper(gp[, 1], m, k, n_alleles) * 0.5;
        } else {
            w2 = phyper(gp[, 1], m, k, n_alleles);
        }
        whichtype <- which(w2<optimal_threshold)
        for(lin in 1:length(whichtype)){
          posi <- c(posi,which(ptr==1)[which(idx==whichtype[lin])])
        }
    }
    return(posi)
}

