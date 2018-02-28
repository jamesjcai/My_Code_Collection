


%macro random(nvarX, f);
 
        random  %do i = 1 %to &nvarX*&f;
         X&i
  %end; / g type=toep(1);
%mend random;

%macro model_XZ(Y, f);
 
        model &Y=  %do i = 1 %to &f;
         XZ&i
  %end; 
%mend random;

%macro random_tau(nvarX, f, taudata);
 
        random  %do i = 1 %to &nvarX;
         X&i
  %end; /g gdata=&taudata;
%mend random_tau;


%macro numvar(number,tdata);
  %local i;
  %do i = 1 %to &number;
    &tdata&i
  %end;
%mend numvar;

options nomprint nosymbolgen nomlogic;
%macro HM_approach(dir_input=, dir_output=, Ydata=, Xdata=, Zdata=, IGprior=, nvarW=,  nvarZ=, maf_rare= );

data pdata;	set IG.&IGprior;run;

filename results "&dir_output.out_results_vps13b.txt";
filename logfile "&dir_output.log_results_vps13b.txt";



%global tau_g posterior_mean;
proc printto print=results log=logfile new;run;

 proc iml;

 nerror_delta=0;nblowout=0;n_tau_g=0;nstderr_neg=0;
 
	%let p1 = %eval(&nvarW + 2);
    %let p2 = %eval(&nvarZ + &nvarW + 1);
    %let d1 = %eval(&nvarZ + &p1);
	%let p1m1=%eval(&p1-1);
	%let nvarZ_nointercept=%eval(&nvarZ-1);

	
store _all_;

   
   data bioinf_allvar;
   infile "&dir_input&Zdata..txt" lrecl=1000 truncover dlm="09"x dsd;
   input poly gerp func $;
   if func in ('NON_SYNO','SPLICE_S', 'STOP_GAI') then funct=1;
     else funct=0;
   run;


   data case_contr;
   infile "&dir_input&Ydata..txt" lrecl=1000 truncover dlm="09"x dsd;
   input case;
   run;

proc iml;load _all_;	

    use bioinf_allvar; 
    read all var _num_ into bioinf_allvar;
	
	nvars_f = nrow(bioinf_allvar);
    call symput('nvarsf',left(char(nvars_f)));print &nvarsf ;print "Number of variants:" &nvarsf; 

	beta_estimate=j(&nvarsf,1,.);stderr_est=j(&nvarsf,1,.);
    zvalues=j(&nvarsf,1,.);

store _all_;

  data Xdata_allvar0;
   infile "&dir_input&Xdata..txt" LRECL= 1000 truncover dlm="09"x dsd;
     input %numvar(&nvarsf, X);
	
   run;

   data Xdata_allvar;
   infile "&dir_input&Xdata..txt" LRECL= 1000 truncover dlm="09"x dsd;
     input %numvar(&nvarsf, X);
	 array v X1--X&nvarsf;
 do over v; v=min(v,1);end;
   run;

proc iml;load _all_;	

use case_contr; 
    read all var _num_ into Y;

use Xdata_allvar0; 
    read all var _num_ into X0;

	maf=X0[,][:,]/2;

use Xdata_allvar; 
    read all var _num_ into X1;

	NS = nrow(X1);
    call symput('N',left(char(NS)));

	var_rare=maf<=&maf_rare & maf>0& t(bioinf_allvar[,3])=1;/*third column of bioinf_allvar is the nonsynonymous status*/
	X=X1[,loc(maf<=&maf_rare & maf>0& t(bioinf_allvar[,3])=1)];
	/*restricting only to nonsynonymous variants that have MAF positive and smaller 
       that the threshold &maf_rare*/; 

 nvarrare = sum(var_rare);
 call symput('nvar_rare',left(char(nvarrare)));/*number of rare variants remaining in the model*/
	print "Number of rare variants with MAF below the specified threshold:" &nvar_rare;


    cases=X1[loc(Y=1),][+,];
    controls=X1[loc(Y=0),][+,];

	bio_no0=bioinf_allvar[loc(maf<=&maf_rare & maf>0& t(bioinf_allvar[,3])=1),1:2];
    Z=j(&nvar_rare,1,1)||bio_no0;/*first and second columns of bioinf_allvar are polyphen and gerp*/

	print "Higher level covariates"; print Z;

	/*create XZ matrix; */
    XZ=X*Z;cnameXZ='XZ1':"XZ&nvarZ.";

	cnameY='Y';cnameX='X1':"X&nvar_rare.";

	X_W_Y_XZ=Y||X||XZ;cname=cnameY||cnameX||cnameXZ;

	create dataset from X_W_Y_XZ[colname=cname];
     append from X_W_Y_XZ;

store _all_;

%glimmix(data=dataset, procopt=mmeqsol,
  stmts=%model_XZ(Y,&nvarZ);
				   %random(&nvar_rare,1);
    %str( make 'mmeqsol' out=mmeqsol;make 'g' out=g;
        ),
        options=noprint, error=binomial
             );
      

data covariance;set _cov;
if covparm ne 'Variance' then delete; else call symput('tau_g', left(estimate));
run;

%macro dummy; 
%do i=1 %to &nvar_rare;
if X&i=1 then variant=&i; 
%end;  
%mend dummy;

data ds1;set _ds;
%dummy;
run; 

ods listing close;
proc mixed data=ds1;
class variant;
  %model_XZ(_z,&nvarZ);                             
   random   variant;  
   weight _w;  
   prior data=pdata/out=post alg=rwc nsample=10000; 	
   ods output covparms=g1;
run;   

proc univariate data=post round=0.01;
var covp1;
ods output basicmeasures=bayesian;
run;

ods listing;


data bayesian;set bayesian;
if locmeasure ne 'Mean' then delete; else call symput('posterior_mean', left(locvalue));
run;

data tau_posterior;        
        %do i=1 %to &nvar_rare;
        row = &i;
        col = &i;value=&posterior_mean;output;%end;keep row col value;
      run;

%glimmix(data=dataset, procopt=mmeqsol,
  stmts=%model_XZ(Y,&nvarZ);
				   %random_tau(&nvar_rare,1,tau_posterior);
    %str( make 'mmeqsol' out=mmeqsol;
        ),
        options=noprint, error=binomial
             );

proc IML; load _all_;

    %let d2 = %eval(&nvarZ + &nvarW + &nvar_rare + 1);
	%let nmsol = %eval(&d2+1); 

    use mmeqsol;
    read all var { %numvar(&nmsol, col) } into mmeqsol;/*print mmeqsol;*/
	/*read all into mmeqsol;print mmeqsol;*/

	varpi = mmeqsol[&p1:&p2, &p1:&p2];/* Covariance of second stage coefficient estimates */
    varD = mmeqsol[&d1:&d2, &d1:&d2];/* Covariance of random coefficient estimates*/
    cov_pd = mmeqsol[&p1:&p2, &d1:&d2];/* Covariance of second stage and random coeff.     */
    pi = mmeqsol[&p1:&p2, &nmsol];/* Second stage coefficient estimates (i.e., of XZ) *//*print pi;*/
    delta = mmeqsol[&d1:&d2, &nmsol];/* Random coefficient estimates		     */

    b=Z*pi + delta; /* calculating the estimates for beta  */

    varb=Z*varpi*t(Z) + vard + Z*cov_pd + t(cov_pd)*t(Z);
    /* Covariance matrix for the betas, equation (5) 	     */

    stderr=vecdiag(varb);/*print stderr;*/
    do i=1 to nrow(varb);
    if stderr[i]<0 then do; stderr[i]=0;nstderr_neg=nstderr_neg+1;end;
    end;if max(stderr)>1000 then do;nblowout=nblowout+1;end;
    stderr=sqrt(stderr);
    lower=b-1.96*stderr; /* upper and lower confidence bounds of beta */
    upper=b+1.96*stderr;/*print upper;*/

    varpi=vecdiag(varpi);/*print varpi;*/
    /* variances of the second stage covariate estimates */
	do i=1 to nrow(varpi);
    if varpi[i]<0 then varpi[i]=0;
    end;
    stdpi=sqrt(varpi);
	lower_pi = pi-1.96*stdpi; /* upper and lower confidence bounds of pi */
    upper_pi = pi + 1.96*stdpi;/*print upper_pi;*/

 beta_estimate[loc(maf<=&maf_rare & maf>0& t(bioinf_allvar[,3])=1)]=b;
    stderr_est[loc(maf<=&maf_rare & maf>0& t(bioinf_allvar[,3])=1)]=stderr;
   zvalues[loc(maf<=&maf_rare & maf>0& t(bioinf_allvar[,3])=1)]=b/stderr;

	cnameBetaEst='beta_est';
	cnameStderr='stderr';
    cnameZval='zvalues';
	cnameBioinf='bioinf1':"bioinf&nvarZ_nointercept.";
	cnameCase='cases';
	cnameContr='controls';
	  
    cname_caus=cnameZval||cnameBetaEst||cnameStderr||cnameBioinf||cnameCase||cnameContr;
 
Causal_Zval=Zvalues||beta_estimate||stderr_est||bioinf_allvar[,1:2]||t(cases)||t(controls);

	create Rdataset from Causal_Zval[colname=cname_caus];
     append from Causal_Zval;

	store _all_;


data _null_;
  set Rdataset;
  FILE  "&dir_output.HM_results_vps13b.txt" ; 
  /*names in the put statement should be the names as they appear in the cname_caus*/;
  put zvalues beta_est stderr %numvar(&nvarZ_nointercept,bioinf) cases controls;
run;


quit;quit;

  proc printto;run;
%mend HM_approach;


%inc 'G:\GEM\iuliana\vps13b\corrected\code for journal\glimmix.sas' / nosource;
libname IG 'G:\GEM\iuliana\vps13b\corrected\code for journal\';

%HM_approach(dir_input=,dir_output=,
Ydata=vps13b_case_contr, Xdata=vps13b_ped, Zdata=vps13b_PolyGerp, IGprior=igparameters01,nvarW=0, 
nvarZ=3, maf_rare=0.05); 

/*the macro outputs a dataset "vps13b_results.txt" which lists the zvalues, the estimated
log relative risks of the variants, their standard errors, followed by the higher level covariates
used in the model (in this case an intercept, Polyphen, and Gerp) as well as the corresponding 
case control frequencies for each variant*/