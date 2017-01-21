#include "p2sweep.h"
#include "hualib.h"
#define maf 1

void getCount(int n1, int n2,int rawData[][10000],double* f1, double *f2, int *x1, int *x2,double *position, double *rawPosition, int nsegsites, int *nnsegsites,double rho) 
{
  int i,j;
  int k;
  int tempCount1,tempCount2;
  extern int N;
  // extern int rawData[1000][1000];
  //extern int **data;

  for(i=0,k=0;i<nsegsites; i++)
    {
      tempCount1=0;
      tempCount2=0;
      for(j=0;j<n1;j++)
	{
	  tempCount1=tempCount1+rawData[j][i];
	  //  data[j][k]=rawData[j][i];
	}
      for(j=n1;j<n1+n2;j++){
	tempCount2+=rawData[j][i];
	//	data[j][k]=rawData[j][i];
      }
      if(tempCount2>maf && tempCount2<n2-maf/* && tempCount1>0 && tempCount1<n1*/)
	{
	  f1[k]=tempCount1/((double)n1);
	  f2[k]=tempCount2/((double)n2);
	  x1[k]=tempCount1;
	  x2[k]=tempCount2;
	  position[k]=(rawPosition[i]-rawPosition[0])*(rho/(4.0*N));
	  //  printf("%d %f %f %f\n", i,rawPosition[i],rawPosition[0],position[k]);
	  k++;
	}      
    }
  *nnsegsites=k;
  // printf("nnsegsites: %d\n",*nnsegsites);
  // for(i=0;i<*nnsegsites;i++)
  //    printf("%d %f %f %f %f %d %d\n",i,position[i],rho,q1[i],q2[i],x1[i],x2[i]);
}


void getCount2(int n1, int n2,double* rawq1, double* rawq2, double* q1, double *q2, int *x1, int *x2,double *rawposition,double *position,int nsegsites, int *nnsegsites) 
{
  int i,j;
  int k;
  int tempID;
  double tempF[3];

  for(i=0,k=0;i<nsegsites; i++)
    {
      if(rawq2[i]>0.0 && rawq2[i]<1.0)
	{
	  q1[k]=rawq1[i];
	  q2[k]=rawq2[i];
	  x1[k]=(int)(q1[k]*n1);
	  x2[k]=(int)(q2[k]*n2);
	  position[k]=(rawposition[i]-rawposition[0]);
	  k++;
	}      
    }
  if(k<MINSITES)
    {//printf("k <50\n");
    exit(0);
    }
  else if(k==MINSITES)
    {
      *nnsegsites=k;
    }
  else
    {

      /*
      for(j=k;j>=MINSITES;j--)
	{
	  tempID=(int)(j*ran1(&idum));
	  q1[tempID]=q1[j-1];
	  q2[tempID]=q2[j-1];
	  x1[tempID]=x1[j-1];
	  x2[tempID]=x2[j-1];
	  position[tempID]=position[j-1];
	  //  printf("k %d,tempID %d,j=%d, %f %f \n",k,tempID,j, position[tempID], position[j-1]);
	}     
      *nnsegsites=j+1;
      */
      *nnsegsites=k;
      //  printf("nnsegsites %d\n",*nnsegsites);
    }
  // printf("nnsegsites: %d\n",*nnsegsites);
  //  for(i=0;i<*nnsegsites;i++)
  //    printf("%d %f %f %f %d %d\n",i,position[i],q1[i],q2[i],x1[i],x2[i]);
}


double getW(double* q2, double* q1, int nsegsites)
{
  int i;
  double tt=0.0;
  double w;
  for(i=0;i<nsegsites;i++)
    {
      tt=tt+(q1[i]-q2[i])*(q1[i]-q2[i])/(q2[i]*(1-q2[i]));
    } 
  w=tt/nsegsites;
  return w;
} 


double estimatec_s(double p0, double c, double s)
{
  int i;
  double value=0.0;
  if(s>0.0)
    {//printf("c/s*log(1/p0) %e \n",c/s*log(1/p0));
    return (1.0-pow(p0,(c/s)));
    }
  
  else 
    {return 1.0;
    }
}



double func10(double t)
{
  extern double sig2;
  extern double q2_I;
  if(t<1.0 && t>0.0)
    return 1.0/sqrt(2.0*PI*sig2)*exp(-pow((t-q2_I),2)/(2.0*sig2));
  else if(t==0.0)
    return gaussiancdf(0.0,q2_I,sqrt(sig2));
  else if(t==1.0)
    return (1.0-gaussiancdf(1.0,q2_I,sqrt(sig2)));
  else
    perror("mistake in function10\n");
}

double func1(double t)
{
  extern double sig2;
  extern double q2_I;
  extern int n1_I;
  extern int x1_I;   
  // printf("t=%f, binopdf=%f, sig2=%f, %f, firstterm=%f,func1= %f\n",t,binopdf(n1_I,x1_I,t),sig2, PI,1/sqrt(2*PI*sig2),1/sqrt(2*PI*sig2)*exp(-pow((t-q2_I),2)/(2*sig2))*binopdf(n1_I,x1_I,t));
  if(t<1.0 && t>0.0)
    return 1/sqrt(2.0*PI*sig2)*exp(-pow((t-q2_I),2.0)/(2.0*sig2))*binopdf(n1_I,x1_I,t);
  else if(t==0.0)
    return gaussiancdf(0.0,q2_I,sqrt(sig2))*binopdf(n1_I,x1_I,0.0);
  else if(t==1.0)
    return (1.0-gaussiancdf(1.0,q2_I,sqrt(sig2)))*binopdf(n1_I,x1_I,1.0);
  else 
    perror("mistake in func1\n");
 
}

double func1_seg2p(double t)
{
  extern double sig2;
  extern double q2_I;
  extern int n1_I;
  extern int x1_I; 
  // printf("t=%f, binopdf=%f, sig2=%f, %f, firstterm=%f,func1= %f\n",t,binopdf(n1_I,x1_I,t),sig2, PI,1/sqrt(2*PI*sig2),1/sqrt(2*PI*sig2)*exp(-pow((t-q2_I),2)/(2*sig2))*binopdf(n1_I,x1_I,t));
  return 1/sqrt(2*PI*sig2)*exp(-pow((t-q2_I),2)/(2*sig2))*binopdf(n1_I,x1_I,t);
 }


double func2(double t)
      {
	double value;
	extern double cc_I;
	extern double sig2;
	extern double q2_I;

	if(t>0.0 && t<1.0)
	value=(GEi(t,1.0-cc_I)*1.0/sqrt(2*PI*sig2)*exp(-1*pow((t-1+cc_I-cc_I*q2_I),2)/(2.0*pow(cc_I,2)*sig2))*(t-1+cc_I)/pow(cc_I,2))
	  +((1-GEi(t,cc_I))*1/sqrt(2*PI*sig2)*(cc_I-t)/pow(cc_I,2)*exp(-pow((t-cc_I*q2_I),2)/(2*pow(cc_I,2)*sig2)));
	else if(t==0.0)
	  value=gaussiancdf(0.0,q2_I,sqrt(sig2));
	else if(t==1.0)
	  value=(1.0-gaussiancdf(1.0,q2_I,sqrt(sig2)));
	else 
	  perror("mistake in func2\n");
	return value;
      }

double func2_seg2p(double t)
      {
	double value;
	extern double cc_I;
	extern double sig2;
	extern double q2_I;

	value=(GEi(t,1.0-cc_I)*1.0/sqrt(2*PI*sig2)*exp(-1*pow((t-1+cc_I-cc_I*q2_I),2)/(2.0*pow(cc_I,2)*sig2))*(t-1+cc_I)/pow(cc_I,2))+((1-GEi(t,cc_I))*1/sqrt(2*PI*sig2)*(cc_I-t)/pow(cc_I,2)*exp(-pow((t-cc_I*q2_I),2)/(2*pow(cc_I,2)*sig2)));
	return value;
      }

double func3(double t)
     {
       double value1,value2,value;
       extern double cc_I;
       extern double sig2;
       extern double q2_I;
       extern int n1_I;
       extern int x1_I;
       if(t>0.0 && t<1.0)
	value1=(GEi(t,1.0-cc_I)*1/sqrt(2*PI*sig2)*exp(-pow((t-1+cc_I-cc_I*q2_I),2)/(2.0*pow(cc_I,2)*sig2))*(t-1+cc_I)/pow(cc_I,2)
		+(1-GEi(t,cc_I))*1/sqrt(2*PI*sig2)*(cc_I-t)/pow(cc_I,2)*exp(-pow((t-cc_I*q2_I),2)/(2*pow(cc_I,2)*sig2)));
       else if(t==0.0)
	 value1=gaussiancdf(0.0,q2_I,sqrt(sig2));
       else if(t==1.0)
	 value1=(1.0-gaussiancdf(1.0,q2_I,sqrt(sig2)));
       else 
	 {
	 perror("mistake in func3\n");
	 printf("%f\n",t);
	 }
	value2=binopdf(n1_I,x1_I,t);
        value=value1*value2;
	//	printf("func3, t:%f value1:%f value2:%f %d %d %f\n",t, value1,value2,n1_I, x1_I,q2_I);
	return value;
     }

double func3_seg2p(double t)
     {
       double value1,value2,value;
       extern double sig2;
       extern double cc_I;
       extern double q2_I;
       extern int x1_I;
       extern int n1_I;
	value1=(GEi(t,1.0-cc_I)*1/sqrt(2*PI*sig2)*exp(-pow((t-1+cc_I-cc_I*q2_I),2)/(2.0*pow(cc_I,2)*sig2))*(t-1+cc_I)/pow(cc_I,2) +(1-GEi(t,cc_I))*1/sqrt(2*PI*sig2)*(cc_I-t)/pow(cc_I,2)*exp(-pow((t-cc_I*q2_I),2)/(2*pow(cc_I,2)*sig2)));
	value2=binopdf(n1_I,x1_I,t);
        value=value1*value2;
	// printf("func3, t:%f value1:%f value2:%f %d %d %f\n",t, value1,value2,n1_I, x1_I,q2_I);
	return value;
     }

double logPNeutrality(double w,double p0,double *q2, int *x1,int n1,int nnsegsites)
      {
	int i,j;
	double value=0.0,tempValue,tempValue2;
	extern int n1_I;
	extern double cc_I;
	
	extern int x1_I;
	extern double q2_I;
	extern double sig2;
	n1_I=n1;
	cc_I=1.0;

	for(i=0;i<nnsegsites;i++)
	  {
	    sig2=w*q2[i]*(1-q2[i]);
	    x1_I=x1[i];
	    q2_I=q2[i];
	    //	    tempValue2=qromb(func10,0.001,0.999);
            tempValue=log(qromb(func1_seg2p, 0.001, 0.999))/*+func1(0.0)+func1(1.0)*/ /*-log(qromb(func10,0.001,0.999))*/;
	   
	    value+=tempValue;
	    //	    printf("neutrality:%d %f\n",i,trapzd(func1,0.00001,0.9999,20));
	  }
	return value;
      } 

double logPNeutrality2(double w,double p0,double *q2, int *x1,int n1,int *snpList,int SnpListLength)
      {
	int i,j;
	double value=0.0,tempValue,tempValue2;
	int Oposition;
	extern int n1_I;
	extern double cc_I;
	extern int x1_I;
	extern double q2_I;
	extern double sig2;
	n1_I=n1;
	cc_I=1.0;

	for(i=0;i<SnpListLength;i++)
	  {
	    Oposition=snpList[i];
	    sig2=w*q2[Oposition]*(1-q2[Oposition]);
	    x1_I=x1[Oposition];
	    q2_I=q2[Oposition];
       	    tempValue2=log(qromb(func10,0.001,0.999));
            tempValue=log(qromb(func1_seg2p, 0.001, 0.999))/*+func1(0.0)+func1(1.0)*/ /*-log(qromb(func10,0.001,0.999))*/;
	   
	    value+=tempValue-tempValue2;
	    //	    printf("neutrality:%d %f\n",i,trapzd(func1,0.00001,0.9999,20));
	  }
	return value;
      } 




double prob_sweep_fun(int n1, int x1,double q2,double c,double s, double p0,double w)
     {
       double value,value1,value2;
       extern double cc_I;
       extern double sig2;
       extern int n1_I;
       extern int x1_I;
       extern double q2_I;

       cc_I=estimatec_s(p0,c,s);
       sig2=w*q2*(1-q2);
       n1_I=n1;
       x1_I=x1;
       q2_I=q2;
       //     printf("cc %f sig2 %f n1 %d x1 %d\n",cc_I, sig2, n1_I, x1_I);
       double tempValue1=0.0;       
       tempValue1=qromb(func2_seg2p,0.001,0.999);
       if(tempValue1>=0)
	 value1=tempValue1;
       else
	 {
	 printf("%d %d %f %f %f %f %f\n",n1,x1,q2,c,s,p0,w);
	 value1=0.0;
	 }
       //      value1=log(qromb(func2_seg2p,0.001,0.999)/*+func2(0.0)+func2(1.0)*/);
       tempValue1=qromb(func3_seg2p,0.001,0.999);
       if(tempValue1>=0)
	 value2=tempValue1;
       else
	 {
	   value2=0.0;
	   printf("func3: %d %d %f %f %f %f %f\n",n1,x1,q2,c,s,p0,w);
	 //       value2=log(qromb(func3_seg2p,0.001,0.999)/*+func3(0.0)+func3(1.0)*/);
	 }
       value=value2-value1;
       //   printf("value1 %f value2 %f %f %f %f\n",value1,value2,qromb(func3,0.001,0.999),func3(0.0),func3(1.0)); 
       //  value=value2;
       /* if(cc_I<0.00001){
	 printf("%d %d %d %f %f\n", flag1,flag2,flag3,c,s);
              printf("cc %f, value1 %f, value2 %f\n",cc_I,value1,value2);
	      } */       
       return value;
     }

double logPSelection(double w,double p0,int *x1, double *q2, double *pos,double *sValue, int sLength,int n1, int indicator,int nnsegsites,double* logPLocation, double* maxS)

      {
	int i,j,l;
	double value=0.0;
	double tempLogP,tempValue;
	double pos0;   
	int mutPos; 
	double tempMaxP,midDouble;
	double s,repos,maxSS[MINSITES];
	int tempMaxPosition, midInt;
        double logPCurve[1000];

	for(i=0/*nnsegsites/2*/;i<nnsegsites /* nnsegsites/2+1 */;i++)
	  {
            mutPos=i;
	    pos0=pos[mutPos]+0.00001;
	    //added for debug
	    // pos0=0.0025;
	    tempMaxP=-1800;
	    maxSS[i]=0.0;
	    for(j=0;j<sLength;j++)
	      {
		s=sValue[j];
		tempLogP=0.0;
                for(l=0;l<nnsegsites;l++)
	         {
                   repos=max(fabs(pos[l]-pos0),0.00001);
		   // if(mutPos==5)
		   // printf("mupos, %d, j %d, l %d,pos0 %f, repos %f\n",mutPos,j,l, pos0, repos);
		   tempLogP+=prob_sweep_fun(n1,x1[l],q2[l],repos,s,p0,w);
		   // if(i==18)
		   //  printf("%d,%d,%d,%f\n",i,j,l,prob_sweep_fun(n1,x1[l],q2[l],repos,s,p0,w));
	         }      
 
       
		//		printf("tempLogP: %f, tempMaxP,%f\n", tempLogP, tempMaxP);
		if(tempLogP>=tempMaxP)
		  {
		   
		    tempMaxP=tempLogP;
		    maxSS[i]=s;
		  }
	      }
	    //  return tempMaxP;

	    logPCurve[i]=tempMaxP;
	    //   printf("%f\n",tempMaxP);
	  }
	tempMaxP=logPCurve[0];
	tempMaxPosition=0;
	for(i=0;i<nnsegsites;i++)
	  {
	    //  tempMaxP=GEi(tempMaxP,logPCurve[i])*tempMaxP+(1-GEi(tempMaxP,logPCurve[i]))*logPCurve[i];
	    //  tempMaxPosition=GEi(tempMaxP,logPCurve[i])*tempMaxPosition+(1-GEi(tempMaxP,logPCurve[i]))*i;
	    if(tempMaxP<=logPCurve[i])
	      {
		tempMaxP=logPCurve[i];
		*maxS=maxSS[i];
		tempMaxPosition=i;
	      }
	    else
	      continue;
	      
	  }
	value=tempMaxP;
       
	*logPLocation=pos[tempMaxPosition]+0.00001;
		printf("%f\n",*logPLocation);

		/*	for(i=0;i<nnsegsites;i++)
	  printf("%f ",pos[i]);
	  printf("\n"); */
	return value;    
  }

double logPSelection2(double w,double p0,int *x1, double *q2, double *pos,double *sValue, int sLength,int n1, int indicator,int nnsegsites,double* logPLocation, double* maxS, int *snpList,int SnpListLength,double gridsize,int gridnumber,double gridstart,double totalGD)

      {
	int i,j,l;
	double value=0.0;
	double tempLogP,tempValue;
	double pos0;   
	int mutPos; 
	double tempMaxP,midDouble;
	double s,repos,maxSS[MINSITES];
	int tempMaxPosition, midInt;
        double logPCurve[1000];
	int Oposition;
	for(i=0;i<gridnumber;i++)
	  {
	    if(1>=gridnumber)
	      pos0=(gridstart+0.5)*totalGD;
	    else
	      pos0=(i*gridsize+gridstart)*totalGD;
	    //added for debug
	    // pos0=0.0025;
	    tempMaxP=-1800;
	    maxSS[i]=0.0;
	    for(j=0;j<sLength;j++)
	      {
		s=sValue[j];
		tempLogP=0.0;
                for(l=0;l<SnpListLength;l++)
	         {
		   Oposition=snpList[l];
                   repos=max(fabs(pos[Oposition]-pos0),0.00001);
		   // if(mutPos==5)
		   // printf("mupos, %d, j %d, l %d,pos0 %f, repos %f\n",mutPos,j,l, pos0, repos);
		   tempLogP+=prob_sweep_fun(n1,x1[Oposition],q2[Oposition],repos,s,p0,w);
		   // if(i==18)
		   //  printf("%d,%d,%d,%f\n",i,j,l,prob_sweep_fun(n1,x1[l],q2[l],repos,s,p0,w));
	         }      
       
		//		printf("tempLogP: %f, tempMaxP,%f\n", tempLogP, tempMaxP);
		if(tempLogP>=tempMaxP)
		  {
		   
		    tempMaxP=tempLogP;
		    maxSS[i]=s;
		  }
		//		printf("%f %f %f\n",s,tempLogP,tempMaxP);
	      }
	    if(gridnumber==1)
	      {
		*maxS=maxSS[i];
		return tempMaxP;
	      }

	    logPCurve[i]=tempMaxP;
	    //   printf("%f\n",tempMaxP);
	  }
	tempMaxP=logPCurve[0];
	tempMaxPosition=0;
	for(i=0;i<gridnumber;i++)
	  {
	    //  tempMaxP=GEi(tempMaxP,logPCurve[i])*tempMaxP+(1-GEi(tempMaxP,logPCurve[i]))*logPCurve[i];
	    //  tempMaxPosition=GEi(tempMaxP,logPCurve[i])*tempMaxPosition+(1-GEi(tempMaxP,logPCurve[i]))*i;
	    if(tempMaxP<=logPCurve[i])
	      {
		tempMaxP=logPCurve[i];
		*maxS=maxSS[i];
		tempMaxPosition=i;
	      }
	    else
	      continue;
	      
	  }
	value=tempMaxP;
       
	*logPLocation=(gridstart+gridsize*tempMaxPosition)*totalGD+0.00001;
	//	printf("%f\n",*logPLocation);

		/*	for(i=0;i<nnsegsites;i++)
	  printf("%f ",pos[i]);
	  printf("\n"); */
	return value;    
  }


double logPNeutralityScan(double w,double p0,double *q2, int *x1,int *n1,int startP,int endP)
      {
	int i,j;
	double value=0.0,tempValue,tempValue2;
	extern double sig2;
	extern double q2_I;
	extern double cc_I;
	extern int n1_I;
	extern int x1_I;

	cc_I=1.0;
	for(i=startP;i<endP;i++)
	  {
	    sig2=w*q2[i]*(1-q2[i]);
	    x1_I=x1[i];
	    q2_I=q2[i];         
	    n1_I=n1[i];
	    //  tempValue=log(qromb(func1, 0.001, 0.999))-log(qromb(func10,0.001,0.999));//this part is for seg2pop case.
	    tempValue2=log(qromb(func10,0.001,0.999)/*+func10(0.0)+func10(1.0)*/);
	    tempValue=log(qromb(func1_seg2p,0.001,0.999)/*+func1(0.0)+func1(1.0)*/);
	    value+=(tempValue-tempValue2);
	    // printf("%f %f %f %f\n",func10(0.0),func10(1.0),tempValue2,func10(0.0)+func10(1.0)+tempValue2);
	    // if(fabs(tempValue-log(qromb(func3,0.001,0.999)))>1.0)
	    //  printf("mistake, %f %f\n",tempValue, log(qromb(func3,0.00001,0.99999)));
	    //   printf("neutrality:%d n1 %d x1 %d %f\n",i,n1_I, x1_I,tempValue);
	  }
	return value;

      }  

double logPSelectionScan(double w,double p0,int *x1, double *q2, double *pos,int *n1, int mutPos,int startP, int endP,double *maxS,double *sValue,int sLength)

      {
	int j,l;
	double value=0.0;
	double tempLogP,tempValue;
	double pos0;   
	double s,repos,maxs;
	double maxLogP=-1800.0;
      
	pos0=pos[mutPos]+0.00001;
        for(j=0;j<sLength/*sN*/;j++)
         {
	   //	s=sMin+j*sInc;
	   s=sValue[j];
	   tempLogP=0.0;
           for(l=startP;l<endP;l++)
	     {
               repos=max(fabs(pos[l]-pos0),0.00001);
		   // if(mutPos==5)
		   //   printf("mupos, %d, j %d, l %d,pos0 %f, repos %f\n",mutPos,j,l, pos0, repos);
	       tempValue=prob_sweep_fun(n1[l],x1[l],q2[l],repos,s,p0,w);
	       tempLogP+=tempValue;
	       //   printf("selection,%d %f %f\n",l,s,tempValue);	
	   // if(i==18)
		   // printf("%d,%d,%d,%f\n",i,j,l,prob_sweep_fun(n1,x1[l],q2[l],repos,s,p0,w));
	     }      
		//		printf("tempLogP: %f, tempMaxP,%f\n", tempLogP, tempMaxP);)
   
		//	printf("%d %d\n",i,j);
	   if(tempLogP>=maxLogP)
	     {
	       maxLogP=tempLogP;
	       maxs=s;
	     }

		// 	printf("s:%e tempLogP %e\n",s,tempLogP);
         }
	  
       value=maxLogP;
       *maxS=maxs;
       //       printf("logP: %e maxS %e \n",maxLogP,maxS);
       return value;    
  }

int getWeight(int **data2,int hapN,int *indexx,int *snpList,int snpListLength,double *wt,double corrLevel)
 {
  int i,j,k;
  int snpI, weigthI;
  double tempValue;
  int cAB,cA,cB;
  double pAB,pA,pB;
  double corr; 
  if(corrLevel<0.0001)
    {
      for(i=0;i<snpListLength;i++)
	{
	  wt[i]=1;
	}
    }
  else
    {
      for(i=0;i<snpListLength;i++)
	{
	  tempValue=0.;
	  for(j=0;j<snpListLength;j++)
	    {
	      cAB=0;
	      cA=0;
	      cB=0;
	      for(k=0;k<hapN;k++)
		{
		  cA+=data2[k][indexx[snpList[i]]];
		  cB+=data2[k][indexx[snpList[j]]];
		  cAB+=data2[k][indexx[snpList[i]]]*data2[k][indexx[snpList[j]]];
		}
	      pAB=((double)cAB)/((double)hapN);
	      pA=((double)cA)/((double)hapN);
	      pB=((double)cB)/((double)hapN);
	      corr=fabs((pAB-pA*pB)/sqrt(pA*(1.0-pA)*pB*(1.0-pB)));
	      if(corr>=corrLevel)
		tempValue+=1.0;
	      else
		tempValue+=0.;
	    }
	  wt[i]=1/tempValue;
      //  printf("%f ",tempValue);
	}
//  printf("\n");
    }  
return 0;
 }


int getWeight(int **data2,int hapN,int *indexx,int *snpList,int snpListLength,double *wt,double corrLevel,int phaseFlag)
{
  int i,j,k;
  int snpI, weigthI;
  double tempValue;
  int cAB,cA,cB;
  double pAB,pA,pB;
  double corr; 
  int *vec1,*vec2;


  if(phaseFlag==1)
    {
      for(i=0;i<snpListLength;i++)
       {
         tempValue=0.;
         for(j=0;j<snpListLength;j++)
	   {
	    cAB=0;
	    cA=0;
	    cB=0;
	    for(k=0;k<hapN;k++)
	     {
	      cA+=data2[k][indexx[snpList[i]]];
	      cB+=data2[k][indexx[snpList[j]]];
	      cAB+=data2[k][indexx[snpList[i]]]*data2[k][indexx[snpList[j]]];
	     }
	    pAB=((double)cAB)/((double)hapN);
	    pA=((double)cA)/((double)hapN);
	    pB=((double)cB)/((double)hapN);
	    corr=fabs((pAB-pA*pB)/sqrt(pA*(1.0-pA)*pB*(1.0-pB)));
	    if(corr>=corrLevel)
	      tempValue+=1.0;
	    else
	      tempValue+=0.;
	  }
        wt[i]=1/tempValue;
         //  printf("%f ",tempValue);
       }
    }//end of loop(unphaseFlag==0)
  else 
    {
      vec1=ivector(0,hapN);
      vec2=ivector(0,hapN);
      for(i=0;i<snpListLength;i++)
	{
	  tempValue=0;
	  for(j=0;j<snpListLength;j++)
	    {
	      for(k=0;k<hapN;k++)
		{
		  vec1[k]=data2[k][indexx[snpList[i]]];
		  vec2[k]=data2[k][indexx[snpList[j]]];
		}
	      phase2loci_em(vec1,vec2,hapN,&pAB,&pA,&pB);
	      corr=fabs((pAB-pA*pB)/sqrt(pA*(1.0-pA)*pB*(1.0-pB)));
	      if(corr>=corrLevel)
		tempValue+=1.0;
	      else
		tempValue+=0.0;
	    }
	  wt[i]=1/tempValue;
	}
      free_ivector(vec1,0,hapN);
      free_ivector(vec2,0,hapN);
    }//end of (unphaseFlag==1)
  return 0;
}


int getSnpWindow(double *pos,int snpN,int *snpList, int *snpListLength, int windowSnpNumber,double mutGpos, double geneticWindowSize)
{
  int i,j;
  int mutP;
  int tempLeftN=0,tempRightN=0,leftEnd=0, rightEnd=0;
  double halfGeneticWindow=geneticWindowSize/2.0;
  for(i=0;i<snpN;i++)
    {
      if(pos[i]>mutGpos)
	{
	  mutP=i-1;
	  break;
	}
    }
  for(i=mutP;i>=0;i--)
    {
      if(fabs(pos[i]-mutGpos)<=halfGeneticWindow)
	tempLeftN++;
      else
	{
	  leftEnd=i+1;
	  break;
	}
    }
  leftEnd=MAX(0,leftEnd);


  for(i=mutP+1;i<=snpN;i++)
    {
      if(fabs(pos[i]-mutGpos)<=halfGeneticWindow)
	{
	  tempRightN++;
	}
      else
	{
	  rightEnd=i-1;
	  break;
	}
    }
      rightEnd=MIN(rightEnd,snpN);

      if((tempLeftN+tempRightN)<windowSnpNumber)
	{
	  *snpListLength=(tempLeftN+tempRightN);
	  for(i=0;i<tempLeftN;i++)
	    {
	      snpList[i]=leftEnd+i;
	    }
	  for(i=tempLeftN,j=0;i<*snpListLength;i++,j++)
	    {
	      snpList[i]=j+mutP;
	    }
	  return 0;
	}
      else
	{
	  *snpListLength=(tempLeftN+tempRightN);
	  /*  if(*snpListLength>5000)
	 { printf("snpListL: %d\n",*snpListLength);
	  exit(0);
	  }*/
	  for(i=0;i<tempLeftN;i++)
	    {
	      snpList[i]=leftEnd+i;
	    }
	  for(i=tempLeftN,j=0;i<*snpListLength;i++,j++)
	    {
	      snpList[i]=j+mutP;
	    }
	  return 1;
	}

}
/*
int getSnpWindowPPGW(double *pos,int snpN,int *snpList, int *snpListLength, int windowSnpNumber,double mutGpos, double geneticWindowSize)
{
  int i,j;
  int mutP;
  int tempLeftN=0,tempRightN=0,leftEnd=0, rightEnd=0;
  double halfGeneticWindow=geneticWindowSize/2.0;
  for(i=0;i<snpN;i++)
    {
      if(pos[i]>mutGpos)
	{
	  mutP=i-1;
	  break;
	}
    }
  for(i=mutP;i>=0;i--)
    {
      if(fabs(pos[i]-mutGpos)<=halfGeneticWindow)
	tempLeftN++;
      else
	{
	  leftEnd=i+1;
	  break;
	}
    }
  leftEnd=MAX(0,leftEnd);


  for(i=mutP+1;i<=snpN;i++)
    {
      if(fabs(pos[i]-mutGpos)<=halfGeneticWindow)
	{
	  tempRightN++;
	}
      else
	{
	  rightEnd=i-1;
	  break;
	}
    }
      rightEnd=MIN(rightEnd,snpN);

      if((tempLeftN+tempRightN)<windowSnpNumber)
	{
	  *snpListLength=(tempLeftN+tempRightN);
	  for(i=0;i<tempLeftN;i++)
	    {
	      snpList[i]=leftEnd+i;
	    }
	  for(i=tempLeftN,j=0;i<*snpListLength;i++,j++)
	    {
	      snpList[i]=j+mutP;
	    }
	  return 0;
	}
      else
	{
	  *snpListLength=(tempLeftN+tempRightN);
	  //  if(*snpListLength>5000)
	  // { printf("snpListL: %d\n",*snpListLength);
	  // exit(0);
	  // }
	  for(i=0;i<tempLeftN;i++)
	    {
	      snpList[i]=leftEnd+i;
	    }
	  for(i=tempLeftN,j=0;i<*snpListLength;i++,j++)
	    {
	      snpList[i]=j+mutP;
	    }
	  return 1;
	}

}

*/

int getSnpWindow_ms(int snpN, int *snpList, int *snpListLength, int windowSnpNumber)
{
  int i,j;
  for(i=0;i<snpN;i++)
    {
      snpList[i]=i;
    }
  *snpListLength=snpN;

  if(snpN<windowSnpNumber)
    return(0);
  else
    return(1);
}

int getSubsetSnpWindow(int windowSnpNumber, double geneticWindowSize,int *snpList, int *snpListLength, int *idum)
{
  int i,j;
  int tempID;
  int snpNN=*snpListLength;
  for(i=snpNN;i>windowSnpNumber;i--)
    {
      tempID=(int)(i*ran1(idum));
      snpList[tempID]=snpList[snpNN];
      snpNN--;
    }
  *snpListLength=snpNN;
return 0;
}

int getSubsetSnpWindow_ms(int windowSnpNumber,int *snpList,int *snpListLength,int *idum)
{
  int i,j;
  int tempID;
  int snpNN=*snpListLength;
  for(i=snpNN;i>windowSnpNumber;i--)
    {
      printf("begin here\n");
      tempID=(int)(i*ran1(idum));
      snpList[tempID]=snpList[snpNN];
      printf("\ntempID: %d %d \n",tempID,snpList[tempID]);
      snpNN--;
    }
  *snpListLength=snpNN;
return 0;
}


/* adjust for fixed windowsize and fixed snp number */
double logPSelectionScan3(double w,double p0,int *x1, double *q2, double *pos,int *n1, double mutGpos,double *maxS,double *sValue,int sLength, int *snpList, int SnpListLength,double* wt)

      {
	int j,l;
	double value=0.0;
	double tempLogP,tempValue;
	double pos0;   
	int Oposition;
	double s,repos,maxs;
	double maxLogP=-1800.0;
	/*	for(j=0;j<SnpListLength;j++)
	  {
	    printf("mutGpos:%f snplistLength:%d j:%d snpList[j]:%d pos:%f n1:%d x1:%d q2:%f\n",mutGpos, SnpListLength,j,snpList[j],pos[snpList[j]],n1[snpList[j]],x1[snpList[j]],q2[snpList[j]]);

	    }
	*/      
	pos0=mutGpos;
        for(j=0;j<sLength/*sN*/;j++)
         {
	   //	s=sMin+j*sInc;
	   s=sValue[j];
	   tempLogP=0.0;
           for(l=0;l<SnpListLength;l++)
	     {
	       Oposition=snpList[l];
               repos=max(fabs(pos[Oposition]-pos0),0.000001);
	       //  printf("\n n1:%d x1:%d q2:%f repos:%f s:%f p0:%f w:%f\n",n1[Oposition],x1[Oposition],q2[Oposition],repos,s,p0,w);
	       tempValue=prob_sweep_fun(n1[Oposition],x1[Oposition],q2[Oposition],repos,s,p0,w);
	       //    printf("s %d, snp %d; st %f tempValue %f \n",j,l,wt[l],tempValue);
	       tempLogP+=wt[l]*tempValue;
	       //   printf("selection,%d %f %f\n",l,s,tempValue);	
	     }      
		//		printf("tempLogP: %f, tempMaxP,%f\n", tempLogP, tempMaxP);)
		//	printf("%d %d\n",i,j);
	   if(tempLogP>=maxLogP)
	     {
	       maxLogP=tempLogP;
	       maxs=s;
	     }
		// 	printf("s:%e tempLogP %e\n",s,tempLogP);
         }
	  
       value=maxLogP;
       *maxS=maxs;
       //       printf("logP: %e maxS %e \n",maxLogP,maxS);
       return value;    
  }


/* adjust for fixed windowsize and fixed snp number */
double logPSelectionScan2(double w,double p0,int *x1, double *q2, double *pos,int *n1, double mutGpos,double *maxS,double *sValue,int sLength, int *snpList, int SnpListLength)

      {
	int j,l;
	double value=0.0;
	double tempLogP,tempValue;
	double pos0;   
	int Oposition;
	double s,repos,maxs;
	double maxLogP=-1800.0;
	/*	for(j=0;j<SnpListLength;j++)
	  {
	    printf("%d %d %d %f %f\n",SnpListLength,j,snpList[j],pos[snpList[j]],(pos[mutPos]-pos[snpList[j]]));

	    }*/      
	pos0=mutGpos;
        for(j=0;j<sLength/*sN*/;j++)
         {
	   //	s=sMin+j*sInc;
	   s=sValue[j];
	   tempLogP=0.0;
           for(l=0;l<SnpListLength;l++)
	     {
	       Oposition=snpList[l];
               repos=max(fabs(pos[Oposition]-pos0),0.000001);
	       tempValue=prob_sweep_fun(n1[Oposition],x1[Oposition],q2[Oposition],repos,s,p0,w);
	       tempLogP+=tempValue;
	       //   printf("selection,%d %f %f\n",l,s,tempValue);	
	     }      
		//		printf("tempLogP: %f, tempMaxP,%f\n", tempLogP, tempMaxP);)
		//	printf("%d %d\n",i,j);
	   if(tempLogP>=maxLogP)
	     {
	       maxLogP=tempLogP;
	       maxs=s;
	     }
		// 	printf("s:%e tempLogP %e\n",s,tempLogP);
         }
	  
       value=maxLogP;
       *maxS=maxs;
       //       printf("logP: %e maxS %e \n",maxLogP,maxS);
       return value;    
  }


double logPNeutralityScan3(double w,double p0,double *q2, int *x1,int *n1,int *snpList, int SnpListLength,double *wt)
      {
	int i,j;
	double value=0.0,tempValue,tempValue2;
	int Oposition;
	extern double sig2;
	extern double q2_I;
	extern double cc_I;
	extern int n1_I;
	extern int x1_I;

	cc_I=1.0;
	for(i=0;i<SnpListLength;i++)
	  {
	    Oposition=snpList[i];
	    sig2=w*q2[Oposition]*(1-q2[Oposition]);
	    x1_I=x1[Oposition];
	    q2_I=q2[Oposition];         
	    n1_I=n1[Oposition];
	    //  tempValue=log(qromb(func1, 0.001, 0.999))-log(qromb(func10,0.001,0.999));//this part is for seg2pop case.
	    tempValue2=log(qromb(func10,0.001,0.999)/*+func10(0.0)+func10(1.0)*/);
	    tempValue=log(qromb(func1_seg2p,0.001,0.999)/*+func1(0.0)+func1(1.0)*/);
	    value+=(wt[i])*(tempValue-tempValue2);
	    // printf("%f %f %f %f\n",func10(0.0),func10(1.0),tempValue2,func10(0.0)+func10(1.0)+tempValue2);
	    // if(fabs(tempValue-log(qromb(func3,0.001,0.999)))>1.0)
	    //  printf("mistake, %f %f\n",tempValue, log(qromb(func3,0.00001,0.99999)));
	    //   printf("neutrality:%d n1 %d x1 %d %f\n",i,n1_I, x1_I,tempValue);
	  }
	return value;

      }  

double logPNeutralityScan2(double w,double p0,double *q2, int *x1,int *n1,int *snpList, int SnpListLength)
      {
	int i,j;
	double value=0.0,tempValue,tempValue2;
	int Oposition;
	extern double sig2;
	extern double q2_I;
	extern double cc_I;
	extern int n1_I;
	extern int x1_I;

	cc_I=1.0;
	for(i=0;i<SnpListLength;i++)
	  {
	    Oposition=snpList[i];
	    sig2=w*q2[Oposition]*(1-q2[Oposition]);
	    x1_I=x1[Oposition];
	    q2_I=q2[Oposition];         
	    n1_I=n1[Oposition];
	    //  tempValue=log(qromb(func1, 0.001, 0.999))-log(qromb(func10,0.001,0.999));//this part is for seg2pop case.
	    tempValue2=log(qromb(func10,0.001,0.999)/*+func10(0.0)+func10(1.0)*/);
	    tempValue=log(qromb(func1_seg2p,0.001,0.999)/*+func1(0.0)+func1(1.0)*/);
 value+=(tempValue-tempValue2);
	    // printf("%f %f %f %f\n",func10(0.0),func10(1.0),tempValue2,func10(0.0)+func10(1.0)+tempValue2);
	    // if(fabs(tempValue-log(qromb(func3,0.001,0.999)))>1.0)
	    //  printf("mistake, %f %f\n",tempValue, log(qromb(func3,0.00001,0.99999)));
	    //   printf("neutrality:%d n1 %d x1 %d %f\n",i,n1_I, x1_I,tempValue);
	  }
	return value;
      }

  
/*
struct snpNode
{  
  int position;
  snpNode *next;
};

typedef struct snpNode ELEMENT;
typedef ELEMENT *LINK;

LINK position2list(int position1, int position2, int mutPosition)
{
  LINK head;
  int i,j;
  head=malloc(sizeof(ELEMENT));
  head->position=position1;

  for(i=position1+1;i<=position2;i++)
    {
      if(i!=mutPosition&&i<position2)
	{
	  head->next=malloc(sizeof(ELEMENT));
	  head=head->next;
	  head->position=i;
	}
      else if(i==position2)
	{
	  head->next=NULL;
	  head->position=i;
	}
      else
	continue;
    }
  return head;
}

void delMinEntry(LINK head,int count,double *pos)
{
  int i,j;
  double minDist=0., upScore=0., downScore=0.;
  double tempScore;
  LINK thead=head;
  int flag;

  for(i=0;i<count;i++)
    {
      downScore= fabs(pos[thead->position]-pos[thead->next->position]);
      tempScore= upScore +downScore;
      upScore=downScore;
      if(tempScore<=minDist)
	{flag=i;
	minDist=tempScore;
	}
      thead=head->next;
    }
  j=0;
  thead=head;
  for(i=0;i<count;i++)
    {
      if(i==flag-1)
	{
	  thead->next=thead->next->next;
	}
      else
	thead=head->next;
    }

 

}

int pickSNP4Window(int mutPos,int snpN, double *pos,double geneticWindowSize, int windowSnpNumber, int *snpIndex)
  {
    int tempLeftN=0, tempRightN=0,leftEnd,rightEnd;
    int i,j,k,l;
    double halfGeneticWindow=geneticWindowSize/2.0;
    //  int halfSNPN=windowSnpNumber/2;
    int totalN;
    int *tempIndex;
    LINK head;

    for(i=mutPos-1;i>=0;i--)
      {
	if(fabs(pos[i]-pos[mutPos])<=halfGeneticWindow)
	  {
	    tempLeftN++;
	  }
	else
	  {
	    leftEnd=i+1;
	    break;
	  }
      }

    for(i=mutPos+1;i<=snpN;i++)
      {
	if(fabs(pos[i]-pos[mutPos])<=halfGeneticWindow)
	  {tempEightN++;
	  }
	else
	  {
	    rightEnd=i-1;
	    break;
	 
      }
    totalN=tempLeftN+tempRight fabs(pos[thead->position]-pos[thead->next->position]N;

    if(totalN<windowSnpNumber)
      return 0;
    else
      {

	head=position2list(leftEnd,rightEnd);

	for(i=0;i<=leftN;i++)
	  tempIndex[i]=leftN+i;
	for(i=0;i<=rightN;i++)
	  tempIndex[i+tempLftN]=musPos+i;
	for(i=totalN;i> windowSnpNumber; i--)
	  {          

      }
    return 1;
  }

*/


 double fst(double *q2,double *q1,int n2, int n1,int nnsegsites)
    {

      int i,n,tempCount=0;
      double fst,p1,p2;
      double nc,MSP,MSG,p_hat;
      double *fst_v;
      double mean_d(double *v,int l);
      n=n1+n2;
      double s=2.0;
      nc= (1/(s-1.0))*((n1+n2)-(double)(pow((double)n1,2)+pow((double)n2,2))/(n1+n2));
      if( ! ( fst_v = ( double*) malloc( (unsigned)(nnsegsites*sizeof(double )) ) ) )
      perror("alloc error in fst_v") ;
      for(i=0;i<nnsegsites;i++)
      {
	p1=q1[i];
	p2=q2[i];
	p_hat=((double)n1/n)*p1+((double)n2/n)*p2;
	MSP=(1/(s-1.0))*(n1*pow((p1-p_hat),2)+n2*pow((p2-p_hat),2));
	MSG=(1.0/((n1-1)+(n2-1)))*(n1*p1*(1.0-p1)+n2*p2*(1.0-p2));
	if(fabs(MSP+(nc-1.0)*MSG)>0.000001){
	  fst_v[i]=(MSP-MSG)/(MSP+(nc-1.0)*MSG);}
	else
	  fst_v[i]=0.0;
      }
      fst=mean_d(fst_v,nnsegsites);
      return(fst);
    }

double mean_d(double *v,int l)
{

  int i,count=0;
  double tempValue=0.0;
  for(i=0;i<l;i++)
    {
      if(v[i]>0.0)
	{
	  tempValue+=v[i];
	  count+=1;
	}
    }
  if(count>0)
  tempValue=tempValue/count;
  else
    tempValue=0.0;
  return tempValue;
}

 
//EM algorithm for inferring correlation of two SNP loci.

int phase2loci_em(int *snpv1,int *snpv2,int hapN,double* pAB,double* pA, double* pB)
{
  int i,j;
  int cAB,cA,cB,nA,nB,n1,n2, xhet=0;
  double **f;
  int **data_count;
  double logL,tolerance=0.00001;  
  double dx=100*tolerance;
  double theta=-999999.0;
  double theta_old;
  f=dmatrix(0,1,0,1);
  data_count=imatrix(0,1,0,1);
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      data_count[i][j]=0;

  nA=nB=0;
  n1=n2=0;
  for(i=0;i<2;i++)
    {
      for(j=0;j<2;j++)
	{
	  f[i][j]=0.25;
	}
    }

  for(i=0;i<hapN;i++)
    {
      nA+=1*EQUAL(snpv2[i],1);
      n1+=1*EQUAL(snpv2[i],1)+1*EQUAL(snpv2[i],0);
      nB+=1*EQUAL(snpv1[i],1);
      n2+=1*EQUAL(snpv1[i],1)+1*EQUAL(snpv1[i],0);
      //      printf("%d %d %d %d %d %d\n",snpv1[i],snpv2[i],nA,n1,nB,n2);
    }
  //for cases only one locus segregating.
  if(nA==n1)
    {
      *pA=1;
      *pB=((double)nB)/((double)hapN);
      *pAB=*pB;
      return 1;
    }
  else if(nA==0)
    {
      *pA=0;
      *pB=((double)nB)/((double)hapN);
      *pAB=0;
      return 1;
    }
  else if(nB==n2)
    {
      *pB=1.0;
      *pA=((double)nA)/((double)hapN);
      *pAB=*pA;
      return 1;
    }
  else if(nB==0)
    {
      *pB=0.0;
      *pA=((double)nA)/((double)hapN);
      *pAB=0.0;
      return 1;
    }
  else 
    {
      ;
    }


  for(i=0;i<hapN/2;i++)
    {
      if((EQUAL((snpv2[i*2+0]+snpv2[i*2+1]),1))&&EQUAL((snpv1[i*2+0]+snpv1[i*2+1]),1))
	{
	  xhet+=1;
	}
      else
	{
	  data_count[0][0]+=(EQUAL(snpv2[i*2+0],1)*EQUAL(snpv1[i*2+0],1)+EQUAL(snpv2[i*2+1],1)*EQUAL(snpv1[i*2+1],1));
	  data_count[0][1]+=(EQUAL(snpv2[i*2+0],1)*EQUAL(snpv1[i*2+0],0)+EQUAL(snpv2[i*2+1],1)*EQUAL(snpv1[i*2+1],0));
	  data_count[1][0]+=(EQUAL(snpv2[i*2+0],0)*EQUAL(snpv1[i*2+0],1)+EQUAL(snpv2[i*2+1],0)*EQUAL(snpv1[i*2+1],1));
	  data_count[1][1]+=(EQUAL(snpv2[i*2+0],0)*EQUAL(snpv1[i*2+0],0)+EQUAL(snpv2[i*2+1],0)*EQUAL(snpv1[i*2+1],0));
	}
    }
 
  // printf("datacout: %d %d %d %d %d\n",data_count[0][0],data_count[0][1],data_count[1][0],data_count[1][1],xhet);

  int kk=0;
  while((dx>tolerance)&&(kk<1000))
    {
      theta_old=theta;

      for(i=0;i<2;i++)
	{
	  for(j=0;j<2;j++)
	    {
	      f[i][j]=(data_count[i][j]+(EQUAL(i,j))*(1-theta)*((double)xhet) + (1-EQUAL(i,j))*(theta*((double)xhet)))/((double)hapN);
	    }
	}
      theta=(f[0][1]*f[1][0])/(f[0][0]*f[1][1]+f[0][1]*f[1][0]);
      logL=logL_em(f,data_count);
      kk++;
      dx=fabs(theta-theta_old);
      //      printf("p00 %f p01 %f p10 %f p11:%f\n",f[0][0],f[0][1],f[1][0],f[1][1]);
    }
  *pAB=f[0][0];
  *pA=nA/((double)hapN);
  *pB=nB/((double)hapN);
  free_dmatrix(f,0,1,0,1);
  free_imatrix(data_count,0,1,0,1);
  return 1;
}

//esitmate the complete log-likelihood of the data for 2-loci-em.

double logL_em(double **f, int **count)
{
  double logL=0.0;
  int i,j;
  for(i=0;i<=1;i++)
    {
      for(j=0;j<=1;j++)
	{
	  //	  printf("p:%f count %d ;",f[i][j],count[i][j]);
	  logL+=count[i][j]*log(f[i][j]);
	}
    }
  //  printf("\nlogL %f\n",logL);
  return logL;
}

void get_haplotypeFreq(double *pB2,double *pAB2,double *pAb2,double *pB1,double *pAB1,double *pAb1,int *snpv1_1,int *snpv2_1,int hapN1,int *snpv1_2,int *snpv2_2,int hapN2,int *x1_n)
{
  double PB1,PB2,PAB1,PAb1,PaB1,Pab1,PAB2,PAb2,PaB2,Pab2,PA1,PA2;
  double Ptemp;
  int i,j;

  phase2loci_em(snpv1_2,snpv2_2,hapN2,&PAB2,&PA2,&PB2);
  phase2loci_em(snpv1_1,snpv2_1,hapN1,&PAB1,&PA1,&PB1);
  PAb1=PA1-PAB1;
  PAb2=PA2-PAB2;
  PaB1=PB1-PAB1;
  PaB2=PB2-PAB2;
  Pab1=1-PAB1-PAb1-PaB1;
  Pab2=1-PAB2-PAb2-PaB2;
  if(PB1>=PB2)
    {
      if(PAB2<PaB2)
	{
	  Ptemp=PAB2;
	  PAB2=PaB2;
	  PaB2=Ptemp;
	  Ptemp=PAb2;
	  PAb2=Pab2;
	  Pab2=Ptemp;
	  Ptemp=PAB1;
	  PAB1=PaB1;
	  PaB1=Ptemp;
	  Ptemp=PAb1;
	  PAb1=Pab1;
	  Pab1=Ptemp;
	  *x1_n=(int)((PAB1+PAb1)*hapN1);
	}
      else
	{
	  *x1_n=(int)((PAB1+PAb1)*hapN1); 
	}
    }
  else
    {
      if(PAb2<Pab2)
	{
	  PB2=1.0-PB2;
	  PB1=1.0-PB1;
	  PAB2=Pab2;
	  PAb2=PaB2;
	  PAB1=Pab1;
	  PAb1=PaB1;
	  *x1_n=(int)((PAB1+PAb1)*hapN1);
	}
      else
	{
	  PB2=1.0-PB2;
	  PB1=1.0-PB1;
	  Ptemp=PAB2;
	  PAB2=PAb2;
	  PAb2=Ptemp;
	  Ptemp=PAB1;
	  PAB1=PAb1;
	  PAb1=Ptemp;
	  *x1_n=(int)((PAB1+PAb1)*hapN1);
	}
    }
 
  *pB1=PB1;
  *pB2=PB2;
  *pAB1=PAB1;
  *pAB2=PAB2;
  *pAb1=PAb1;
  *pAb2=PAb2;
}  


/******************************************************************************************************************
Subroutines for read in different format files
*******************************************************************************************************************/

int readinGenoEigenFormat(FILE *eigenGenoInput,char **genoEigen,int eigenIndN)
{
  int i,j,k;
  char line[5000];
  int count=0;
  while(fgets(line,5000,eigenGenoInput)!=NULL)
    {
      sscanf(line,"%s\n",genoEigen[count]);
      count++;
    }
  return count;
}

int readinSNPEigenFormat(FILE *eigenSNPInput,char **snpNameEigen,int *snpChrEigen,double *snpGPEigen,int *snpPosEigen, char **snpAlleleEigen)
{
  int eigenSNPN=0;
  int count=0;
  char line[5001];
  char tempS[100];
  while(fgets(line,5000,eigenSNPInput)!=NULL)
    {
      sscanf(line,"%s %d %s %d %c %c\n", snpNameEigen[count],&snpChrEigen[count],tempS,&snpPosEigen[count],&snpAlleleEigen[count][0],&snpAlleleEigen[count][1]);
      snpGPEigen[count]=atof(tempS);
      count++;
    }

  eigenSNPN=count;
  return eigenSNPN;
}

int readinIndEigenFormat(FILE *eigenIndInput,char **indIDEigen,char *indSexEigen,char **indPopEigen)
{
  int eigenIndN=0;
  int i,j;
  char line[5001];
  int count=0;
  while(fgets(line,5000,eigenIndInput)!=NULL)
    {
      sscanf(line,"%s %c %s\n",indIDEigen[count],&indSexEigen[count],indPopEigen[count]);
      count++;
    }
  eigenIndN=count;
  return eigenIndN;
}


int readinGenoXpclrFormat(FILE *xpclrGenoInput,char **genoXpclr, int xpclrIndN)
{
  int i,j;
  char line[5001];
  int count=0;
  while(fgets(line,5000,xpclrGenoInput)!=NULL)
    {
      sscanf(line,"%s\n",genoXpclr[count]);
      count++;
    }
  return count; 
}
/*

int readinSNPXpclrFormat(FILE *xpclrSNPInput, char *snpNameXpclr,int *snpChrXpclr,double *snpGPXpclr,int *snpPPXpclr,char**snpAlleleXpclr)
{
  int xpclrSNPN=0;
  int i,j;
  char line[5001];
  char tempString[100];
  while(fgets(line,5000,xpclrSNPInput)!=NULL)
    {
      sscanf(line,"%s %d %s %d %c %c\n",snpNameXpclr[xpclrSNPN],&snpChrXpclr[xpclrSNPN],tempString,&snpPPXpclr[xpclrSNPN],&snpAlleleXpclr[xpclrSNPN][0],&snpAlleleXpclr[xpclrSNPN][1]);
      snpGPXpclr[xpclrSNPN]=atof(tempString);
      xpclrSNPN++;
    }

}
*/

				     


/*******************************************
The following subroutines are for XPCLR4S.   

*******************************************/
 		     

double logPSelection4ss_power(int **data1,int **data2,int hapN1,int hapN2,double w,double *pos,int *n1, int mutpos,double *maxS,double *sValue,int sLength,int *snpList,int SnpListLength,double *wt,int *snpv1_1,int *snpv2_1,int *snpv1_2,int *snpv2_2,int *indexHC,double *p2,double *p1,int *x1)
     {
       int j,l,k,i,x1_1;
       double value=0.0;
       double tempLogP,tempValue;
       double pos0;
       int Oposition;
       double s,repos,maxs;
       double maxLogP=-1800.0;
       double pB1,pB2,pAB1,pAB2,pAb1,pAb2;
       pos0=pos[mutpos];
       for(i=0;i<hapN1;i++)
	 {
	   snpv1_1[i]=data1[i][indexHC[mutpos]];
	 }
       for(i=0;i<hapN2;i++)
	 {
	   snpv1_2[i]=data2[i][indexHC[mutpos]];
	 }

       for(j=0;j<sLength;j++)
	 {
	   s=sValue[j];
	   tempLogP=0.0;
	   for(l=0;l<SnpListLength;l++)
	     {
	       if(snpList[l]==mutpos)
		 continue;
	       else
		 {
		   Oposition=snpList[l];
		   repos=fabs(pos[Oposition]-pos0);
		   for(i=0;i<hapN1;i++)
		     {
		       snpv2_1[i]=data1[i][indexHC[Oposition]];
		     }
		   for(i=0;i<hapN2;i++)
		     {
		       snpv2_2[i]=data2[i][indexHC[Oposition]];
		     }
		 
		   get_haplotypeFreq(&pB2,&pAB2,&pAb2,&pB1,&pAB1,&pAb1,snpv1_1,snpv2_1,hapN1,snpv1_2,snpv2_2,hapN2,&x1_1);
		   //	   printf("snp1,%d snp2, %d pB2,%f pB1,%f  pA2,%f pA1,%f p2A,%f p1A,%f p2B %f, p1B %f\n",mutpos,Oposition,pB2,pB1,(pAb2+pAB2),(pAb1+pAB1),p2[Oposition],p1[Oposition],p2[mutpos],p1[mutpos]);
		   //		   x1_1=(int)((pAB1+pAb1)*n1[Oposition]);
		   //		   x1_1=x1[Oposition];
		   //		   printf("mutpos: %d pB2 %f pB1 %f; p2: %f p1: %f\n",mutpos,p2[mutpos],p1[mutpos],p2[snpList[l]],p1[snpList[l]]);
		   //		   printf("freq:pB2,%f pB1,%f pAB2,%f pAb2,%f pAB1,%f pAb1,%f\n",pB2,pB1,pAB2,pAb2,pAB1,pAb1);
		   //		   printf("repos %f %f %f\n",repos,pos0,pos[Oposition]);
		   tempValue= prob_softSweep_fun(n1[Oposition],x1_1,pB2,pB1,pAB2,pAb2,repos,s,w);
		   tempLogP+=wt[l]*tempValue;
		 }
	     }
	   if(tempLogP>maxLogP)
	     {
	       maxLogP=tempLogP;
	       maxs=s;
	     }
	 }
       value=maxLogP;
       *maxS=maxs;
       return value;
     }

double logPNeutrality4ss_power(double w,double *p2, int *x1,int *n1, int *snpList, int SnpListLength, double *wt, int mutPos)
     {
       int i,j;
       double value=0.0,tempValue,tempValue2;
       extern double w_I;
       extern double sig2;
       extern double q2_I;
       extern int n1_I;
       extern int x1_I;
       w_I=w;
       for(i=0;i<SnpListLength;i++)
	 {
	   if(snpList[i]==mutPos)
	     continue;
	   else
	     {
	       sig2=w*p2[snpList[i]]*(1-p2[snpList[i]]);
	       x1_I=x1[snpList[i]];
	       q2_I=p2[snpList[i]];
	       n1_I=n1[snpList[i]];
	       tempValue2=log(qromb(func10,0.001,0.999));
	       tempValue=log(qromb(func1,0.001,0.999));
	       //	       printf("neutra:%d %d %d %f %f %f %f\n",i,x1_I,n1_I,tempValue, tempValue2,wt[snpList[i]],(tempValue-tempValue2));
	       value+=(wt[i])*(tempValue-tempValue2);
	     }
	 }
       return value;
     }

/*estimate the duration based on the starting and ending allele freq*/
double get_t(double s,double p0,double pt)
     {
       double tValue;
       pt=MIN(pt,0.999);
       if(pt<p0)
	 perror("Mistake in get_t 1\n");
       else if(pt<=0)
	 perror("Mistake in get_t 2\n");
       else if(p0<=0)
	 {
	 perror("Mistake in g et_t 3\n");
	 printf("%f %f %f\n",p0,pt,s);
	 }
       else
	 {
	   if(s<0.0)
	     perror("Mistake in get_t function\n");
	   else if(s<=0.0)
	     {
	       tValue=0.0;
	       return tValue;
	     }
	   else
	     tValue=-1/s*log((p0*(1-pt)/(pt*(1-p0))));
	 }
       return tValue;
     }


double partA_pAB(double t)
    {
      double value,value1,value2;
      extern double r_I;
      extern double s_I;
      extern double pB0_I;

      value1=(1-pB0_I)*exp(-t*(r_I+s_I));
      value2=(pB0_I+(1-pB0_I)*exp(-t*s_I));
      if(value2<=0.0)
	perror("mistake in partA_pAB\n");
      else
	value=value1/value2;
      return value;
    }

double partA_pAb(double t)
    {
      double value,value1,value2;
      extern double r_I;
      extern double s_I;
      extern double pB0_I;

      value1=(pB0_I)*exp(-t*r_I);
      value2=(pB0_I+(1-pB0_I)*exp(-t*s_I));
      if(value2<=0.0)
	perror("mistake in partA_pAb\n");
      else
	value=value1/value2;
      return value;
    }

double func4ss_pA(double p)
     {
       double value;
       extern double muu_I;
       extern double varr_I;
       extern double pAB0_I;
       extern double pAb0_I;
       extern double pB0_I;
       extern double pBt_I;
       extern double s_I;
       extern double r_I;
       extern double tau_I;
       extern double w_I;
       extern int n1_I;
       extern int x1_I;
       if(p>0.0 && p<1.0)
	 {
	   value=1/sqrt(2.0*PI*varr_I)*exp(-pow((p-muu_I),2.0)/(2.0*varr_I))*binopdf(n1_I,x1_I,p);
	 }
       else if(p<=0.0)
	 {
	   value=gaussiancdf(0.0,muu_I,sqrt(varr_I))*binopdf(n1_I,x1_I,p);
	 }
       else if(p>=1.0)
	 {
	   value=(1.0-gaussiancdf(1.0,muu_I,sqrt(varr_I)))*binopdf(n1_I,x1_I,p);
	 }
       else
	 {
	   perror("miskate in soft sweep function");
	 }
       return value;
     }

double func4ss_pA_all(double p)
     {
       double value;
       extern double muu_I;
       extern double varr_I;
       extern double pAB0_I;
       extern double pAb0_I;
       extern double pB0_I;
       extern double pBt_I;
       extern double s_I;
       extern double r_I;
       extern double tau_I;
       extern double w_I;
       extern int n1_I;
       extern int x1_I;
       if(p>0.0 && p<1.0)
	 {
	   value=1/sqrt(2.0*PI*varr_I)*exp(-pow((p-muu_I),2.0)/(2.0*varr_I));
	 }
       else if(p<=0.0)
	 {
	   value=gaussiancdf(0.0,muu_I,sqrt(varr_I));
	 }
       else if(p>=1.0)
	 {
	   value=(1.0-gaussiancdf(1.0,muu_I,sqrt(varr_I)));
	 }
       else
	 {
	   perror("miskate in soft sweep function");
	 }

       return value;
     }

double prob_softSweep_fun(int n1_1,int x1_1,double pB0,double pBt,double pAB0,double pAb0,double c,double s,double w)
     {
       double value,value2;
       double pABt,pAbt,pA0;
       extern int n1_I;
       extern int x1_I;
       extern double pB0_I;
       extern double pBt_I;
       extern double s_I;
       extern double r_I;
       extern double tau_I;
       extern double w_I;
       extern double pAB0_I;
       extern double pAb0_I;
       extern double muu_I;
       extern double varr_I;

       n1_I=n1_1;
       x1_I=x1_1;
       pB0_I=pB0;
       pBt_I=pBt;
       pAB0_I=pAB0/pB0_I;
       pAb0_I=pAb0/(1.0-pB0_I);
       pA0=pAb0+pAB0;
       s_I=s;
       r_I=c;
       w_I=w;
       //   printf("s: %f pB0: %f pBt: %f\n",s_I,pB0_I,pBt_I);
       tau_I=get_t(s_I,pB0_I,pBt_I);
       // printf("t is: %f\n",tau_I);
       if(tau_I<=0.0)
	 {
	   pABt=pAB0_I;
	   pAbt=pAb0_I;
	 }
       else
	 {
	   pABt=pAB0_I-r_I*(pAB0_I-pAb0_I)*qromb(partA_pAB,0.0,tau_I);
	   pAbt=pAb0_I+r_I*(pAB0_I-pAb0_I)*qromb(partA_pAb,0.0,tau_I);
	 }
       if(tau_I<=0.0)
	 {
	   muu_I=pA0;
	 }
       else
	 {
	   muu_I=MIN(pABt*pBt_I+pAbt*(1-pBt_I),1.0);
	 }
       //   printf("muu %f %f %f %f\n",pBt,pABt,pAbt,muu_I);
       varr_I=w_I*pA0*(1-pA0);
       value=qromb(func4ss_pA,0.001,0.999);
       //   printf("func:%f %f %f %f %f %f %d %d\n ",muu_I,varr_I,s_I,r_I,tau_I,w_I,n1_I,x1_I);
       value2=qromb(func4ss_pA_all,0.001,0.999);
       //       printf("value:%d %d %f value2: %f\n",x1_I,n1_I,value,value2);
       value=log(value)-log(value2);
       
 //    if(value>=0)
	 return value;
	 //    else
	 //	 return 0;
     } 

