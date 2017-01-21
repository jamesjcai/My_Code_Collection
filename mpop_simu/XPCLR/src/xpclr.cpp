/***************************************************************************
SOFTWARE COPYRIGHT NOTICE AGREEMENT
This software and its documentation are copyright(2009) by Harvard Medical
School. All rights are reserved.

This software is supplied without any warranty or guaranteed support
what so ever. HMS cannot be responsible for its use, misuse, or functionality.

***************************************************************************/

/**************************************************************************
  Cross-Population Composite Likelihood Ratio test

      Hua Chen  hchen@genetics.med.harvard.edu
                 2009-1-30
 
***************************************************************************/
#include "xpclr.h"
#define MAXSNPLENGTH 400000
#define MAXSNPLENGTH2 800000
#define MAXHAP 3000
#define MAXSNPINWINDOW 10000


int determineCentromere(double *position, int length)
{
  double maxInterval=0;
  int i,j, Iposition=0;

  for(i=1;i<length;i++)
    {
      if((position[i]-position[i-1])>=maxInterval)
	{
	  maxInterval=(position[i]-position[i-1]);
	  Iposition=i;
	}
    }
  return Iposition;
}


int **data1, **data2;
char dataString[MAXSNPLENGTH2];
int hapN1,hapN2;
double *position;
double *pposition;
int *indexHC;
double *rawPosition;
double *rawPposition;
double *gridPosition;
double *phGridPosition;
double gridSize;
int gridNumber;
double *q1,*q2;
double *rawq1,*rawq2;
int *x1;
int centroPosition;
double q2_I,sig2;
int x1_I,n1_I,n2_I;
double w,fstValue;
double cc_I;
double re; //for transforming from unit of cM to 100cM (recombination rate).
int N=10000;
double sMin=0.00,sMax=0.1,sInc=0.002;
double sValue[16]={0.0,0.00001,0.00005,0.0001,0.0002,0.0004,0.0006,0.0008,0.001,0.003,0.005,0.01,0.05,0.08,0.1,0.15};
int sLength=16;
int sN=(int)((sMax-sMin)/sInc+1); 
double logPRatio,logPS,logPN, maxS;
int logPLocation;
double logPCurve[1000];
int snpN;
int startP,endP;
int tempInt1, tempInt2;
char outFile3[150];
int sparseIndicator=0;
int phaseIndicator=0;
double gWindowSize;
int windowSnpNumber;
int cchrn;
double ccorr=0.95;
double setStartPos;
int tempI;

//extern variables declared for xpclr4s
double w_I,r_I,s_I,tau_I,pB0_I,pBt_I,pAB0_I,pABt_I,pAb0_I,pAbt_I,varr_I,muu_I; 


int main(int argc, char* argv[])
  {
    int inputFileStyle=0;
    char popLabel1[20],popLabel2[20];
    if(argc==10 && ((strcmp("-h",argv[1])==0)))
      {
	inputFileStyle=1;
	strcpy(popLabel1,argv[4]);
	strcpy(popLabel2,argv[5]);
	gridSize=atof(argv[9]);
	if(strcmp("-w1",argv[6])==0)
	  {
	    sparseIndicator=1;
	    gWindowSize=atof(argv[7]);
	    windowSnpNumber=atoi(argv[8]);
	    //   printf("%f %d\n",gWindowSize,windowSnpNumber);
	  }
	else
	  {
	    sparseIndicator=0;
	    windowSnpNumber=atoi(argv[8]);
	  }
      }
    else if(argc==9 && (strcmp("-f",argv[1])==0)) 
      {
	inputFileStyle=2;
	gridSize=atof(argv[8]);
	printf("the -f input still need modification\n");
	printf("Choose to use -c or -h input format!\n");
	exit(0);
      }
    else if(argc==13 && (strcmp("-xpclr",argv[1])==0))
      {
	inputFileStyle=4;
	gridSize=atof(argv[9]);
	cchrn=atoi(argv[10]);
	ccorr=atof(argv[12]);
	if(ccorr>1.0 || ccorr<0.0)
	  perror("Mistake: Correlation level sould be <= 1 or > 0 !\n");
        if(strcmp("-w1",argv[6])==0)
 	 {
	   sparseIndicator=1;
	   gWindowSize=atof(argv[7]);
	   windowSnpNumber=atoi(argv[8]);
	 }
       else
	 {
	   sparseIndicator=0;
	   windowSnpNumber=atoi(argv[8]);
	 }
	if(strcmp("-p1",argv[11])==0)
	  {
	    phaseIndicator=1;
	  }
	else if(strcmp("-p0",argv[11])==0)
	  {
	    phaseIndicator=0;
	  }
	else
	  {
	    perror("Mistake in command line: check -p\n");
	  }
      }
    else if(argc==14 && (strcmp("-xpclr",argv[1])==0))
      {
	inputFileStyle=4;
	gridSize=atof(argv[9]);
	cchrn=atoi(argv[10]);
	ccorr=atof(argv[12]);
	if(ccorr>1.0 || ccorr<0.0)
	  perror("mistake, corr level sould be less than 1 or greater than 0\n");
        if(strcmp("-w1",argv[6])==0)
 	 {
	   sparseIndicator=1;
	   gWindowSize=atof(argv[7]);
	   windowSnpNumber=atoi(argv[8]);
	 }
       else
	 {
	   sparseIndicator=0;
	   windowSnpNumber=atoi(argv[8]);
	 }
	if(strcmp("-p1",argv[11])==0)
	  {
	    phaseIndicator=1;
	  }
	else if(strcmp("-p0",argv[11])==0)
	  {
	    phaseIndicator=0;
	  }
	else
	  {
	    perror("Mistake in command line: check -p\n");
	  }
	setStartPos=atof(argv[13]);
      }
    else if(argc==8 && (strcmp("-c",argv[1])==0))
      {inputFileStyle=3;
       gridSize=atof(argv[7]);
       if(strcmp("-w1",argv[5])==0)
 	 {
	   sparseIndicator=1;
	   gWindowSize=atof(argv[6]);
	   windowSnpNumber=atoi(argv[7]);
	 }
       else
	 {
	   sparseIndicator=0;
	   windowSnpNumber=atoi(argv[7]);
	 }
      }
    else 
      {
      printf("Usage:\n XPCLR -xpclr hapmapInput1 hapmapInput2 mapInput outFile -w gWin(Morgan) snpWin gridSize(bp) chrN -p corrLevel\n");
	//printf("or: XPCLR -c fInputFile outputFile -w1 gWin sWin gridsize(100cM)\n");
	//printf("or: XPCLR -f freqInputFile positionInputFile outputFile -w gWin sWin gridsize(100cM)\n");
	printf("-w1: gWin sets the size of a sliding window(units: 100cM),sWin sets # of SNPs in a window. otherwise, no sliding window\n");
	printf("-p1:the input genotpe is already phased. -p0: the input genotype is not phased\n");
	printf("corrLevel: the value is on (0,1], set corrLevel equal to 0 if no correction is needed\n");
	exit(1);
      }

    char inFile1[100],inFile2[100], inFile3[100],outFile1[150],outFile2[150];
    FILE *fIn1, *fIn2, *fIn3, *fOut1, *fOut2;   

    int i,j,k,l;
    int *n1,*n2,nsegsites;
    double p0;
    char string0[30],string1[30],string2[30],string3[30],string4[30],string5[30],string6[30],string7[30],string8[30],string9[30],string10[30],string11[30],string12[30],string13[30],string14[30],string15[30],string16[30];
    double tempFN1, tempFN2,tempFN3;
    double **f4;
    int **c4;
    int **tc4;
    int *snpList;
    double *wt;
    int snpListLength;
    int idum=100;
    int denseIndicator=0;

    if(sparseIndicator==1)
      {
	snpList=ivector(0,MAXSNPINWINDOW);
	wt=dvector(0,MAXSNPINWINDOW);
      }
    position=dvector(0,MAXSNPLENGTH);
    pposition=dvector(0,MAXSNPLENGTH);
    q1=dvector(0,MAXSNPLENGTH);
    q2=dvector(0,MAXSNPLENGTH);
    x1=ivector(0,MAXSNPLENGTH);
    n1=ivector(0,MAXSNPLENGTH);
    n2=ivector(0,MAXSNPLENGTH);

    p0=1/(2.0*N);

    if(inputFileStyle!=3)
      {
	rawPosition=dvector(0,MAXSNPLENGTH);
	rawPposition=dvector(0,MAXSNPLENGTH);
	rawq1=dvector(0,MAXSNPLENGTH);
	rawq2=dvector(0,MAXSNPLENGTH);
	indexHC=ivector(0,MAXSNPLENGTH);
      }  
    if(inputFileStyle==4)
      {
	strcpy(inFile1,argv[2]);
	strcpy(inFile2,argv[3]);
	strcpy(inFile3,argv[4]);
	strcpy(outFile1,argv[5]);
	strcat(outFile1,".xpclr.txt");
	strcpy(outFile2,argv[5]);
	strcat(outFile2,".xpclr.log");
	fIn1=fopen(inFile1,"r");
	if(fIn1==NULL)
	  perror("Mistake in the input file 1\n");
	fIn2=fopen(inFile2,"r");
	if(fIn2==NULL)
	  perror("Mistake in the input file 2\n");
	fIn3=fopen(inFile3,"r");
	if(fIn3==NULL)
	  perror("Mistake in the input file 3\n");
	fOut1=fopen(outFile1,"w");
	if(fOut1==NULL)
	  perror("Mistake in the output file 1\n");
	/*	fOut2=fopen(outFile2,"w");
	if(fOut2==NULL)
	  perror("mistake in the out2file\n");
	*/
	data1=imatrix(0,MAXHAP,0,MAXSNPLENGTH);
	data2=imatrix(0,MAXHAP,0,MAXSNPLENGTH);
	int count1=0, count2=0;
	char *pp;
	char ss1[]=" ";
	while(fgets(dataString,MAXSNPLENGTH2,fIn1)!=NULL)
	  {
	    count2=0;
	    pp=strtok(dataString," ");
	    data1[count2][count1]=atoi(pp);
	    while((pp=strtok(NULL," "))!=NULL)
	      {
		count2++;
		data1[count2][count1]=atoi(pp);
	      }
	    count1++;
	  }
	snpN=count1;
	hapN1=count2;
	//	printf("snpN1 %d hapN1 %d \n",snpN,hapN1);

	count1=0;
	while(fgets(dataString,MAXSNPLENGTH2,fIn2)!=NULL)
	  {
	    count2=0;
	    pp=strtok(dataString,ss1);
	    data2[count2][count1]=atoi(pp);
	    while((pp=strtok(NULL,ss1))!=NULL)
	      {
		count2++;
		data2[count2][count1]=atoi(pp);
	      }
	    count1++;
	  }
	hapN2=count2;

	count1=0;
	int tempCount1,tempCount2,tempSize;
	fscanf(fIn3,"%s %s %s %s %s %s\n",string1,string6,string2,string3,string4,string5);
	rawPosition[count1]=atof(string2);
	rawPposition[count1]=atof(string3);
	while((fscanf(fIn3,"%s %s %s %s %s %s\n",string1,string6,string2,string3,string4,string5))!=EOF)
	  {
	    count1++;
	    rawPosition[count1]=atof(string2);
	    rawPposition[count1]=atof(string3);
	  
	  }
	for(i=0,tempCount1=0;i<snpN;i++)
	  {
	    tempSize=0;
	    for(j=0;j<hapN2;j++)
	      {
		if(data2[j][i]<2)
		tempSize+=data2[j][i];
	      }
	    if(tempSize>0 && tempSize<hapN2)
	      {
		q2[tempCount1]=((double)tempSize)/((double)hapN2);
		n1[tempCount1]=hapN1;
		n2[tempCount1]=hapN2;
		indexHC[tempCount1]=i;
		tempSize=0;
		for(j=0;j<hapN1;j++)
		  {
		    if(data1[j][i]<2)
		      tempSize+=data1[j][i];
		  }
		x1[tempCount1]=tempSize;
		q1[tempCount1]=((double)x1[tempCount1])/((double)hapN1);
		position[tempCount1]=rawPosition[i];
		pposition[tempCount1]=rawPposition[i];
		tempCount1++;
	      }
	  }
	//	printf("snpN2 %d hapN2 %d \n", snpN,hapN2);
	snpN=tempCount1;

	fclose(fIn1);
	fclose(fIn2);
	fclose(fIn3);
      }
    else if(inputFileStyle==2)
      {
	re=0.01;
        n1[0]=510; //56 
        n2[0]=100; //42
	strcpy(inFile1,argv[2]);
	strcpy(inFile2,argv[3]);
	strcpy(outFile1,argv[4]);
	strcat(outFile1,".compScan.txt");
	strcpy(outFile2,argv[4]);
	strcat(outFile2,".logPCurveScan.txt");
	strcpy(outFile3,argv[4]);
	strcat(outFile3,".screen");
    
	fIn1 = fopen(inFile1,"r");   
	if(fIn1 ==NULL)
	  perror("mistake in the in file1\n");
	fIn2 =fopen(inFile2, "r");
	if(fIn2==NULL)
	  perror("mistake in the in file2\n");
	fOut1=fopen(outFile1,"w");
	if(fOut1==NULL)
	  perror("mistake in the outfile\n");
	fOut2=fopen(outFile2,"w");
	if(fOut2==NULL)
	  perror("mistake in the out2file\n");
	i=0;
	while(fscanf(fIn1,"%s\t%s\t%s\t%s", string1,string2,string3,string4)!=EOF)
	  {
	    rawq1[i]=atof(string1);
	    rawq2[i]=atof(string2);
	    tempFN1=atof(string3);
	    tempFN2=atof(string4);
	    i++;
	  }

	nsegsites=i;

	fscanf(fIn2,"%s\n",string1);
	tempFN1=atof(string1);
	i=1;
	while(fscanf(fIn2,"%s\n",string2)==1)
	  {tempFN2=atof(string2);
	  rawPosition[i]=(tempFN2-tempFN1)*re;
	  i++;
	  }
	rawPosition[0]=0.0;
	if(i!=nsegsites)
	  perror("mistake in readin file\n");
    //pick subset of SNPs.
	for(i=0,l=0;i<nsegsites;i++)
	  {
	    if(/*rawq1[i]>0 && rawq1[i]<1.0 && */rawq2[i]>0.0 && rawq2[i]<1.0)
	      {
		q1[l]=rawq1[i];
		q2[l]=rawq2[i];
		position[l]=rawPosition[i];
		x1[l]=(int)(q1[l]*n1[0]);
		l++;
	      }
	  }
	snpN=l;
	//	printf("snpN: %d\n", snpN);
	fclose(fIn1);
	fclose(fIn2);
   
      }//if inputFileStyle==2
    else if(inputFileStyle==3)
      {
	strcpy(inFile1,argv[2]);
	strcpy(outFile1,argv[3]);
	strcat(outFile1,".clrgs.txt");
	fIn1=fopen(inFile1,"r");
	if(fIn1==NULL)
	  perror("mistake in the in file 1\n");
	fOut1=fopen(outFile1,"w");
	if(fOut1==NULL)
	  perror("mistake in the out file 1\n");
	int tempSize,tempF;
	int count=0;
	while(fscanf(fIn1,"%s %s %s %s %s %s %s",string1,string2,string3,string4,string5,string6,string7)!=EOF)
	  {
	    tempF=atoi(string7);
	    tempSize=atoi(string6);
            if((tempF>0) && (tempF<tempSize))
	      {
		position[count]=atof(string2);
		pposition[count]=atof(string3);
		n1[count]=atoi(string4);
		n2[count]=atoi(string6);
		x1[count]=atoi(string5);
		q1[count]=(double)x1[count]/((double)n1[count]);
		q2[count]=(double)(atoi(string7))/((double)n2[count]);
		count++;
	      }
	   }
	snpN=count;
	fclose(fIn1);
      }
    else
      {printf("mistake in inputFileStyle\n");
      exit(1);
      }

    centroPosition=determineCentromere(pposition,snpN);

    gridNumber=int((pposition[snpN-1]-pposition[0])/gridSize);
    //  printf("gridnumber: %d\n",gridNumber);
    gridPosition=dvector(0,gridNumber+1);
    phGridPosition=dvector(0,gridNumber+1);
    gridPosition[0]=position[0];
    phGridPosition[0]=pposition[0];
    int tempStartt=0,tempEndd=snpN;
    

    for(i=1;i<gridNumber;i++)
      {
	phGridPosition[i]=phGridPosition[i-1]+gridSize;
	for(k=tempStartt;k<tempEndd;k++)
	  {
	    if(pposition[k]>=phGridPosition[i])
	      {
		gridPosition[i]=(phGridPosition[i]-pposition[max((k-1),0)])/(pposition[k]-pposition[max((k-1),0)])*(position[k]-position[max((k-1),0)])+position[max((k-1),0)];
		tempStartt=max(0,(k));
		break;
	      }
	  }
      }


   //estimate global w value.
    w=getW(q2,q1,snpN);
    fstValue=fst(q2,q1,n2[0],n1[0],snpN);
    //    printf("w:\t%lf\tfst:\t%lf\n\n\n", w,fstValue);

    if(sparseIndicator==0)
      { //this part is old code,without setting grid for searching causal mutant positions.
        for(i=0;i<snpN;i++)
         {
          	  printf("%d\t",i);
	  startP=MAX((i-windowSnpNumber/2),0);
	  endP=MAX((i+windowSnpNumber/2),(windowSnpNumber-i));
	  endP=MIN(endP,snpN);
          logPN=logPNeutralityScan(w,p0,q2,x1,n1,startP,endP);
	  logPS=logPSelectionScan(w,p0,x1, q2, position,n1,i,startP,endP,&maxS,sValue,sLength);
	  logPRatio=2.0*(logPS-logPN);
	  fprintf(fOut1,"%d %d %f %f %f %f\n",cchrn,i,pposition[i],position[i],logPRatio,maxS);
       	  printf("%f %f %f %f %f %f\n",pposition[i],position[i],logPS,logPN,logPRatio,maxS);
         }
      }
    else
      {
	if(argc==14 && inputFileStyle==4)
	  {printf("set the start Position: %f\n",setStartPos);
	    tempI=(int)((setStartPos-pposition[0])/((double)gridSize));
	    tempI=MAX(0,tempI);
	  }
	else
	  tempI=0;
	for(i=tempI;i<gridNumber;i++)
	  {
	    printf("\rprocess:%2.2f%%",(double)i/((double)gridNumber)*100);
            denseIndicator= getSnpWindow(position,snpN,snpList,&snpListLength,windowSnpNumber,gridPosition[i],gWindowSize);
	   if(1==denseIndicator)
	     {
              getSubsetSnpWindow(windowSnpNumber,gWindowSize,snpList,&snpListLength,&idum);
	     }  
	   if(inputFileStyle==4)
	     {
	       getWeight(data2,hapN2,indexHC,snpList,snpListLength,wt,ccorr,phaseIndicator);
	       logPN=logPNeutralityScan3(w,p0,q2,x1,n1,snpList,snpListLength,wt);
	       logPS=logPSelectionScan3(w,p0,x1,q2,position,n1,gridPosition[i],&maxS,sValue,sLength,snpList,snpListLength,wt);
	     }
	   else
	     {
	       logPN=logPNeutralityScan2(w,p0,q2,x1,n1,snpList,snpListLength);
	       logPS=logPSelectionScan2(w,p0,x1,q2,position,n1,gridPosition[i],&maxS,sValue,sLength,snpList,snpListLength);
	     }
	    logPRatio=2.0*(logPS-logPN);

	    if((phGridPosition[i]>=(pposition[centroPosition-1])) &&(phGridPosition[i]<=pposition[centroPosition]))
	      {
              logPS=0;
	      logPN=0;
	      logPRatio=0;
	      }
	    fprintf(fOut1,"%d %d %d %f %f %f %f\n",cchrn,i,snpListLength,phGridPosition[i],gridPosition[i],logPRatio,maxS);
	  }
      }
    free_dvector(position,0,MAXSNPLENGTH);
    free_dvector(pposition,0,MAXSNPLENGTH);
    free_dvector(q1,0,MAXSNPLENGTH);
    free_dvector(q2,0,MAXSNPLENGTH);
    free_ivector(x1,0,MAXSNPLENGTH);
    free_ivector(n1,0,MAXSNPLENGTH);
    free_ivector(n2,0,MAXSNPLENGTH);
    free_dvector(gridPosition,0,gridNumber+1);
    free_dvector(phGridPosition,0,gridNumber+1);
    if(sparseIndicator==1)
      {
	free_ivector(snpList,0,MAXSNPINWINDOW);
	free_dvector(wt,0,MAXSNPINWINDOW);
      }

    if(inputFileStyle!=3)
      {
	free_dvector(rawPosition,0,MAXSNPLENGTH);
	free_dvector(rawPposition,0,MAXSNPLENGTH);
	free_dvector(rawq1,0,MAXSNPLENGTH);
	free_dvector(rawq2,0,MAXSNPLENGTH);
	free_ivector(indexHC,0,MAXSNPLENGTH);
	//	fclose(fOut2);
	if(inputFileStyle==4)
	  {
	    free_imatrix(data1,0,MAXHAP,0, MAXSNPLENGTH);
	    free_imatrix(data2,0,MAXHAP,0, MAXSNPLENGTH);
	  }
      }

    fclose(fOut1);
    return 0; 
  }


