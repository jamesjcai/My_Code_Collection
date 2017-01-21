
/******************************************************************************
   Subroutines for generating random number, memory allocation etc.
       Hua Chen   2008-6-20
       hchen@genetics.med.harvard.edu

*******************************************************************************/
 
#include "hualib.h"

/******************************************************************************
  
   subroutines for generating random numbers.

*******************************************************************************/

double ran1(int *idum)  {
	int j;
	int k;
	static int iy=0;
	static int iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
	  /*
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
	  */
	  *idum = MAX(-*idum,*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
       	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;

}


 double ran1()
  {
    double rand1;
    // jamescai rand1 =  drand48();
	rand1 = (double(rand()) / RAND_MAX);
    //printf("%f\n", rand1);
    return(rand1);
    //  return((double)rand1/((double)RAND_MAX+1));
  }
 
double randOpen01(int *idum)
  {
    double r;
    double ran1(int *idum);
    do{
      r=ran1(idum);
    }while(r<=0.0 || r>=1.0);
    return r;
  }


/***************************************************************/

#include <math.h>

double binopdf(int n,int k, double p)
{
  double bico(int n,int k);
  double prob;
  if(p>0.0 && p<1.0)
    {
    prob=bico(n,k)*pow(p,(double)k)*pow((1-p),(double)(n-k));
    }
  else if(p==0.0)
    {
      prob=EQUAL(k,0)*1.0;
    }
  else if(p==1.0)
    {
      prob=EQUAL(k,n)*1.0;
    }
  else 
    {
      perror("mistake in binopdf\n");
    }
  return prob;
    
}

double factln(int n)
{
	double gammln(double xx);
	//	void nrerror(char error_text[]);
	static double a[101];

	if (n < 0) 
	  {
	    printf("Negative factorial in routine factln: %d\n",n);
	    exit(0);
	  }
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
}

double bico(int n, int k)
{
	double factln(int n);

	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}


// rate =  1/mu, *idum is the seed for ran1.

double exprnd(double rate, int *idum)
    {
      double rdum,ernd;
      double ran1(int *idum);
      while((rdum=ran1(idum))==1.0);
      ernd=-1.0/(rate)*log(1.0-rdum);
      return ernd;
    }

/*generate rv of gamma distribution */
double gamdev(int ia, int *idum)
   {
     double ran1(int *idum);
     int j;
     double am,e,s,v1,v2,x,y;
     if(ia<1) nrerror((char *)"error in routine gamdev");
     if(ia <6){
       x=1.0;
       for(j=1;j<=ia;j++) x+=ran1(idum);
       x=-log(x);
     }else{
       do{
	 do{
	   do{
	     v1=ran1(idum);
	     v2=2.0*ran1(idum)-1.0;
	   }while (v1*v1+v2*v2>1.0);
	   y=v2/v1;
	   am=ia-1;
	   s=sqrt(2.0*am+1.0);
	   x=s*y+am;
	 }while (x<=0.0);
	 e=(1.0+y*y)*exp(am*log(x/am)-s*y);
       }while (ran1(idum)>e);
     }
     return x;
   }  

/*generate a gamma distribution */
double gamma(double shape, double scale, int *idum)
{
  double ran1(int *idum);
  double gasdev(int *idum);
  if (shape < 1.)
    return gamma(shape + 1., scale, idum) * pow(randOpen01(idum), 1.0 / shape);

  const double d = shape - 1. / 3.;
  const double c = 1. / sqrt(9. * d);
  double x, v;
  for (;;) {
    do {
      x = gasdev(idum);
      v = 1.0 + c * x;
    } while (v <= 0.0);
    v = v * v * v;
    const double u = randOpen01(idum);
    const double x2 = x * x;
    if (u < 1.0 - 0.0331 * x2 * x2)
      return (d * v / scale);
    if (log(u) < 0.5 * x2 + d * (1.0 - v + log(v)))
      return (d * v / scale);
  }

}

/*generate rv of normal distribution N(0,1)*/
double gasdev(int *idum)
   {
     double ran1(int *idum);
     static int iset=0;
     static double gset;
     double fac, rsq, v1, v2;
     if(*idum <0) iset=0;
     if(iset==0){
       do{
	 v1=2.0*ran1(idum)-1.0;
	 v2=2.0*ran1(idum)-1.0;
	 rsq=v1*v1+v2*v2;
       } while(rsq>=1.0 || rsq==0.0);
       fac=sqrt(-2.0*log(rsq)/rsq);
       //use the box-muller transformation
       gset=v1*fac;
       iset=1;
       return v2*fac;
     }
     else
       {
	 iset=0;
	 return gset;
       }
   }


int poisso(double u)
{
	double  cump, ru,p;
	int i=1;

	if( u > 30. ) return( (int)(0.5 + gasdev(u,u)) );
	ru = ran1();
	p = exp(-u);
	if( ru < p) return(0);
	cump = p;
	
	while( ru > ( cump += (p *= u/i ) ) )
		i++;
	return(i);
}


/* a slight modification of crecipes version */
/* copied from ms(Hudson2002) */
double gasdev(double m,double v)
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	double ran1();

	if  (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset= v1*fac;
		iset=1;
		return( m + sqrt(v)*v2*fac);
	} else {
		iset=0;
		return( m + sqrt(v)*gset ) ;
	}
}


/*generate rv of beta distribution. */
double betadev(double a, double b, int *idum)
   {
     double ran1(int *idum);
     double bnldev(double pp, int n, int *idum);    
     double gamma(double shape, double scale, int *idum);
     double g1, g2, r,tg,p;
     g1=gamma(a,1,idum);
     g2=gamma(b,1,idum);
     r=g1/(g1+g2);
     if(g1==0.0 && g2==0.0) tg=1;
     if(tg==1)
       { p=a/(a+b);
       r=bnldev(p,1,idum);
       }
     return r;
   }
 


/* generate rv of binomial distribution.  */    
double bnldev(double pp, int n, int *idum)
      {  int j;
         static int nold = (-1);
         double am,em,g,angle,p,bnl,sq,tt,y;
         static double pold = (-1.0), pc, plog, pclog, en, oldg;
	
         p=(pp <= 0.5 ? pp : 1.0-pp);
         am=n*p;
         if( n< 25) {
           bnl = 0.0;
           for(j=1; j<=n; j++)
              if(ran1(idum)<p) ++bnl;
         }
         else if(am<1.0){
            g=exp(-am);
            tt=1.0;
            for(j=0; j<=n; j++) {
               tt*=ran1(idum);
               if(tt<g)break;
            }
            bnl = (j<= n ? j:n);
         }
         else {
            if(n!=nold) {
               en=n;
               oldg=gammln(en+1.0);
               nold=n;
             }
            if(p!=pold) {
               pc =1.0-p;
               plog=log(p);
               pclog=log(pc);
               pold=p;
            }
	 
            sq=sqrt(2.0*am*pc);
            do {
               do {
                   angle =PI*ran1(idum);
                   y=tan(angle);
                   em=sq*y+am;
                  } while(em <0.0 || em >=(en+1.0));
                em = floor(em);
                tt=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)-gammln(en-em+1.0)
                       + em*plog+(en-em)*pclog);
              //  cout <<"tt:"<<am<<endl;
              } while(ran1(idum) >tt);
               bnl =em;

           }
           if(p!=pp) bnl =n-bnl;
           return bnl;
        }



#define ITMAX 100
#define EPSG 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPSG) break;
	}
	if (i > ITMAX) nrerror((char *)"a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPSG
#undef FPMIN
       
/* take log of gamma function. from NR for C.*/           
 double gammln(double xx)
       {
          double x,y,tmp, ser;
          static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
      }


#define ITMAX 100
#define EPSG 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
	  if (x < 0.0) nrerror((char *)"x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPSG) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror((char *)"a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPSG


/*incomplete gamma function */
double gammp(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror((char *)"Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}


/* the error function */
double erff(double x)
{
	double gammp(double a, double x);

	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}

/* the cdf function of normal distribution */
/*x is the variable, m is mean and s is the sd */

#define SQRT2 1.4142135623730950488016887242097

double gaussiancdf(double x, double m, double s)
  {
    return 0.5+0.5*erff((x-m)/s/SQRT2);
  }


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


#define NRANSI

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) perror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}
#undef NRANSI


#define EPSQ 1.0e-2
#define JMAX 30
#define JMAXP (JMAX+1)
#define K 5

double qromb(double (*func)(double), double a, double b)
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	double trapzd(double (*func)(double), double a, double b, int n);
	
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPSQ*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	perror("Too many steps in routine qromb");
	return 0.0;
}
#undef EPSQ
#undef JMAX
#undef JMAXP
#undef K



#define EPSS 1.0e-4
#define JMAX 30

double qtrap(double (*func)(double),double a,double b)
{
  double trapzd(double (*func)(double),double a,double b, int n);
  int j;
  double s,olds=0.0, temps;
  for(j=6;j<=JMAX;j++)
    {
      s=trapzd(func,a,b,j);
      if(j>5)
	if(fabs(s-olds)<EPSS*fabs(olds) || (s==0.0 && olds==0.0))
	  {return s;
	  }
      
      //   printf("%e %e\n",olds,s);
      temps=olds;
      olds=s;
    }
  perror("too many steps in routine qtrap!\n");
  printf("%e %e \n",temps, s);
  return s;
}



#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC




/*****************************************************************************

            subroutines for allocation dynamic memories.

*****************************************************************************/

double *vector(long nl, long nh)
{ double *v;
 v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
 if(!v) perror("allocation failure in vector");
 return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror((char *)"allocation failure in vector()");
	return v-nl+NR_END;
}


char *cvector(long nl, long nh)
{char *v;
 v=(char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char)));
 if(!v)nrerror((char *)"allocation failure is cvector()");
 return v-nl+NR_END;
}

//allocate an int vector with the subscript range v[nl....nh]
int *ivector(long nl, long nh)
{int *v;
 v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
 if(!v)nrerror((char *)"allocation failure is ivector()");
 return v-nl+NR_END;
}

//allocate a double matrix with subscript range m[nrl...nrh][ncl...nch]
double **matrix(long nrl, long nrh, long ncl, long nch)
{ long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;
  
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror((char *)"allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror((char *)"allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

char **cmatrix(long nrl, long nrh, long ncl, long nch)
{ long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  char **m;
  
	/* allocate pointers to rows */
	m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char*)));
	if (!m) nrerror((char *)"allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
	if (!m[nrl]) nrerror((char *)"allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void matrixvector(double *v1, double **m, double *v2, int n)
    {
      int i, j;

      for (i=0; i<=n; i++)
         for (v1[i]=0.0, j=0; j<=n; j++) v1[i] += m[i][j] * v2[j];
    }
    
 void zerovec(double *v, int n)  {
	int i;
	
	for(i=0; i<=n; i++) v[i] = 0.0;
	}
	
double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror((char *)"allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror((char *)"allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror((char *)"allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror((char *)"allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	double **m;
	/* allocate array of pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror((char *)"allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror((char *)"allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;
	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror((char *)"allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror((char *)"allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror((char *)"allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

/****************************************************************************

*****************************************************************************/

char *itoa(int value, char *digits, int base)
        {
          char *d;
	  unsigned u;
 /*assume unsigned is big enough to hold all the *unsigned values -x
   could possibly be -- don't know how well this assumption hold on the
   DeathStation 9000, so beware of nasal demons.*/

          d = digits;
	  if(base==0)
	    base=10;
	  if(digits==NULL || base<2 || base>36)
	    return NULL;
	  if(value<0){
	    *d++ = '-';
	    u=-((unsigned)value);
	  }else
	    u=value;
	  utoa(u,d,base);
	  return digits;
	}

char *utoa(unsigned value, char *digits, int base)
        { 
          char *s, *p;
	  s=(char *)"0123456789abcdefghijklmnopqrstuvxwyz";
	  if(base==0)
	    base=10;
	  if(digits==NULL ||base<2||base>36)
	    return NULL;
	  if(value<(unsigned)base){
	    digits[0] = s[value];
	    digits[1]='\0';
	  }else {
	    for(p=utoa(value /((unsigned)base),digits,base);
		*p;
		p++);
	    utoa(value % ((unsigned)base), p, base);
	  }
	  return digits;
	}


/*the quick sort */
 void qs(long *items,long *rank, int left, int right)
    {

       register int i,j;
       long x,y;
       i=left, j=right;
       x=items[(int)(left+right)/2];
       do{
            while((items[i]<x) && (i<right)) i++;
            while((x <items[j]) && (j >left)) j--;

            if(i<=j)
              {
                y=items[i];
                items[i]=items[j];
                items[j]=y;

                y=rank[i];
                rank[i]=rank[j];
                rank[j]=y;
                i++; j--;
              }
       }while(i<=j);
       if(left<j) qs(items,rank, left,j);
       if(i<right)qs(items,rank,i,right);

    }




int ranvec(int n,double pbuf[])
{
	int i;
	for(i=0; i<n; i++)
		pbuf[i] = ran1();

	return(0);
}

int order(int n,double pbuf[])
{
        int gap, i, j;
        double temp;

        for( gap= n/2; gap>0; gap /= 2)
           for( i=gap; i<n; i++)
                for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
                   temp = pbuf[j];
                   pbuf[j] = pbuf[j+gap];
                   pbuf[j+gap] = temp;
                   }
return 0;
}


int ordran(int n,double pbuf[])	
{
	ranvec(n,pbuf);
	order(n,pbuf);
	return(0);
}

/* */

int mnmial(int n,int nclass,double p[],int rv[])
{
	double ran1();
	double x, s;
	int i, j;

	for(i=0; i<nclass; i++) rv[i]=0;
	for(i=0; i<n ; i++) {
	   x = ran1();
	   j=0;
	   s = p[0];
	   while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
	   rv[j]++;
	   }
	return(j);
}


/***************************************************************************
     dynamic allocation of memory in C++
**************************************************************************/
/*
template < typename T >
T **Allocate2DArray( int nRows, int nCols)
{
    //(step 1) allocate memory for array of elements of column
    T **ppi = new T*[nRows];

    //(step 2) allocate memory for array of elements of each row
    T *curPtr = new T [nRows * nCols];

    // Now point the pointers in the right place
    for( int i = 0; i < nRows; ++i)
    {
        *(ppi + i) = curPtr;
         curPtr += nCols;
    }
    return ppi;
}


template < typename T >
void Free2DArray(T** Array)
{
    delete [] *Array;
    delete [] Array;
}

*/


/***************************************************************************
                   statistical subroutines.
****************************************************************************/
    /*Given an array of data[1..n], this routine returns its mean ave, 
       average deviation adev,standard deviation sdev, variance var, 
       skewness skew, and kurtosis curt.*/

  void moment(double data[], int n, double *ave, double *adev, double *sdev,
              double *var, double *skew, double *curt)

   {
     int j;
     float ep=0.0,s,p;
     if (n <= 1) nrerror((char *)"n must be at least 2 in moment");
     s=0.0; //First pass to get the mean.
     for (j=1;j<=n;j++) s += data[j];
     *ave=s/n;
     *adev=(*var)=(*skew)=(*curt)=0.0; 
     /*Second pass to get the rst (absolute), second, third, and fourth 
       moments of the deviation from the mean. */
     for (j=1;j<=n;j++) 
       {
         *adev += fabs(s=data[j]-(*ave));
         ep += s;
         *var += (p=s*s);
         *skew += (p *= s);
         *curt += (p *= s);
       }
     *adev /= n;
     *var=(*var-ep*ep/n)/(n-1); //Corrected two-pass formula.
     *sdev=sqrt(*var); //Put the pieces together according to the contional denitions.
     if (*var) 
       { 
         *skew /= (n*(*var)*(*sdev));
         *curt=(*curt)/(n*(*var)*(*var))-3.0;
       } 
     else 
       {
         WARNNING("No skew/kurtosis when variance = 0 (in moment)");
	 *var=0.0;
       }
   }



/*
/*     chi-square test   
void chsone(float bins[], float ebins[], int nbins, int knstrn, float *df,
	float *chsq, float *prob)
{
	float gammq(float a, float x);
	void nrerror(char error_text[]);
	int j;
	float temp;

	*df=nbins-knstrn;
	*chsq=0.0;
	for (j=1;j<=nbins;j++) {
		if (ebins[j] <= 0.0) nrerror("Bad expected number in chsone");
		temp=bins[j]-ebins[j];
		*chsq += temp*temp/ebins[j];
	}
	*prob=gammq(0.5*(*df),0.5*(*chsq));
}


               Kolmogorov-Smirnov Test               
 Given an array data[1..n], and given a user-supplied function of a single
   variable func which is a cumulative distribution function ranging from 0
   (for smallest values of its argument) to 1(for largest values of its
   argument), this routine returns the K¨CS statistic d, and the significance
   level prob. Small values of prob showtha t the cumulative distribution
   function of data is significantly different from func. The array data is
   modified by being sorted into ascending order. 

#include <math.h>
#define NRANSI
#include "nrutil.h"

void ksone(float data[], unsigned long n, float (*func)(float), float *d,
	float *prob)
{
	float probks(float alam);
	void sort(unsigned long n, float arr[]);
	unsigned long j;
	float dt,en,ff,fn,fo=0.0;

	sort(n,data);
	en=n;
	*d=0.0;
	for (j=1;j<=n;j++) {
		fn=j/en;
		ff=(*func)(data[j]);
		dt=FMAX(fabs(fo-ff),fabs(fn-ff));
		if (dt > *d) *d=dt;
		fo=fn;
	}
	en=sqrt(en);
	*prob=probks((en+0.12+0.11/en)*(*d));
}
#undef NRANSI


*/




