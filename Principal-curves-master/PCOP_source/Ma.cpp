// contiene la dimension del espacio original, la matriz y el xo  necesaria para
// transformar los puntos al sistema de coordenadas original


#include "Ma.h"

M_a::M_a(int d,int p,float **M,float *x){
			 Dim = d;
			 profundidad = p;
			 Ma = M;
			 xa = x;
}

M_a::~M_a(){
  int i;
  for (i=0;i<Dim;i++)  free(Ma[i]);
  free(Ma);
  free(xa);
}

float	 *M_a::aplicar_Ma_punt(float *punt){
  float *p2;
  float *p3;


  p2 = Mxv(Ma,punt-profundidad); // habremos creado todos los puntos de Dim+profundidad con las 1eras pos. =0
  p3 = sum_v(p2,xa);

  free(p2);
  return p3;

}



float  *M_a::aplicar_Ma_vect(float *vect){
  float *v3;

  v3 = Mxv(Ma,vect-profundidad); // la profundidad sera maximo 2

  return v3;                    // vector en coordenades originals
}


M_a *M_a::donar_M_a(float **Mbopt,float *xo){
  int i,j,n_prof;
  float **n_Ma, **n_Ma2;
  float *n_xa;

  /* new prof */
  n_prof = profundidad +1;

  /* new Ma */
  n_Ma  = (float**)malloc(Dim*sizeof(float *));
  for (i=0;i<Dim;i++)  n_Ma[i] = (float*)calloc(Dim,sizeof(float));

  for (i=0;i<profundidad;i++) n_Ma[i][i] = 1;

  for (i=0;i<Dim-profundidad;i++)
	for (j=0;j<Dim-profundidad;j++)
		 n_Ma[i+profundidad][j+profundidad] = Mbopt[i][j];

  n_Ma2 = MxM(Ma,n_Ma);

  
  for (i=0;i<Dim;i++) free(n_Ma[i]);
  free(n_Ma);

  /* new xa */
  n_xa = aplicar_Ma_punt(xo);

  return new M_a(Dim,n_prof,n_Ma2,n_xa);   // li pasem la profunditat del subspai
}


//////Private
//// vect ops


float *M_a::Mxv(float **M1,float *v){
// vxM, trabajamos con vectores fila.
 int i,j;
 float sum;
 float *v3 = (float *) malloc(Dim*sizeof(float));

 for(i=0;i<Dim;i++){
	 sum = 0;
	 for(j=0;j<Dim;j++){
		sum += v[j]*M1[j][i];
	 }
	 v3[i] = sum;
 }
 return v3;
}


float **M_a::MxM(float **M1,float **M2){
// vxM, trabajos con vectores fila.
 int i,ii,j;
 float sum;
 float **M3;

 M3 = (float**)malloc(Dim*sizeof(float *));

 for (i=0;i<Dim;i++)
	  M3[i] = (float *)calloc(1,Dim*sizeof(float));

 for(i=0;i<Dim;i++){
  for(ii=0;ii<Dim;ii++){
	 sum = 0;
	 for(j=0;j<Dim;j++){
		sum += M1[i][j]*M2[j][ii];
	 }
	 M3[i][ii] = sum;
  }
 }
 return M3;
}

float *M_a::sum_v (float *v1,float *v2){
 int i;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  v3[i] = v1[i]+v2[i];
 return v3;
}



