

#include "Mb.h"


M_b::M_b(int d,float *b){
  int i;
  int j;
  float *v_dif,*v_acum1,*v_acum2,*mesc;

  Dim =d;
  /* creamos e inicializamos Mb y MId */
  Mb    = (float**)malloc(Dim*sizeof(float *));
  MId   = (float**)malloc(Dim*sizeof(float *));
  MInv  = NULL;

  for (i=0;i<Dim;i++){
	  Mb[i] = (float*)calloc(Dim,sizeof(float));
	  MId[i] = (float*)calloc(Dim,sizeof(float));
  }
  for (i=0;i<Dim;i++){
	  Mb[i][i] = 1;
	  MId[i][i] = 1;
  }

  /* insertamos b_act en matriz identidad */

  i = 0;
  while(!b[i]) i++;
  if (i){
    for(j=i-1;j;j--)Mb[j+1]=Mb[j]; 
    for(j=Dim-2;j>i+1;j--)Mb[j+1]=Mb[j];
  }  

  Mb[0] = b;   // reusamos el vector b pasado por parametro.


  /* proceso de ortonormalización de Gram-Schmidt */

  v_acum1 = (float *)calloc(Dim,sizeof(float));
  Mb[0] = norma_v(Mb[0]);
  for (i=1;i<Dim;i++){
	  delete v_acum1;
      v_acum1 = (float *)calloc(Dim,sizeof(float)); 
	  for (j=0;j<i;j++){
		  mesc = mult_esc(mult_v(Mb[i],Mb[j]),Mb[j]);
		  v_acum2 = sum_v(v_acum1,mesc);
		  delete v_acum1;
		  v_acum1 = v_acum2;
      }
	  v_dif = dif_v(Mb[i],v_acum1);
	  Mb[i]= norma_v(v_dif);
  }
   
}

M_b::M_b(int d,float **n_M,float *n_xo){
  int i;
  Dim = d;
  Mb = n_M;
  xo = n_xo;

  /* creamos e inicializamos  MId */
  MId   = (float**)malloc(Dim*sizeof(float *));
  MInv  = NULL;

  for (i=0;i<Dim;i++){
	  MId[i] = (float*)calloc(1,Dim*sizeof(float));
  }
  for (i=0;i<Dim;i++){
	  MId[i][i] = 1;
  }

}

M_b::~M_b(){
  int i;

  if (MInv)
	 for (i=0;i<Dim;i++){
	  free(Mb[i]);
	  free(MInv[i]);
	  free(MId[i]);
	 }
  else
	 for (i=0;i<Dim;i++){
	  free(Mb[i]);
	  free(MId[i]);
	 }

  free(Mb);
  free(MInv);
  free(MId);

  // free(xo) xo  es gestiona desde espai (calcular_corba_en_sentit(), calcular_corba_en_sentit_contrari()).
}


M_b *M_b::girar(int eix,float angle){
	// 0 < eix < Dim                    giro bo i  eix en sentit positiu
 
 float **Mbaux;

	// insertamos el giro en la matriz identidad
	MId[0][0] = cos(angle);
	MId[0][eix] = sin(angle);
	MId[eix][eix] = cos(angle);
	MId[eix][0] = -1*sin(angle);

	Mbaux = MxM(Mb,MId);

	// corregimos la matriz identidad
	MId[0][0]  = 1;
	MId[0][eix]= 0;
	MId[eix][eix] = 1;
	MId[eix][0]= 0;

	return new M_b(Dim,Mbaux,xo);
}

M_b *M_b::replicar(){
	float **c_Mb;
    int i,j;

    c_Mb = (float **) malloc(Dim*sizeof(float *));
    for (i=0;i<Dim;i++)
	  c_Mb[i] = (float *)malloc(Dim*sizeof(float));
  
	
	  for (i=0;i<Dim;i++) memmove(c_Mb[i],Mb[i],Dim*sizeof(float));
	//for (j=0;j<Dim;j++) c_Mb[i][j] = Mb[i][j];  
      
    
   return new M_b(Dim,c_Mb,xo);
}


void M_b::calcular_la_inversa(){
	int i,j;
	float **c_Mb;
	// calculamos la inversa  realizados los giros de Mb y antes de aplicarla sobre puntos. Si el Mb resulta optimo se calculara mas de 1 vez

	if (MInv){
	  for (i=0;i<Dim;i++)  free(MInv[i]);
	  free(MInv);
	}

	/* copiar la Mb  */ // ya que sera modificada por inv()

    c_Mb = (float **) malloc(Dim*sizeof(float *));
    for (i=0;i<Dim;i++)
	  c_Mb[i] = (float *)malloc(Dim*sizeof(float));
  
     
	for (i=0;i<Dim;i++)  memmove(c_Mb[i],Mb[i],Dim*sizeof(float));
	//for (j=0;j<Dim;j++) c_Mb[i][j] = Mb[i][j];  

	/* calcular inversa */

	MInv = inv(c_Mb);

	/* borrar c_Mb  */
	for (i=0;i<Dim;i++)  free(c_Mb[i]);
	free(c_Mb);
}


float  *M_b::aplicar(float *punt){  /* aplica Mb al punt y dona el punt per l'espai inferior */
  float *p2;
  float *p3;
  p2 = dif_v(punt,xo);
  p3 = Mxv(MInv,p2);     /* se pasa la inversa */

// la primera coordenada sera la distancia al pla, les seguents, les del punt
// al pla inferior.

  free(p2);
  return p3;
}

float  *M_b::desaplicar(float *punt){  /* operación inversa a aplicar */
  float *p2;
  float *p3;
  
  p2 = Mxv(Mb,punt);     
  p3 = sum_v(p2,xo);

// la primera coordenada sera la distancia al pla, les seguents, les del punt
// al pla inferior.

  free(p2);
  return p3;
}


void M_b::rebre_xo(float *punt){
// se ejecutara una vez por avance_cluster() y todas las veces que el pop candidato cruce el hiperplano del pop anterior.
//free xo;    //  no fa falta alliverar l'espai. O bé ha estat eliminat al crear el nou xo, o es un xo compartir amb la matriu del pop anterior.
  xo = punt;  

}

	 /* Ma es la matriu que l'espai inferior  necesacitará per pasar*/
	 /* els seus punts i vectors finals a les coordenades originals    */

M_a *M_b::donar_M_a(M_a *Ma){
// se ejecutara como máximo una vez si resulta ser el M_b optimo.

  return Ma->donar_M_a(Mb,xo);  // xo, en aquest punt haurá de ser de tamany Dim+prof.

}


float *M_b::donar_bopt(){
	return Mb[0];
}


// PRIVATE //////////////////////////////////////////////////////////////////

float **M_b::inv(float **M){
  int i,j,ii,aux;
  float **Inv,Mji;

  Inv = (float**)malloc(Dim*sizeof(float *));
  for (i=0;i<Dim;i++)
	  Inv[i] = (float *)calloc(Dim,sizeof(float));
  for (i=0;i<Dim;i++)
	  Inv[i][i] = 1;


  for (i=0;i<Dim;i++){
	 for(j=(i+1)%Dim;j!=i;j=(j+1)%Dim){
		Mji = M[j][i];
		for(ii=0;ii<Dim;ii++){
			Inv[j][ii]= Inv[j][ii]*M[i][i]-Inv[i][ii]*Mji;
			M[j][ii]= M[j][ii]*M[i][i]-M[i][ii]*Mji;
		}
	 }
  }


/*  for (i=0;i<Dim;i++){
	 j=(i+1)%Dim;
	 aux = j;            // si no se vuelca la  j, el % no funciona correctamente borlan c++ !!¿¿???
	 while (j!=i){
		Mji = M[j][i];
		for(ii=0;ii<Dim;ii++){
			Inv[j][ii]= (Inv[j][ii]*M[i][i])-(Inv[i][ii]*Mji);   // no pot modificarse M avans de Inv
			M[j][ii]= (M[j][ii]*M[i][i])-(M[i][ii]*Mji);
		}
		j=(j+1)%Dim;
	 }
  }
*/
  for (j=0;j<Dim;j++){
	 for(ii=0;ii<Dim;ii++)
		  Inv[j][ii]= Inv[j][ii]/M[j][j];
  }
  return Inv;
}

float *M_b::Mxv(float **M1,float *v){
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


float **M_b::MxM(float **M1,float **M2){
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

float *M_b::mult_esc(float e,float *v){
 int i;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  v3[i] = v[i]*e;
 return v3;
}

float M_b::mult_v(float *v1,float *v2){
 int i;
 float sum = 0.0;
 for(i=0;i<Dim;i++){
		sum += v1[i]* v2[i];
 }
 return sum;
}

float *M_b::sum_v (float *v1,float *v2){
 int i;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  v3[i] = v1[i]+v2[i];
 return v3;
}

float *M_b::dif_v (float *v2,float *v1){
 int i;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  v3[i] = v2[i]-v1[i];
 return v3;
}

float *M_b::norma_v(float *v){
 // devuelve el vector normalizado
 int i;
 float nrm =0.0;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  nrm += pow(v[i],2);
 nrm = sqrt(nrm);
 for(i=0;i<Dim;i++)  v3[i] = v[i]/nrm;
 return v3;
}





