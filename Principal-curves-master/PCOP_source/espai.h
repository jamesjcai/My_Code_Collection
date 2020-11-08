#include "ll_pnt.h"
#include "ll_flt.h"
// #include "Ma.h"
#include "Mb.h"
#include "ll_p.h"

#include <memory.h>
#include "math.h"
#include "stdlib.h"
#include "stdio.h" 



//#define PROF_REQ 1  // 0= calculamos curva sobre STV, 1= sobre las curvas de los subespacios,2= sobre las superficies de los subspacios.....



#define C_EPS 0.05   // Sin recursividad, el valor 0.01 fue muy bien
#define NPTMIN 50    // n puntos mínimo para calcular la curva
#define LD   0.5

#define PINZA_MAX PI/4
#define PINZA_MIN PI/((4*NPARTS)+1)  // 90/4 

class espai {
	 private :

	 int   Dim;
	 int   profundidad;
//###	 M_a   *Ma;
	 ll_p  *ll_pt;

     float suma_d;  // suma distancies del mstree.
	 float h_tail;  // distancia al pla.
	 float delta;   // advance
	 float diam;    // longitud de la diagonal del cluster. 
	 float *eps_x;  // eps_x conte les distancies minimes entre el xo del cluster i el xmean del cluster perque aquest sigui valid. Es unic per tots els cluster de l'espai actual.
	 float *xomig;  // l'xomig es necesitara per obtenir_STV() i per obtenir el xo i bo_opt inicial de una corba.
                    // quan calculem la corba en sentit contrari no partirem d'aquest xomig si no del corresponent pop.   

     int bficorba;  // boolean que controla si hi ha aparegut un nou punt en el calcul del nou pop.
	 
	 typedef struct m_d_s{
		 float *xmean;
		 float span;
		 float density;
	 };

	 struct opt{
		 float  VTG;
		 M_b    *Mb;
		 M_b    *Mb_ant;
		 espai  *espai;
		 m_d_s  mds;
	 }optims;

	 struct x{
		 float *act;
		 float *ant;
	 }xo;

	 struct pop{
		 float *alpha;
		 float I;
		 float *b_ast;
		 float var_k;
		 float span;
		 float density;
         espai *espai;         
     };

//	 ll_pnt *ll_pop; ###

	 float Var_PC;
	 float Var_res;
	 float GTV;

     int    dist_al_pla(float *n_punt);
	 void   calcular_htail_delta_xomig_epsx();
	 espai  *obtenir_cluster(M_b *Mb,m_d_s *mds);
	 float  *treure_coord(float *n_pnt);  // fem la projeccio sobre el pla per pasar al subspai de dimensio inferior
	 int    fi_corba(float *n_pnt);
	 int    no_creua_corba(float *pop);
	 void   calcular_Mb(int ejegir,M_b *Mb,float porcion_pinza);
	 float  calcular_corba_en_un_sentit();
	 float  calcular_corba_en_sentit_contrari();
	 float  finalitzacio();       /* retorna la VTG de la corba */
	 float  *allargar(float *bopt);
	 float  obtenir_STV();
	 float  *obtenir_bo_inicial(float *alfa);
     float  Bmst();
	 float  gammln(float xx);
	 float  kernel(float d);
	 int    major(float *v1,float *v2);
     void   eigens(float *A,float *RR,float *E,int N );  //Copyright 1973, 1991 by Stephen L. Moshier
	 
	 // vect ops
	 float  distancia(float *pnt1,float *pnt2);
	 float *mult_esc(float e,float *v);
	 float  mult_v(float *v1,float *v2);
	 float *sum_v(float *v1,float *v2);
	 float *dif_v (float *v2,float *v1);
     float *norma_v(float *v);
	 
	 // modificado 16/4/2002 declaramos las variables de forma statica para que no varien para los diferentes espais
	 static int  PROF_REQ;
	 static int  NPARTS;
	 static float C_H;  
	 static float C_D;     // siempre menor que 0.5
	 // fin

	 public:

	 espai(ll_p *ll_punts,int d,int p);
	 ~espai();
	 float   obtenir_VTG(float **xm);
	 void    rebre_M_a(M_a *n_Ma);      // li pasem el nou Ma al subespai.
     

	 M_a   *Ma;      //###
	 ll_pnt *ll_pop; //###



	 // modificado 16/4/2002 inicializamos las variables
	 // operacions inicialització i extracció d'informació només per el 1er espai.
	 void  inicializar_nparts_ch_cd(int profreq,int nparts,float c_h,float c_d);		
	 void  obtenir_data(FILE *fd_s);
	 // fin

};
