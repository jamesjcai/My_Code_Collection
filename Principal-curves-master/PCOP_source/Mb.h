

#include <memory.h>
#include <string.h>

#include "math.h"
#include "stdlib.h"
#include <memory.h>

#include "stdio.h" //###

#include "Ma.h"

#define PI 3.1415926535897932

class M_b{
	 private :

	 int      Dim;
	 float    *xo;
	 float    **Mb;
	 float    **MId;
	 float    **MInv;

	 // vect ops

	 float **inv(float **M);
	 float *Mxv(float **M,float *v);
	 float **MxM(float **M1,float **M2);
	 float *mult_esc(float e,float *v);
	 float  mult_v(float *v1,float *v2);
	 float *sum_v(float *v1,float *v2);
	 float *dif_v(float *v2,float *v1);
	 float *norma_v(float *v);

	 public:

	 // constructoras

	 M_b(int Dim,float *b);   // reusamos el vector b pasado por parametro.
	 M_b(int Dim,float **n_M,float *n_xo);
	 ~M_b();

	 M_b  *girar(int eix,float angle); // s'aplicara un gir a Mb per l'eix o dimensio donat. L'angle d'aquests girs sera el que diferenci un objecte matriu resultant d'un altre
	 M_b  *replicar();                 // fará una copia de la M_b; 
	 void  calcular_la_inversa(); 	// calculamos la inversa  realizados los giros de Mb y antes de aplicarla sobre los puntos. Si el Mb resulta optimo se calculara mas de 1 vez, sino una sola vez
	 M_a   *donar_M_a(M_a *Ma);     // Ma es la matriu que els espais inferiors  necesitaran(una per cada subspai)
									 // per pasar els punts a les coordenades originals
	 void  rebre_xo(float *punt);
	 float *aplicar(float *punt);  // aplica Mb al punt y dona el punt per l'espai inferior
     float *desaplicar(float *punt); // extraer para las coordenadas originales. op inversa a aplicar
	 float *donar_bopt();

};
