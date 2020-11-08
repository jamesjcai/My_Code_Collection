// contiene la dimension del espacio original, la matriz y el xo  necesaria para
// transformar los puntos al sistema de coordenadas original

#include "stdlib.h"

#include "stdio.h" //###

class M_a{
	private:
	 int 		 Dim;
	 int      profundidad;
	 float	 **Ma;
	 float    *xa;

	 // vect ops

	 float    *Mxv(float **M1,float *v);
	 float    **MxM(float **M1,float **M2);
	 float    *sum_v(float *v1,float *v2);

	public:
	 // constructores
	 M_a(int d,int p,float **M,float *x);
	 ~M_a();

	 float *aplicar_Ma_punt(float *punt);
	 float *aplicar_Ma_vect(float *vect);
	 M_a   *donar_M_a(float **Mbopt,float *xo); // proporciona el M_a  per posar els punts en les cordenades de l'espai inicial, al nou subspai
};
