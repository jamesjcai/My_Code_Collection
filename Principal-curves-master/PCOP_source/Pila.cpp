#include "stdlib.h"

#include "pila.h"

pila::pila(){top = NULL;}

void pila::apilar (void *pt){
  n_top = (node *)malloc(sizeof(node));

  n_top->pt = pt;
  n_top->seg = top;
  top = n_top;
}

void *pila::desapilar(){
  void  *tp=top->pt;              // no comprovem si la pila es buida, es suposa un control extern

  n_top = top->seg;
  free(top);
  top = n_top;

  return tp;
}

int pila::pila_buida(){
	 return !top;
}





