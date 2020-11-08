

 #include "stdlib.h"

 #include "ll_flt.h"

 ll_flt::ll_flt(){
	Topleft = (node *) malloc(sizeof(node));
	Topright = Topleft;
 }

 ll_flt::~ll_flt(){
	node *act_pt,*aux_pt;

	act_pt = Topleft;
	while (act_pt){
	  aux_pt = act_pt->seg;
	  delete act_pt;
	  act_pt = aux_pt;
	}
 }


 void ll_flt::add(float info){
	Topright->info = info;
	Topright->seg  = (node *)calloc(1,sizeof(node));
	Topright = Topright->seg;
 }

 void ll_flt::addrev(float info){
	newTopleft = (node *) malloc (sizeof(node));
	newTopleft->info = info;
	newTopleft->seg  = Topleft;
	Topleft = newTopleft;
 }

 void ll_flt::resetpt(void **pt){
	 *pt=Topleft;
 }

 void *ll_flt::noend(void *pt){
	return (((node *) pt)->seg);   // a l'ultim sera NULL perque creem amb calloc.
							  // si pt->seg == NULL el camp info de pt estara buit
 }

 float ll_flt::llpt(void *pt){
	return (((node *)pt)->info);
 }

 void ll_flt::advpt(void **pt){
	 *pt = ((node *)*pt)->seg;
 }

 void ll_flt::modpt(void *pt,float info){
	((node*)pt)->info = info;
 }



