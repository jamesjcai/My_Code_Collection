
 #include "stdlib.h"

 #include "ll_pnt.h"


 ll_pnt::ll_pnt(){
	Topleft = (node *) calloc(1,sizeof(node));
	Topright = Topleft;
 }

 ll_pnt::~ll_pnt(){
	node *act_pt,*aux_pt;

	act_pt = Topleft;
	while (act_pt){
	  aux_pt = act_pt->seg;
	  delete act_pt->info;
	  delete act_pt;
	  act_pt = aux_pt;
	}
 }


 void ll_pnt::add(void *info){
	Topright->info = info;
	Topright->seg  = (node *)calloc(1,sizeof(node));
	Topright = Topright->seg;
 }

 void ll_pnt::addrev(void *info){
	newTopleft = (node *) malloc (sizeof(node));
	newTopleft->info = info;
	newTopleft->seg  = Topleft;
	Topleft = newTopleft;
 }

 void ll_pnt::resetpt(void **pt){
	 *pt=Topleft;
 }

 void *ll_pnt::noend(void *pt){
	return (((node *) pt)->seg->seg);   // a l'ultim elem. contindrà info, s'haurà de tractar apart.
							  // si pt->seg == NULL el camp info de pt estara buit
 }

 void *ll_pnt::llpt(void *pt){
	return (((node *)pt)->info);
 }

 void ll_pnt::advpt(void **pt){
	 *pt = ((node *)*pt)->seg;
 }

 void ll_pnt::modpt(void *pt,void *info){
	((node*)pt)->info = info;
 }



