#include "ll_q.h"


ll_q::ll_q(int np){
	npunts = np;
	Topleft = (node *) calloc(1,sizeof(node));// topleft no contendra info
}

ll_q::~ll_q(){    
	node *pt = Topleft;
	node *auxpt;

	while(pt!=NULL){
		auxpt = pt;
		pt = pt->seg;
		free(auxpt);
	}
}

void ll_q::add_ord(float info){
	node *pt=Topleft;
	node *newpt;

	while(pt->seg && pt->seg->info> info){
		  pt = pt->seg;
	}
	newpt = (node *) malloc(sizeof(node));
	newpt->info = info;
	newpt->seg = pt->seg;
	pt->seg = newpt;

}

float ll_q::dmax(){
	int limite,i;
	node  *nd;



	limite = npunts*0.00;  // a un dmax més petit -> calculs més rapids, pero un cluster més restrictiu.

	nd = Topleft->seg;
	for (i=1;i<limite;i++){ 
		nd = nd->seg;
	}
	return (nd->info);

}


