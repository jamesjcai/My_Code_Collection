#include "Ll_p.h"

ll_p::ll_p(int d){
	 int i;  
	 Dim = d;
	 orcluster = ESQUERRA;
	 numcl = 0;
     sum_w = 0; // suma dels pesos dels punts
	 suma_d = 0;

	 topleft = (node *) malloc(sizeof(node));
	 topleft->coord = (float *) calloc(Dim+1,sizeof(float));topleft->coord++;
     topleft->coord[X] = -1*INF; 
	 topright = (node *) malloc(sizeof(node));
	 topright->coord = (float *) calloc(Dim+1,sizeof(float));topright->coord++; // Dim +1 para incluir los pesos
	 topright->coord[X] = INF; 

	 
	 topleft->seg[DRETA] = topright;
	 topleft->noin[DRETA]= (node *)topright;
	 topright->seg[ESQUERRA] = topleft;
	 topright->noin[ESQUERRA]= (node *)topleft;

	 
     topright->seg[DRETA]   = NULL;
	 topright->noin[DRETA] = NULL;
	 topleft->seg[ESQUERRA] = NULL;
	 topleft->noin[ESQUERRA] = NULL;

	 topright->marca = -1; // marca de ja insertats al spanning tree
	 topleft->marca = -1;
	 vn_punts = 0;
	 max = (float *) malloc(Dim*sizeof(float));
	 for(i=0;i<Dim;i++) max[i] = -1*INF; 
	 min = (float *) malloc(Dim*sizeof(float));
	 for(i=0;i<Dim;i++) min[i] = INF; 
	 x_mean = (float *) calloc(Dim,sizeof(float));  //si es una llista de un espai de prof 0, es calculara l'xmean,
                                                    //sino, sera l'origen de coordenades, equivalent al ppp del espai superior.   
}

 ll_p::~ll_p(){
	 node *pt= topleft;
	 node *auxpt;
	 node_satelit *ptst;
	 node_satelit *auxptst;
	 if (pt->seg[DRETA] == pt->noin[DRETA])
        while(pt){ 
  		  auxpt = pt;
		  pt = pt->seg[DRETA];
		  free((auxpt->coord)-1); // per corregir la posicio que deixem lliure al baixar a Dim-1
		  free(auxpt);
        }
	 else	 
	 while(pt){
		ptst = (node_satelit *)pt->noin[DRETA];
		while(ptst){
			 auxptst = ptst;
			 ptst = ptst->seg;
			 free(auxptst);
		}
		ptst = (node_satelit *)pt->noin[ESQUERRA];
		while(ptst){
			 auxptst = ptst;
			 ptst = ptst->seg;
			 free(auxptst);
		}
		auxpt = pt;
		pt = pt->seg[DRETA];
		free((auxpt->coord)-1); // per corregir la posicio que deixem lliure al baixar a Dim-1
		free(auxpt);
	 }

//	free (max);   // les eliminarán l'espai a qui li pasem, quan les hagi utilitzat.
//	free (min);
 }


void ll_p::add_ordX_principal(float *vect){
	node *act_node,*new_node;
	float *x_mean02,*x_mean03;

	int i;
    
	act_node = topright->seg[ESQUERRA];


	if (vect[X] >topleft->seg[DRETA]->coord[X]+(0.5*(act_node->coord[X]-topleft->seg[DRETA]->coord[X]))){   // si es mes gran que el punt mig començem pel canto dret
		
		while (vect[X]<act_node->coord[X]) {
		  act_node = act_node->seg[ESQUERRA];
		}
		new_node =  (node *) malloc(sizeof(node));
		new_node->coord = vect;
		new_node->marca = 0;
		new_node->seg[ESQUERRA] = act_node;
		new_node->noin[ESQUERRA] = (node *)act_node;
		new_node->seg[DRETA] = act_node->seg[DRETA];
		new_node->noin[DRETA] = (node *)act_node->noin[DRETA];
		act_node->seg[DRETA]->seg[ESQUERRA] = new_node;
		act_node->seg[DRETA]->noin[ESQUERRA] = (node *) new_node;
		act_node->seg[DRETA] = new_node;
		act_node->noin[DRETA] = (node *)new_node;
		vn_punts++;		
	}
	else {
		act_node = topleft->seg[DRETA];
		while (vect[X]>act_node->coord[X]) {
		  act_node = act_node->seg[DRETA];
		}
		new_node = (node *) malloc(sizeof(node));
		new_node->coord = vect;
		new_node->marca = 0;
		new_node->seg[DRETA] = act_node;
		new_node->noin[DRETA] = (node *)act_node;
		new_node->seg[ESQUERRA] = act_node->seg[ESQUERRA];
		new_node->noin[ESQUERRA] = (node *)act_node->noin[ESQUERRA];
		act_node->seg[ESQUERRA]->seg[DRETA] = new_node;
		act_node->seg[ESQUERRA]->noin[DRETA] = (node *) new_node;
		act_node->seg[ESQUERRA] = new_node;
		act_node->noin[ESQUERRA] = (node *)new_node;
		vn_punts++;
	}

	 /* buscar max , min, xmean */  // com que ll_p es tot un cluster serán tb el max i mix del cluster
	for (i=0;i<Dim;i++){
		  if (max[i]< vect[i]) max[i] = vect[i];
		  else if (min[i]>vect[i]) min[i] = vect[i];
	}
   /* calcular xmean */
   sum_w += *(vect-1);
   x_mean02 = mult_esc(*(vect-1),vect);
   x_mean03 = sum_v(x_mean,x_mean02);
   delete x_mean; 
   delete x_mean02; 
   x_mean = x_mean03;

}

void ll_p::donar_max_min_xomig(float **mx, float **mn, float **xm,float *s_d){

	 *mx = max;
	 *mn = min;
	 *xm = x_mean;     // punt de ll_p més proper al xmean teóric. xmean ponderat en cas de pertanyer a un spai final.
     *s_d = suma_d;
}


void ll_p::calcular_max_min_cluster(){
 node_satelit *s_act;
 int i;

 numcl++;  // per diferenciar-ho del NULL

 for (i=0;i<Dim;i++)
	 if (max[i]< xorig->coord[i]) max[i] = xorig->coord[i];
	 else if (min[i]> xorig->coord[i]) min[i] = xorig->coord[i];

 for(orcluster=ESQUERRA;orcluster<=DRETA;orcluster++){
  s_act = (node_satelit *)xorig->noin[orcluster];
  
	  for (i=0;i<Dim;i++)
	  if (max[i]< s_act->ptnode->coord[i]) max[i] = s_act->ptnode->coord[i];
	  else if (min[i]>s_act->ptnode->coord[i]) min[i] = s_act->ptnode->coord[i];
  s_act->ptnode->marca = numcl;
  p_n.apilar(s_act->ptnode);

  while (s_act->seg){
	  s_act = s_act->seg;
	  if (s_act->ptnode->marca != numcl){
        for (i=0;i<Dim;i++)
	       if (max[i]< s_act->ptnode->coord[i]) max[i] = s_act->ptnode->coord[i];
	       else if (min[i]>s_act->ptnode->coord[i]) min[i] = s_act->ptnode->coord[i];
        s_act->ptnode->marca = numcl;
		p_n.apilar(s_act->ptnode);
	  }
  }	  

  while(!p_n.pila_buida()){ 
    while (!p_n.pila_buida() && 
 	      !(s_act = ((node_satelit *)((node*) p_n.desapilar())->noin[orcluster]))); // hem de considerar punts situats al limit del cluster sobre la coordenada X	  
	if (s_act){
	  if (s_act->ptnode->marca != numcl){
        for (i=0;i<Dim;i++)
           if (max[i]< s_act->ptnode->coord[i]) max[i] = s_act->ptnode->coord[i];
           else if (min[i]>s_act->ptnode->coord[i]) min[i] = s_act->ptnode->coord[i];
        s_act->ptnode->marca = numcl;
		p_n.apilar(s_act->ptnode);
	  }
	  while (s_act->seg){
		  s_act = s_act->seg;
		  if (s_act->ptnode->marca != numcl){
	         for (i=0;i<Dim;i++)
	           if (max[i]< s_act->ptnode->coord[i]) max[i] = s_act->ptnode->coord[i];
	           else if (min[i]>s_act->ptnode->coord[i]) min[i] = s_act->ptnode->coord[i];
             s_act->ptnode->marca = numcl;        
			 p_n.apilar(s_act->ptnode);
		  }
      }
    }
  }
  
 } 

}


float ll_p::inicialitzacio_principal(){
	ll_q *quartiles;
    float *auxx_mean;
    
	auxx_mean = x_mean;       // x_mean teoric.
	x_mean = mult_esc(1.0/sum_w,auxx_mean);
	delete auxx_mean;
	quartiles = new ll_q(vn_punts); // li pasem el nº de punts de l'espai
	obtener_quartiles(quartiles);   // calculant el min. spanning tree
	dmax = quartiles->dmax();
    x_mean =obtener_satelites();   // nos basamos en dmax para obtener los satelites. Tb buscaremos el xorig.
	delete quartiles;
	return dmax;       // solo necesitamos dmax para separar las diferentes curbas
}

void ll_p::inicialitzacio_final(){
    /* no fem ni spanning tree, ni satelits, donç no calcularem cap corba per aquest spai */
    float *auxx_mean;
    
	auxx_mean = x_mean;       // x_mean teoric.
	x_mean = mult_esc(1.0/sum_w,auxx_mean);  //el xmean el dividim per la suma de pesos dels punts.
	delete auxx_mean;

}

void ll_p::mstinsertar(node *pt){
	  pt->marca = -1; // després tornarem a inicialitzar a 0 a calcular_satelits() per a candidat_visitat().

	  ((node *)pt->noin[ESQUERRA])->noin[DRETA] = ((node *)pt->noin[DRETA]);
	  ((node *)pt->noin[DRETA])->noin[ESQUERRA] = ((node *)pt->noin[ESQUERRA]);

 // netegem per satelits

	  pt->noin[DRETA]=NULL;
	  pt->noin[ESQUERRA]=NULL;
}

int ll_p::mstinsertat(node *pt){
  return pt->marca;
}


void ll_p::obtener_quartiles(ll_q *llqt){  // calculamos el minium spanning tree y guardamos las distancias para obtener los quartiles
	 int   aux_npunts = vn_punts -1;
	 node  *xact;
 	 node  *xpost;
	 node  *xpt;
	 node  *nxtnomstin; 
	 float dpost,dant,d; /* dant : distancia min. a un nodo del tree */
								/* dpost: distancia min. a un nodo fuera del tree */
	 xpt = topleft->seg[DRETA];


     /* preproceso */ // localizamos el segundo punto a insertar en el mst

     xact = xpt->seg[DRETA];
	 xpost = (node *) xpt->noin[DRETA];     /* para el caso : Xdreta = dpost> Xesquerra */
	 dpost = distancia(xpost->coord,xpt->coord);/* cuando xact=topright dpost = Xact-Xpt */

	 while( xact->coord[X] - xpt->coord[X]<dpost ){
		  if ((d= distancia(xact->coord,xpt->coord))<= dpost)                                                          
				dpost = d; 
	      xact = xact->seg[DRETA];   // tots els punts son noinmst
	 }
	  mstinsertar(xpt);   // insertem al mst el xpost anterior, per no desvincular-lo de la llista noin
	  xpt = xpost;        // xpost será el siguiente punto ha insertar   
	  dant = dpost;
	  aux_npunts--;       // restamos uno de más porque la última insercion en ll_q se realiza en el postproceso

	 
     /* cuerpo principal */ 

  while (aux_npunts){
         
	  xact = xpt->seg[DRETA];
	  nxtnomstin = (node *) xpt->noin[DRETA]; 
	  xpost = nxtnomstin;                  /* para el caso : Xdreta = dpost> Xesquerra */
	  dpost = distancia(xpost->coord,xpt->coord);/* cuando xact=topright dpost = Xact-Xpt */

	 /* tratamiento por la derecha */

	 while (xact->coord[X]-xpt->coord[X]<dant && xact->coord[X]-xpt->coord[X]<dpost){
		 if(!mstinsertat(xact)){
			 if ((d= distancia(xact->coord,xpt->coord))< dpost)
				{ dpost = d; xpost = xact;}
			 nxtnomstin = (node *) xact->noin[DRETA];
         }  
		 else if ((d= distancia(xact->coord,xpt->coord))< dant) 
			      dant = d;
		 xact = xact->seg[DRETA];
	 }     

	 if (xact->coord[X]- xpt->coord[X] >= dant){// hemos acotado el dant mín
		 xact = nxtnomstin;  // ens coloquem a l'ultim no insertat
		 while( xact->coord[X] - xpt->coord[X]<dpost && xact!=topright){ // si tb hemos acotado el dpost min, saldremos directamente
			  if ((d= distancia(xact->coord,xpt->coord))<= dpost) 
					{dpost = d; xpost = xact;}
		      xact = (node *)xact->noin[DRETA];
		 }
	 }
	 else{  // hemos acotado el dpost mín
		 while(xact->coord[X]-xpt->coord[X]<dant) {
			  if((mstinsertat(xact)) &&((d= distancia(xact->coord,xpt->coord))< dant))
				dant = d;
			  xact = xact->seg[DRETA];
		 }
	 }

	 /* tratamiento por la izquierda */
	 xact = xpt->seg[ESQUERRA];
     nxtnomstin = (node *) xpt->noin[ESQUERRA]; 
	 while (xpt->coord[X]-xact->coord[X]<dant && xpt->coord[X]-xact->coord[X]<dpost ){
		 if(!mstinsertat(xact)){
			 if ((d= distancia(xact->coord,xpt->coord))< dpost)
				 { dpost = d; xpost = xact;}
			 nxtnomstin = (node *) xact->noin[ESQUERRA];
         }
		 else if ((d= distancia(xact->coord,xpt->coord))< dant) 
			dant = d;
		 xact = xact->seg[ESQUERRA];
	 }

	 if (xpt->coord[X]-xact->coord[X] >= dant){// hemos acotado el dant mín
		 xact = nxtnomstin;   // ens coloquem a l'ultim no insertat
		 while( xpt->coord[X] - xact->coord[X]<dpost ){ // si tb hemos acotado el dpost min, saldremos directamente
			  if ((d= distancia(xact->coord,xpt->coord))< dpost)                                                          
				  {dpost = d; xpost = xact;}
			  xact = (node *)xact->noin[ESQUERRA];
		 }
	 }
	 else {  // hemos acotado el dpost mín
		 while(xpt->coord[X]-xact->coord[X]<dant) {
			  if((mstinsertat(xact)) &&((d= distancia(xact->coord,xpt->coord))< dant))
				dant = d;
			  xact = xact->seg[ESQUERRA];
		 }
	 }

	  /*tratamiento final de iteración */

	  llqt->add_ord(dant); suma_d += dant;
	  mstinsertar(xpt);   // insertem al mst el xpost anterior, per no desvincular-lo de la llista noin
	  xpt = xpost;
	  dant = dpost;
	  aux_npunts--;
  }


      /* postproceso */  // compruevo que no haya un dant mas optimo para el ultimo xpost a insertar en el mst

   xact = xpt->seg[DRETA];
   while(xact->coord[X]-xpt->coord[X]<dant) {
     if((d= distancia(xact->coord,xpt->coord))< dant)  // tots els punts estan mstinsertats
		dant = d;
     xact = xact->seg[DRETA];
   }

   xact = xpt->seg[ESQUERRA];
   while(xpt->coord[X]-xact->coord[X]<dant) {
	  if((d= distancia(xact->coord,xpt->coord))< dant) // tots els punts estan mstinsertats
		  dant = d;
	  xact = xact->seg[ESQUERRA];
   }

	  

   llqt->add_ord(dant); suma_d += dant;
   // netegem per satelits
   xpt->noin[DRETA]=NULL;
   xpt->noin[ESQUERRA]=NULL;
   
   topleft->noin[DRETA]=NULL;
   topright->noin[ESQUERRA]=NULL; 

}


float *ll_p::obtener_satelites(){
  float d,dpost = INF;
  node *xpt = topleft->seg[DRETA];
  node *xact = xpt->seg[DRETA];
  
  while (xpt->seg[DRETA]){
	  while(xact->coord[X]- xpt->coord[X]<dmax){
			if (distancia(xact->coord,xpt->coord)<dmax){
				add_satelit(DRETA,xpt,xact);add_satelit(ESQUERRA,xact,xpt);
                if ((xpt->noin[ESQUERRA]) && ((d=distancia(xpt->coord,x_mean))<dpost))
					{dpost = d;xorig = xpt;}; // busquem xorig. Ha de tindre satelits a dreta i esquerra.
			};
			xact = xact->seg[DRETA];
	  }

	  xpt = xpt->seg[DRETA];
	  xact = xpt->seg[DRETA];
  }
  xoant = xorig;
  delete x_mean;     // borrem l'antic xmean teoric i retornem el xorig practic

  return xorig->coord;
}

void ll_p::add_satelit(int or,node *ptor,node *ptdsti){
	  node_satelit *new_nst;

	  new_nst = (node_satelit *) malloc(sizeof(node_satelit));
	  new_nst->ptnode = ptdsti;
	  new_nst->seg = (node_satelit *)ptor->noin[or];   // s'afegeixen per devant
	  ptor->noin[or] = new_nst;
}

void ll_p::tornar_a_xomig(){
	xoant = xorig;
}

void ll_p::trobar_primer_candidat_clt(float *xm){
  float d,dtop;
  node *xtop;
  node *xodmax;

  orcluster = ( xm[X] > xoant->coord[X]); //true ->DRETA, false ->ESQUERRA  
  while ( fabs(xm[X]-xoant->coord[X]) > dmax)      
	  xoant = xoant->seg[orcluster];  
  xodmax = xoant;   

  while (xoant->marca <1) xoant = xoant->seg[orcluster]; // no exigimos conectado por derecha y izquierda. solo estar conectado
  dtop = distancia(xm,xoant->coord);
  xtop = xoant;
  while (fabs(xoant->coord[X]-xm[X]) < dtop){
      xoant = xoant->seg[orcluster];
	  if(((d = distancia(xm,xoant->coord))<dtop)&& (xoant->marca>0)){ // elim.(xoant->marca>0)
		  dtop =d; xtop = xoant;
      }
  }
  if ( dtop > dmax ) {   // per si hi ha dos clusters a dist>dmax de xm
      xoant = xodmax;
	  orcluster = (orcluster+1) % 2;
      while (fabs(xoant->coord[X]-xm[X]) < dtop){
        xoant = xoant->seg[orcluster];
	    if(((d = distancia(xm,xoant->coord))<dtop)&& (xoant->marca>0)){ // elim.(xoant->marca>0)
	  	  dtop =d; xtop = xoant;
		}
	  }
  }
  xoant = xtop;
}

float *ll_p::primer_candidat_clt(){
  numcl++; 

  if(!(candidat = (node_satelit *) xoant->noin[orcluster])){
      orcluster = (orcluster+1) % 2;
      candidat = (node_satelit *) xoant->noin[orcluster]; // como minim tindra un satelit, sino no pertenyeria al cluster.
  }
  semilla = xoant;  
  return candidat->ptnode->coord;  //### no inserta el primer satelite aunque sea bueno e inserta 2 veces xoant
}

float *ll_p::seguent_candidat_clt(int validacio){

	candidat->ptnode->marca = numcl;
	if (validacio)                          // el objeto espacio nos ha validado el candidato al cluster
			 p_n.apilar(candidat->ptnode);
    /* solo busca puntos dentro del cluster de candidat y en el sentido orcluster sobre la coordenada X*/
	do{
	  if (candidat->seg)  candidat = candidat->seg;
	  else if (semilla){ // podemos venir de tratar la 2nda dir. o  la 1ra.
         orcluster = (orcluster + 1) % 2; 
		 if (!(candidat = (node_satelit *)((node*) semilla)->noin[orcluster])){
 		   do{ 
              if ( p_n.pila_buida()) return NULL;
		      else semilla = (node*) p_n.desapilar();
              if (!(candidat = (node_satelit *)((node*) semilla)->noin[orcluster])){
                  orcluster = (orcluster + 1) % 2; 
			      candidat = (node_satelit *)((node*) semilla)->noin[orcluster];
			      semilla = NULL;
			  }
		   }while(!candidat);  // si el candidat en la 1era dir. es bueno, realizamos la evaluación de candidat 2 veces.
		 }else semilla = NULL;  // nos queda por tratar la 2nda dir.
      }
	  else
		do{ 
           if ( p_n.pila_buida()) return NULL;
		   else semilla = (node*) p_n.desapilar();
           if (!(candidat = (node_satelit *)((node*) semilla)->noin[orcluster])){
                orcluster = (orcluster + 1) % 2; 
			    candidat = (node_satelit *)((node*) semilla)->noin[orcluster];
			    semilla = NULL;
		   }
		}while(!candidat);  // si el candidat en la 1era dir. es bueno, realizamos la evaluación de candidat 2 veces.
	}while (candidat->ptnode->marca == numcl);

	return candidat->ptnode->coord;
}

float *ll_p::canviar_orientacio_clt(){
	 orcluster = (orcluster+1) % 2;    // canviem de sentit, no cambiem xant
//	 candidat = (node_satelit *)xoant->noin[orcluster];  // agafem el primer satelit en sentit contrari al 1er candidat
//     semilla = xoant;
//     if (!candidat) return NULL; // puede que el xoant este al extremo del cluster en la coordenada X
//	 if (candidat->ptnode->marca == numcl) return NULL; // miramos los 2 sentidos por semilla
//	 return candidat->ptnode->coord; // com es en sentit contrari no cal comprobar que estigui previament insertat
     return NULL;
}


int ll_p::n_punts(){
  return vn_punts;
}


///ops lectura llista /////
// necesaris per obtenir la STV final en arribar a la profunditat maxima ( no es poden separar de la classe general)

 void ll_p::resetpt(void **pt){
	*pt = topleft->seg[DRETA];
 }

 void *ll_p::noend(void *pt){
	return (((node *)pt)->seg[DRETA]);      // si es NULL el pt es el centinella
 }

 float *ll_p::llpt(void *pt){
	return (((node *)pt)->coord);
 }

 void ll_p::advpt(void **pt){
	 *pt = ((node *)*pt)->seg[DRETA];
 }


// ops necesaries per la operacio principal

 int  ll_p::llptmarca(void *pt){
	 return ((node *)pt)->marca;
 }

 void ll_p::modpt(void *pt,int info){  // modptmarca
	 ((node *)pt)->marca = info;
 }

 void ll_p::revresetpt(void **pt){
	*pt = topright->seg[ESQUERRA];
 }

 void *ll_p::revnoend(void *pt){
	return (((node *)pt)->seg[ESQUERRA]);   // si es NULL el pt es el centinella
 }

 void ll_p::advrevpt(void **pt){
	 *pt = ((node *)*pt)->seg[ESQUERRA];
 }


////ops_vect///////

float  ll_p::distancia(float *pnt1,float *pnt2){
  int i;
  float sum =0.0;
  for(i =0;i<Dim;i++)
		sum += pow(pnt1[i]-pnt2[i],2);
  return sqrt(sum);
}

float *ll_p::mult_esc(float e,float *v){
 int i;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  v3[i] = v[i]*e;
 return v3;
}

float *ll_p::sum_v (float *v1,float *v2){
 int i;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  v3[i] = v1[i]+v2[i];
 return v3;
}
