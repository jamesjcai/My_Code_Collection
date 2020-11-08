// ### r: cluster rectangular. 
// ### b: bopt inicial último bopt.
// ### a: un pop no puede retroceder.

#include "espai.h"




espai::espai(ll_p *ll_punts,int d,int p){
	ll_pt = ll_punts;
	Dim = d;
	profundidad = p;
    Ma = NULL;
	ll_pop = NULL;

}

espai::~espai(){

	delete Ma;

//  delete ll_pt;           s'esborra en la finalitzacio tant del cas directe com del recursiu
//	free(xomig);            s'esborra en la finalitzacio tant del cas directe com del recursiu       
//	free(eps_x);            s'esborra en la finalitzacio tant del cas directe com del recursiu
//	free(mds.xmean);        s'esborra a construir_corba_sentit_contrari()
//	free(optims.mds.xmean); s'esborra a construir_corba_sentit_contrari()
//	delete (optims.Mb);     estará guardat a les llistes
//	delete (optims.espai);  estará guardat a les llistes

//	free(xo.act);  s'esborren al sortir de calcular_corba_en_un_sentit() i calcular_corba_en_sentit_contrari()
//	free(xo.ant);  s'esborren al sortir de calcular_corba_en_un_sentit() i calcular_corba_en_sentit_contrari()

	/* borrar les llistes de subespais */

	/* borrar les llistes de sortida */

    delete ll_pop;
}

float espai::obtenir_VTG(float **xm){
 float l;
 char c; //###

 if ((profundidad == PROF_REQ)||(Dim==1)||(ll_pt->n_punts()< NPTMIN*Dim)){  //falta verificacio de que son prou punts tenint en compte la dim
   //printf("calcula stv prof %i, np %i \n",profundidad,ll_pt->n_punts()); //###
   
	// Inicialización para obtener el STV
	ll_pt->inicialitzacio_final(); 
    calcular_htail_delta_xomig_epsx();	 

	GTV = obtenir_STV();   
   //printf(" %f \n",GTV); //###
   delete ll_pop;ll_pop =NULL;
 } 
 else {

   ll_pop = new ll_pnt();

   // Inicialización para obtener el GTV
   ll_pt->inicialitzacio_principal();
   calcular_htail_delta_xomig_epsx();	

   l = calcular_corba_en_un_sentit(); // al fer cadascuna de les corbes dels subespais borrem el llpt que li hem pasat.
     /* girem  */
   l -= calcular_corba_en_sentit_contrari();  // pot no donar cap pop. 

   GTV = finalitzacio();  /* d'on obtenim la VTG de la corba */
 }
 delete ll_pt;	// borrarem els punts obtingut el VTG. Tenim els ppp no volem la ll_pt per res.
// free(xomig);  es el punt que li pasem a l'espai superior, es borrarà desde allà.
 free(eps_x);
 *xm = xomig-1;   // xmean que li pasa el ll_p o pop mig de la corba.
 return GTV;

}

espai  *espai::obtenir_cluster(M_b *Mb,m_d_s *mds){
  int   validacio = 1;
  float *v_pnt;
  float *n_pnt = NULL;
  ll_p  *n_ll_pt;
  float  sum_w,w;
  char c;//###
  int i; //file output

  v_pnt =ll_pt->primer_candidat_clt();  
  n_pnt =Mb->aplicar(v_pnt);   // no cal comprobar la distancia al hpla, ja que dmax smpre serà < htail(i els desplaçaments de l'xmean sempre serán menors que htail.
  
  /* calcular densitats del cluster */

  w = kernel(n_pnt[X]/h_tail)*(*(v_pnt-1));// sempre serà abs(n_pnt[X]/h_tail)<1
  if (w<=0) {//###
      w = 0;	  
  } //###
  sum_w = w;   		  				   // necesari per calcular smoth_mean i density
  n_pnt[X] = w;                        // en la coord de la altura guardamos el peso.
  /* crear e insertar primer punt al cluster */

  n_pnt=treure_coord(n_pnt);
  n_ll_pt = new ll_p(Dim-1);
  n_ll_pt->add_ordX_principal(n_pnt);


  v_pnt = ll_pt->seguent_candidat_clt(validacio);
  while (v_pnt &&  fi_corba(v_pnt)){ // se validará como !fi_corba(v_pnt) aunque v_pnt no pertenezca al cluster
	    n_pnt =Mb->aplicar(v_pnt);
		if (validacio=dist_al_pla(n_pnt)){

			  /* calcular densitat del cluster */

			  w =kernel(n_pnt[X]/h_tail)*(*(v_pnt-1));
			  sum_w += w;						// necesari per calcular smoth_mean i density
              n_pnt[X] = w;                        // en la coord de la altura guardamos el peso.

			  /* afegir punt al cluster */

			  n_pnt = treure_coord(n_pnt);
			  n_ll_pt->add_ordX_principal(n_pnt);
		}
		v_pnt =ll_pt->seguent_candidat_clt(validacio);
  }
  if (v_pnt){
	 while (v_pnt){
		 n_pnt =Mb->aplicar(v_pnt);
		 if (validacio=dist_al_pla(n_pnt)){

		  /* calcular densitat del cluster */
			  w =kernel(n_pnt[X]/h_tail)*(*(v_pnt-1));
			  sum_w += w;						// necesari per calcular smoth_mean i density
			  n_pnt[X] = w;                   // en la coord de la altura guardamos el peso.

		 /* insertar nou punt al cluster */
			  n_pnt = treure_coord(n_pnt);
			  n_ll_pt->add_ordX_principal(n_pnt);
		 }
		 v_pnt =ll_pt->seguent_candidat_clt(validacio); // cal demostracio punt fixe
	 }

		 /* calcul final de densitats */
	 mds->span    =  (float)n_ll_pt->n_punts()/(float)ll_pt->n_punts();
	 mds->density =  sum_w/(ll_pt->n_punts()*h_tail);
	 
	 return (new espai(n_ll_pt,Dim-1,profundidad+1));
  }
  else{
		  delete n_ll_pt;
		  return NULL;   // tots els candidats han estat previament insertats en altres clusters, no queden nous punts per tractar.

  }
}

int espai::dist_al_pla(float *n_punt){   // comproba si el punt es a distancia h_tail del pla
	 return (fabs(n_punt[X]) < h_tail);
}

float *espai::treure_coord(float *n_pnt){  // fem la projeccio sobre el pla per pasar al subespai de dimensio inferior
	 return n_pnt+1;       // el float de memoria que es deja lo borra la ll_p al borrar los puntos, tb se utiliza para obtenir_STV;
}


void espai::rebre_M_a(M_a *n_Ma){ // pasarem el M_a al subespai inferior despres d'obtenir l'espai mes optim
	Ma =n_Ma;
}


int espai::fi_corba(float *v_pnt){
	float *Mba_pnt;
	if (bficorba) {
		Mba_pnt =optims.Mb_ant->aplicar(v_pnt);  
	if (Mba_pnt[X]>2*delta) bficorba = FALSE;	//### h_tail subst. 2*delta 	
//	if (Mba_pnt[X]>h_tail) bficorba = FALSE;		
		else return TRUE;
    }            
    return FALSE;
}


int espai::no_creua_corba(float *ncand){
// els creuaments estarán en funció del delta. Amb un menor delta tindrem més presició a tant l'hora de definir la corba com de finalitzar-la.
    void *pt;
    float *punt;

	ll_pop->resetpt(&pt);
	while(ll_pop->noend(pt)){
		if ((distancia(((pop *)ll_pop->llpt(pt))->alpha,ncand))<delta){ // si la corba es un cercle, l'ultim punt estara a dist delta del xomig.
		   punt = optims.Mb_ant->aplicar(((pop *)ll_pop->llpt(pt))->alpha);
           if (punt[X]>0)  return FALSE;
        }
	    ll_pop->advpt(&pt);
	}
	return TRUE;
}


void espai::calcular_Mb(int ejegir,M_b *Mb,float porcion_pinza){
	int i;       /* el ejegir controla el punto fijo, porcion_pinza no variara en las llamadas recursivas */
	float VTG;
	m_d_s mds;
	espai *espai;
	char c; //###

	if (ejegir){
/*		if (profundidad == 0) //###
		printf("ini ejegir %i \n",ejegir); //###
*/				/* cas recursiu */
		  for(i=-NPARTS/2;i<0;i++){
				calcular_Mb(ejegir-1,Mb->girar(ejegir,i*porcion_pinza),porcion_pinza);
		  }
		  for(i=1;i<=NPARTS/2;i++){
				calcular_Mb(ejegir-1,Mb->girar(ejegir,i*porcion_pinza),porcion_pinza);
		  }          
		  calcular_Mb(ejegir-1,Mb,porcion_pinza);  
/*		if (profundidad == 0){ //###
        printf("fi ejegir %i \n",ejegir); //###
//		c = getchar();  //###
		} //###
*/	}
	else {
			 /* cas directe */
          mds.xmean = NULL;   // per si estem a fi de corba.
		  Mb->calcular_la_inversa();
		  espai = obtenir_cluster(Mb,&mds);  // controla si hemos acabado la curva en el sentido actual.
		  if(espai){
                 espai->rebre_M_a(Mb->donar_M_a(Ma)); // li pasem el Ma per situarse a les coordenades originals aplicat al seu pla(Mb)
				 VTG = espai->obtenir_VTG(&mds.xmean);
				 if(VTG < optims.VTG ){
					delete optims.Mb;   // borramos un Mb_act!= Mb, único y ya tratado.
					delete optims.espai;
					optims.VTG   = VTG;
					optims.Mb    = Mb;
					optims.espai = espai;
					delete optims.mds.xmean;
					optims.mds.xmean = Mb->desaplicar(mds.xmean);// obt posible pop
					optims.mds.span  = mds.span;
					optims.mds.density = mds.density;					

				 }
				 else {
					delete Mb;
					delete espai;
                 }
		  } // aunque este cluster nos indique que ya no hay mas puntos por tratar, hemos de comprovarlo tambien para los otros clusters
		  else delete Mb;
		  free(mds.xmean);  
	}
}



float espai::calcular_corba_en_un_sentit(){
// nomes es diferencia de la funcio calcular_corba_en_sentit_contrari per la
//inicialitzacio i per la crida a la insercio de elements a les llistes de
//sortida, add() vs addrev().  Els vectors de b_ast estarán orientats segons 
//el eigenvector, i las distancies de I serán negatives i positives als dos costats de xomig.

  pop   *n_pop;
  float pinza;
  float *n_bo,*b_opt;
  float *v_xact2xm = NULL;
  float *v_xant2xm = NULL;
  int   naux_delta;
  float sum;
  float lambda;
  float alfa;
//  char c; //###
//  int i; //file output

 /* trobar la 1ra comp. principal dels punts  */

 optims.Mb = new M_b(Dim,obtenir_bo_inicial(&alfa)); //###
// optims.Mb = new M_b(Dim,mult_esc(-1,obtenir_bo_inicial()));
 
  /* calcul htail y delta */

 lambda = suma_d/(pow(ll_pt->n_punts(),((float)(Dim-1))/Dim)*Bmst())*
          pow(((float)(Dim-1))/(1-((1-LD)*alfa+LD)),((float)(Dim-1))/(2*Dim))*
          (1/sqrt(2*PI))* pow(((float)(Dim-1))/Dim,Dim/2.);

 h_tail = C_H * 2.214 * pow(4./3,1./5)*lambda*pow(ll_pt->n_punts(),-1./5); // el 2.214 viene de pasar de kernel normal a Epanechnikov
 delta = C_D*h_tail;     // siempre será menor que el h_tail. De esta forma el nuevo xo ha de estar a distancia dmax de un punto de ll_p.

/* trobar el pop mig de la corba */
 
 optims.Mb_ant = optims.Mb->replicar();  // Mb_ant necesari per calcular la operacio ficorba(). Mb es borrarà validat un ficorba(), amb lo que no usarem més el Mb_ant
 optims.Mb_ant->calcular_la_inversa();
 n_bo = mult_esc(-2*h_tail,optims.Mb_ant->donar_bopt()); // es necesari donar xoant al Mb_ant per calcular la operacio ficorba().
 xo.ant = sum_v(xomig,n_bo);
 optims.Mb_ant->rebre_xo(xo.ant);
 delete n_bo;
 
 optims.espai = NULL;
 xo.act = NULL;
 optims.mds.xmean = (float *) malloc(Dim*sizeof(float));
 memmove(optims.mds.xmean,xomig, Dim*sizeof(float)); // ?el xomig passat el necesita la clase llp?.
 do {
    delete optims.espai;
	delete v_xact2xm;
	delete xo.act;
    xo.act = optims.mds.xmean;
	optims.Mb->rebre_xo(xo.act);// aquesta Mb sera inmediatament eliminat ja que te VTG=INF, pero difondra el xo_act entre la resta de matrius candidates a cada gir.

    /* calcular SVG optima per el xo donat */

    optims.VTG = INF;   // el optims.Mb contdrá una replica del optimo del advance_cluster() anterior que es perdrá a calcular_Mb si no hem acabat la corba en aquest sentit.
    optims.espai = NULL;   // l'anterior espai i l'anterior mb optims no es borren, ja que els apunten desde les llistes resultat.
    optims.mds.xmean = NULL;  // l'xo que hem pasat a Mb no es pot eliminar fins que sigui substituit per un altre. 
    pinza = PI/2;
    while (pinza > PI/16){
	    calcular_Mb(Dim-1,optims.Mb->replicar(),pinza/NPARTS); // pasamos el tamaño de una porcion de la pinza
	    pinza = pinza / NPARTS;       // el tamaño de la nueva pinza pasa a ser el de una porcion de la antigua pinza
	}
    v_xact2xm = dif_v(optims.mds.xmean,xo.act);	 
 }while(major(v_xact2xm, eps_x));
 
 delete xo.act;             // borrem l'anterior candidat. A partir d'aquest moment la copia del xomig estarà ja eliminada.
 xo.act = optims.mds.xmean;

 optims.mds.xmean = allargar(optims.mds.xmean); // el nou xmean optim será el que es guardará a les llistes.
 xomig = optims.mds.xmean;  // nou pop mig de la corba. Aquest no el podem perdre.

 optims.Mb->rebre_xo(optims.mds.xmean); // rep la replica allargada per la Mb optima. Amb aquesta es calculará Ma.  
// optims.espai->rebre_M_a(optims.Mb->donar_M_a(Ma)); // li pasem el Ma per situarse a les coordenades originals aplicat al seu pla(Mb)

 n_pop = (pop *) malloc(sizeof(pop));
 sum = 0;
 n_pop->I = 0;  
 delete xo.ant;     // xo de l'anterior Mb_ant.
 xo.ant = xo.act;     // El xoant será eliminat passat el seguent cluster.

 b_opt = optims.Mb->donar_bopt();   // ens servira per actualitzar el cluster
 delete optims.Mb_ant;
 optims.Mb_ant = optims.Mb;   
 optims.Mb = optims.Mb->replicar(); // utilitzarem una replica del Mb per ser la matriu inicial del cluster seguent. El seu xo corresponent al pop que es guarda a alpha, será substituit pero no borrat.

 n_pop->alpha = optims.mds.xmean;	 // los puntos y vectores de salidan serán de dimension Dim+profundidad.
 n_pop->b_ast = allargar(b_opt);      // los puntos y vectores de salidan seran de dimension Dim+profundidad, b_opt no esta normalitzaca
 n_pop->var_k = optims.VTG;
 n_pop->span  = optims.mds.span;
 n_pop->density = optims.mds.density;
 n_pop->espai  = optims.espai;

 ll_pop->add(n_pop);	
    /* actualitzem cluster */

 n_bo  = mult_esc(delta,b_opt);  
 optims.mds.xmean = sum_v(xo.act,n_bo); // ho guardo a xmean pq aquest s'asigna a xo a l'inici del bucle pel seguent pop.

 xo.act = NULL;   // xo.ant apunta ara a xo.act, per lo tant es borrará quan es borri xo.ant.
 delete n_bo;

 /* trobar els seg. pops en el sentit original del bo del cluster */

 xo.act = NULL;
 naux_delta = 1;
 while(no_creua_corba(xo.ant)){
			  /* cerquem el pop */
	 bficorba = TRUE;     // sera cert fins que no es demostri que hi ha algún nou punt per tractar.
	 delete v_xact2xm;
	 delete v_xant2xm;
	 delete xo.act;   	 // asignacio necesaria si l'xmean i el xo han resultat molt allunyats.                 
	 xo.act = optims.mds.xmean;  // asignacio necesaria si l'xmean i el xo han resultat molt allunyats o si l'xmean s'ha situat darrera l'iperpla del pop anterior. Si no es torna a substituir será el valor que guardarem a la llista alpha.
	 optims.Mb->rebre_xo(xo.act);// aquesta Mb sera inmediatament eliminat ja que te VTG=INF, pero difondra el xo_act entre la resta de matrius candidates a cada gir.
     ll_pt->trobar_primer_candidat_clt(optims.mds.xmean);// cerca el primer candidat al cluster.

	 /* calcular SVG optima per el xo donat */

	 optims.VTG = INF;   // el optims.Mb contdrá una replica del optimo del advance_cluster() anterior que es perdrá a calcular_Mb si no hem acabat la corba en aquest sentit.
	 optims.espai = NULL;   // l'anterior espai i l'anterior mb optims no es borren, ja que els apunten desde les llistes resultat.
     optims.mds.xmean = NULL;  // l'xo que hem pasat a Mb no es pot eliminar fins que sigui substituit per un altre. 
	 pinza = PINZA_MAX;
	 while (pinza > PINZA_MIN){
		  calcular_Mb(Dim-1,optims.Mb->replicar(),pinza/NPARTS); // pasamos el tamaño de una porcion de la pinza
		  if (optims.VTG == INF) { // per tots els hiperplans posibles, els punts del cluster ja habian sigut tractats en anteriors clusters, ja hem tractat els 2 sentits de la corba. Surtiriem sempre a la primera pinza
				  delete optims.Mb_ant;
			//	  free(xo.act); se eliminara en calcular_corba_en_sentit_contreri()
				  free(xo.ant);
				  return (sum);   // exit quan no hi han mes punts a tractar
		  }
		  pinza = pinza / NPARTS;       // el tamaño de la nueva pinza pasa a ser el de una porcion de la antigua pinza
	 }
// if (profundidad == 0){
//	 			c = getchar(); //###
// }
     v_xact2xm = dif_v(optims.mds.xmean,xo.act);	 
	 v_xant2xm = dif_v(optims.mds.xmean,xo.ant);
	 if ((mult_v(optims.Mb->donar_bopt(),optims.Mb_ant->donar_bopt())<0) ||
		 (mult_v(v_xant2xm,optims.Mb_ant->donar_bopt())<0)){   // si l'angle es <0 estem cambiant el sentit de la corba(teorema del cosinus). Hem d'evitar-ho
		  /* avancem el xo original una mica mes que al avenç inicial */
		  naux_delta++;
		  delete optims.mds.xmean;
		  n_bo = mult_esc(naux_delta*delta,optims.Mb_ant->donar_bopt());
		  optims.mds.xmean = sum_v(xo.ant,n_bo);  // Per el 1er pop avançariem el centre de la corba. En els seguents pop tindriem un futur xo mes allunyat.
		  delete n_bo;
		  delete optims.Mb;
		  optims.Mb = optims.Mb_ant->replicar(); // no s'asigna l'original pq el boptant esta a b_ast. Será necesari asignar-li el nou xo.		  
	 }
	 else
		 if (major(v_xact2xm, eps_x)){ // comprovem que el xo i el xmean no resultin massa allunyats, si hem tingut que desplaçar el xmean segur que es cumplira la condicio per l'accio del naux_delta. eliminem v_xant2xm
		/* els valors no son valids, hem de buscar un altre cluster partint del xmean */
			 delete optims.espai;       
		}else{
		/* els valors son bons */
		/* actualitzem llistes */
        delete xo.act;
        xo.act = optims.mds.xmean;
		optims.mds.xmean = allargar(optims.mds.xmean); // el nou xmean optim será el que es guardará a les llistes.

        optims.Mb->rebre_xo(optims.mds.xmean); // rep la replica allargada per la Mb optima. Amb aquesta es calculará Ma.  
//		optims.espai->rebre_M_a(optims.Mb->donar_M_a(Ma)); // li pasem el Ma per situarse a les coordenades originals aplicat al seu pla(Mb)

		n_pop = (pop *) malloc(sizeof(pop));
		sum += distancia(xo.act,xo.ant); 
		n_pop->I = sum; 
		delete xo.ant;
		xo.ant = xo.act;      // del xmean no fa falta fer copia com amb Mb ja que no es tornem a utilitzar. El xoant será eliminat passat el seguent cluster.

		b_opt = optims.Mb->donar_bopt();   // ens servira per actualitzar el cluster
		delete optims.Mb_ant;              // es borrara tb el b_opt del pop anterior.
		optims.Mb_ant = optims.Mb;  
        optims.Mb = optims.Mb->replicar(); // utilitzarem una replica del Mb per ser la matriu inicial del cluster seguent. El seu xo corresponent al pop que es guarda a alpha, será substituit pero no borrat.

		n_pop->alpha = optims.mds.xmean;	 // los puntos y vectores de salidan serán de dimension Dim+profundidad.
		n_pop->b_ast = allargar(b_opt);      // los puntos y vectores de salidan seran de dimension Dim+profundidad, b_opt no esta normalitzaca
		n_pop->var_k = optims.VTG;
		n_pop->span  = optims.mds.span;
		n_pop->density = optims.mds.density;
		n_pop->espai = optims.espai;

		ll_pop->add(n_pop);

		/* actualitzem cluster */

		n_bo  = mult_esc(delta,b_opt);   
//		delete optims.mds.xmean;  // apunta al pop guardat a alpha. No es pot borrar.
		optims.mds.xmean = sum_v(xo.act,n_bo); // ho guardo a xmean pq aquest s'asigna a xo a l'inici del bucle pel seguent pop.

		xo.act = NULL;   // xo.ant apunta ara a xo.act, per lo tant es borrará quan es borri xo.ant.
//      delete b_opt;   s'eliminara juntament a la Mb a la que pertany. Es pot eliminar pq el que conté b_ast es una copia allargada.
		delete n_bo;
		naux_delta = 1;
	}
 }
 delete optims.Mb_ant; 
//free(xo.act);   s'eliminara en calcular_corba_en_sentit_contreri()
 free(xo.ant);
 return (sum);   // exit per quan hi ha un creuament amb la propia corba

}


float espai::calcular_corba_en_sentit_contrari(){
// nomes es diferencia de l'anterior funcio per la inicialitzacio i per la crida
// a la insercio de elements a les llistes de sortida, addrev(). Els vectors de b_ast.
// estarán orientats segons el eigenvector, i las distancies de I serán negatives.
  pop   *n_pop;
  float  pinza;
  float *n_bo,*b_opt;
  float *v_xact2xm = NULL;
  float *v_xant2xm = NULL;
  int    naux_delta;
  float  sum = 0;
  void  *pt;
  int i;
  char c; //###  

 ll_pt->tornar_a_xomig();

 ll_pop->resetpt(&pt);

 optims.Mb = new M_b(Dim,mult_esc(-1,((pop *)ll_pop->llpt(pt))->b_ast)); // fara el paper de matriu del ppp del pas anterior
 optims.Mb_ant = optims.Mb->replicar(); 
 optims.Mb_ant->calcular_la_inversa();
 xo.ant = (float *)malloc(Dim*sizeof(float));  
 memmove(xo.ant,xomig,Dim*sizeof(float)); // no podem borrar un element contingut a alpha, pero el necesitem el valor per calcular la distancia al nou pop.
 optims.Mb_ant->rebre_xo(xo.ant);

 /* avançar cluster */
 n_bo  = mult_esc(delta,optims.Mb->donar_bopt()); // avanço en sentit contrari. Agafo un cadidat per al punt seguent al pop del xomig, perque aquest pop no el puc ja modificar.
 optims.mds.xmean = sum_v(xo.ant,n_bo); // ho guardo a xmean pq aquest s'asigna a xo a l'inici del bucle pel seguent pop. L'anterior valor de optims.mds.xmean es NULL.
        
 // xo.act = NULL            el primer xo.act que borrarem será l'últim de la corba en sentit contrari. 
 delete n_bo;
 naux_delta = 1;

 while(no_creua_corba(xo.ant)){
						 /* cerquem el pop */
	 bficorba = TRUE;     // sera cert fins que no es demostri que hi ha algún nou punt per tractar.
	 delete v_xact2xm;
	 delete v_xant2xm;
	 delete xo.act;   	 // asignacio necesaria si l'xmean i el xo han resultat molt allunyats.                 
	 xo.act = optims.mds.xmean;  // asignacio necesaria si l'xmean i el xo han resultat molt allunyats o si l'xmean s'ha situat darrera l'iperpla del pop anterior. Si no es torna a substituir será el valor que guardarem a la llista alpha.
	 optims.Mb->rebre_xo(xo.act);// aquesta Mb sera inmediatament eliminat ja que te VTG=INF, pero difondra el xo_act entre la resta de matrius candidates a cada gir.
     ll_pt->trobar_primer_candidat_clt(optims.mds.xmean);// cerca el primer candidat al cluster.

	 /* calcular SVG optima per el xo donat */

	 optims.VTG = INF;   // el optims.Mb contdrá una replica del optimo del advance_cluster() anterior que es perdrá a calcular_Mb si no hem acabat la corba en aquest sentit.
	 optims.espai = NULL;   // l'anterior espai i l'anterior mb optims no es borren, ja que els apunten desde les llistes resultat.
     optims.mds.xmean = NULL;  // l'xo que hem pasat a Mb no es pot eliminar fins que sigui substituit per un altre. 
	 pinza = PINZA_MAX;
	 while (pinza > PINZA_MIN){
		calcular_Mb(Dim-1,optims.Mb->replicar(),pinza/NPARTS); // pasamos el tamaño de una porcion de la pinza
		if (optims.VTG == INF){ // per tots els hiperplans posibles els punts del cluster ja habian sigut tractats en anteriors clusters hem de fer la corba en sentit contrari. Surtiriem sempre a la primera pinza
			  delete optims.Mb_ant; 
			  free(xo.act);
			  free(xo.ant);
			  return (sum);      // exit per quan no hi han mes punts a tractar
		}
		pinza = pinza / NPARTS;       // el tamaño de la nueva pinza pasa a ser el de una porcion de la antigua pinza
	 }
	 
	 v_xact2xm = dif_v(optims.mds.xmean,xo.act);
	 v_xant2xm = dif_v(optims.mds.xmean,xo.ant);
	 if (mult_v(optims.Mb->donar_bopt(),optims.Mb_ant->donar_bopt())<0 ||
	     mult_v(v_xant2xm,optims.Mb_ant->donar_bopt())<0){   // si l'angle es <0 estem cambiant el sentit de la corba. Hem d'evitar-ho
		/* avancem el xo original una mica mes que l'anterior vegada*/
		naux_delta++;
		delete optims.mds.xmean;
		n_bo = mult_esc(naux_delta*delta,optims.Mb_ant->donar_bopt());
		optims.mds.xmean = sum_v(xo.ant,n_bo);  
		delete n_bo;
		delete optims.Mb;
	    optims.Mb = optims.Mb_ant->replicar(); // no s'asigna l'original pq el boptant esta a b_ast. Será necesari asignar-li el nou xo.
	 }
	 else
	 if (major(v_xact2xm, eps_x)){    // comprovem que el xo i el xmean no resultin massa allunyats, si hem tingut que desplaçar el xmean segur que es cumplira la condicio per l'accio del naux_delta. eliminem v_xant2xm
		/* els valors no son valids, hem de buscar un altre cluster partint del xmean */
		 delete optims.espai;  
	}else{
		/* els valors son bons */
		/* actualitzem llistes */
        delete xo.act;
        xo.act = optims.mds.xmean;
		optims.mds.xmean = allargar(optims.mds.xmean); // el nou xmean optim será el que es guardará a les llistes.

        optims.Mb->rebre_xo(optims.mds.xmean); // rep la replica allargada per la Mb optima. Amb aquesta es calculará Ma.  
//		optims.espai->rebre_M_a(optims.Mb->donar_M_a(Ma)); // li pasem el Ma per situarse a les coordenades originals aplicat al seu pla(Mb)

		b_opt = optims.Mb->donar_bopt();   // ens servira per actualitzar el cluster
		n_bo  = mult_esc(-1,b_opt);        // tots els vectors de b_ast tindran la mateixa orientacio.
		delete optims.Mb_ant;              // es  borra tb el b_opt de pop anterior.  
		optims.Mb_ant = optims.Mb;  
        optims.Mb = optims.Mb->replicar(); // utilitzarem una replica del Mb per ser la matriu inicial del cluster seguent. El seu xo corresponent al pop que es guarda a alpha, será substituit pero no borrat.
 
		n_pop = (pop *) malloc(sizeof(pop));
		sum -= distancia(xo.act,xo.ant); // distancies negatives.
		n_pop->I = sum;
		delete xo.ant;
		xo.ant = xo.act;      // del xmean no fa falta fer copia com amb Mb ja que no es tornem a utilitzar. El xoant será eliminat passat el seguent cluster.
		
		n_pop->alpha = optims.mds.xmean;	 // los puntos y vectores de salidan serán de dimension Dim+profundidad.
		n_pop->b_ast = allargar(n_bo);      // los puntos y vectores de salidan seran de dimension Dim+profundidad, b_opt no esta normalitzada. se guarda en sentido contrario
        delete n_bo;
		n_pop->var_k = optims.VTG;
		n_pop->span  = optims.mds.span;
		n_pop->density = optims.mds.density;
		n_pop->espai = optims.espai;

		ll_pop->addrev(n_pop);


		/* actualitzem cluster */
        
		n_bo  = mult_esc(delta,b_opt);  
//		delete optims.mds.xmean;  // apunta al pop guardat a alpha. No es pot borrar.
		optims.mds.xmean = sum_v(xo.act,n_bo); // ho guardo a xmean pq aquest s'asigna a xo a l'inici del bucle pel seguent pop.

		xo.act = NULL;   // xo.ant apunta ara a xo.act, per lo tant es borrará quan es borri xo.ant.
//      delete b_opt;   s'eliminara juntament a la Mb a la que pertany. Es pot eliminar pq el que conté b_ast es una copia allargada del vector amb el sentit corregit.
		delete n_bo;
		naux_delta = 1;

	 }
 }
 delete optims.Mb_ant; 
 free(xo.act);
 free(xo.ant);
 return (sum);   // exit per quan hi ha un creuament amb la propia corba
}

float *espai::allargar(float *b_opt){
	 int i;
	 float *n_b;
	 n_b = (float *) malloc((Dim+profundidad)*sizeof(float));

	 for(i=0;i<profundidad;i++)
			n_b[i]=0.0;

     memmove(n_b+profundidad,b_opt,Dim*sizeof(float)); //### sense testar 

//	  for(i=profundidad;i<Dim+profundidad;i++)
//	 		n_b[i]=b_opt[i-profundidad];
	  
	 return n_b+profundidad;
}

float *espai::obtenir_bo_inicial(float *alfa){

 float **S = (float **) malloc(Dim*sizeof(float *));
 float *A = (float *) malloc(Dim*(Dim+1)/2*sizeof(float));
 float *EV = (float *) calloc(Dim*Dim,sizeof(float));
 float *E = (float *) malloc(Dim*sizeof(float)); 
 float *b_op = (float *) malloc(Dim*sizeof(float));
 void *pt;
 float *punt;
 float *m = (float *) calloc(Dim,sizeof(float));
 int np = ll_pt->n_punts();
 int i,j;
 float VarTot;

 for (i=0;i< Dim;i++) S[i] = (float *) calloc(Dim,sizeof(float)); 
 
 for(i=0;i< Dim;i++){
   /* calculem l'element de la diagonal */
   ll_pt->resetpt(&pt); 
   while(ll_pt->noend(pt)){
	   punt = ((float *)ll_pt->llpt(pt));
	   m[i] += punt[i];
	   S[i][i]+= pow(punt[i],2); // xomig será la mitja de tots els punts per cada  coordenada.
	   ll_pt->advpt(&pt);
   }
   m[i]=m[i]/np;
   S[i][i] = S[i][i]/np - m[i]*m[i];
   for(j=0;j< i;j++){
	   ll_pt->resetpt(&pt); 
	   while(ll_pt->noend(pt)){
		   punt = (float *)ll_pt->llpt(pt);
		   S[i][j]+= punt[i] * punt[j]; // xomig será la mitja de tots els punts per cada  coordenada.
		   ll_pt->advpt(&pt);
	   }
	   S[i][j] = S[i][j]/np - m[i]*m[j];
	   S[j][i] = S[i][j];
   }
 }

 for(i=0;i< Dim;i++){
   for(j=0;j< Dim;j++){
	   A[((i*i+i)/2)+j] = S[i][j]; /* el orden ij será correcto pues tratamos punto por variable */
   }
 }

 eigens(A,EV,E,Dim);

 i = 0;VarTot = 0;
 for(j=0;j< Dim;j++) {
	 VarTot += E[j];
	 if (E[i] < E[j]) i = j; // buscamos el autovector con mayor autovalor
 } 
 for(j=0;j< Dim;j++) b_op[j] = EV[Dim*i+j];  /* EV[i][j] */

 *alfa = E[i]/VarTot;



 return b_op;	  
}

float espai::Bmst(){
   float cd,lbeta,ubeta,S,fact;
   double par_f,res;
   int i;
   

   par_f =(Dim/2.)+1;
   cd = pow(PI,Dim/2.)/exp(gammln(par_f));
   par_f = 1./Dim;
   lbeta = exp(gammln(par_f))/(Dim*pow(cd,1./Dim));
   ubeta = pow(2,1./Dim)*lbeta;

   S=0;
   fact=1.;
   for(i=1;i<=30;i++) {
	   fact=fact*i;
	   par_f = i+(1./Dim)-1;
	   S+= exp(gammln(par_f))/(fact*pow(i,(1./Dim)+1));
   }
   
   res=S/(Dim*pow(cd,1./Dim));
   return res;
}

float espai::gammln(float xx){
	 double x,y,tmp,ser;
	 static double cof[6]={76.18009172947146,-86.50532032941677,
		                   24.01409824083091,-1.231739572450155,
						   0.1208650973866179e-2,-0.5395239384953e-5};
	 int j;
	 y=x=xx;
	 tmp=x+5.5;
	 tmp -= (x+0.5)*log(tmp);
	 ser = 1.00000000190015;
	 for (j=0;j<=5;j++) ser += cof[j]/++y;
	 return -tmp+log(2.5066282746310005*ser/x);
}



float espai::obtenir_STV(){
// en caso de calcular la stv por no poder generar la curva,
// el xomig es el punto central de la curva.  
	
  int j;
  void *pt;
  float *punt;
  float *c_xomig;
  float sum_w = 0.0;
  float stv=0.0;  

  c_xomig = (float *) malloc((Dim+1)*sizeof(float));
  c_xomig[0] = 0;c_xomig++;     // necesari si es un espai final i l'hem de pasar al espai superior.
  memmove(c_xomig,xomig, Dim*sizeof(float));
  free(xomig);xomig = c_xomig;

  ll_pt->resetpt(&pt);
  while (ll_pt->noend(pt)){
	 punt = ll_pt->llpt(pt);
	 // w = kernel(*(punt-1)/h_tail);    *(punt-1) contendrá ja el pes calculat.
	 sum_w += *(punt-1);
	 for (j=0;j<Dim;j++)
		 stv += pow(punt[j]-xomig[j],2)*(*(punt-1));  // *(punt-1) contindra el pes (distancia al pla ponderat per h_tail i kernel)
	 ll_pt->advpt(&pt);
  }
  return stv/sum_w;
}


void espai::calcular_htail_delta_xomig_epsx(){
	float *range, *min,*max;
	float  sum,length;
	int i;
	 
	 ll_pt->donar_max_min_xomig(&max,&min,&xomig,&suma_d);

	 range = dif_v(max,min);
	 eps_x = mult_esc(C_EPS,range);

	 free (range);
	 free (max);
	 free (min);
}

float espai::kernel(float d){   
      return (0.75*(1- pow(d,2)));
}

int espai::major(float *v1,float *v2){
  int i=0;

  while (i<Dim && fabs(v1[i])<=fabs(v2[i])) i++;  //major o igual
  return (i != Dim);
}

float espai::finalitzacio(){

  ll_flt   w_s;	
  void *pt,*ptant;
  void *ptws;
  float *xomig_0,*xomig_1;
  float Iant,I2ant,dant,wsact;
  float sum = 0;  

  /*  w_s  ::  ws(n) = (I(n+1)-I(n-1))*density(n) */

  ll_pop->resetpt(&pt);
  Iant = ((pop *)ll_pop->llpt(pt))->I;
  dant = ((pop *)ll_pop->llpt(pt))->density;
  ll_pop->advpt(&pt);
  wsact = (((pop *)ll_pop->llpt(pt))->I-Iant)*2*dant;
  sum += wsact;
  w_s.add(wsact);
  while(ll_pop->noend(pt)){
     I2ant = Iant;
     Iant = ((pop *)ll_pop->llpt(pt))->I;
     dant = ((pop *)ll_pop->llpt(pt))->density;
	 ll_pop->advpt(&pt);
     wsact = (((pop *)ll_pop->llpt(pt))->I-I2ant)*dant;
	 sum += wsact;
	 w_s.add(wsact);
  }
  wsact = (((pop *)ll_pop->llpt(pt))->I-Iant)*2*((pop *)ll_pop->llpt(pt))->density;
  sum += wsact;
  w_s.add(wsact);

  /* w_s :: w_s(n)/sum(w_s) */

  w_s.resetpt(&ptws);
  while(w_s.noend(ptws)){
	  w_s.modpt(ptws,w_s.llpt(ptws)/sum);
	  w_s.advpt(&ptws);
  }

  /* density  :: density(n) = 2*w_s(n)/(I(n+1)-I(n-1)) */ 

  sum = 0;
  ll_pop->resetpt(&pt);
  w_s.resetpt(&ptws);
  Iant = ((pop *)ll_pop->llpt(pt))->I;
  sum += Iant * w_s.llpt(ptws);
  ll_pop->advpt(&pt);  
  ((pop *)ll_pop->llpt(pt))->density = (w_s.llpt(ptws))/(((pop *)ll_pop->llpt(pt))->I-Iant); 
  while(ll_pop->noend(pt)){
     I2ant = Iant;
     Iant = ((pop *)ll_pop->llpt(pt))->I;
	 w_s.advpt(&ptws);
	 ll_pop->advpt(&pt);
	 sum += Iant * w_s.llpt(ptws);	 
	 ((pop *)ll_pop->llpt(pt))->density = (2*w_s.llpt(ptws))/(((pop *)ll_pop->llpt(pt))->I-Iant);
  }
  w_s.advpt(&ptws);
  sum += ((pop *)ll_pop->llpt(pt))->I * w_s.llpt(ptws);
  ((pop *)ll_pop->llpt(pt))->density = (w_s.llpt(ptws))/(((pop *)ll_pop->llpt(pt))->I-Iant);


  /*   I  :: I(n) = I(n) - sum(I*ws) */
  
  /*   Var_PB :: curva principal */
  /*   Var_res :: variancia residual */

  Var_PC =0.0;
  Var_res =0.0;
  w_s.resetpt(&ptws);
  ll_pop->resetpt(&pt);
  ((pop *)ll_pop->llpt(pt))->I =((pop *)ll_pop->llpt(pt))->I -sum;
  while (((pop *)ll_pop->llpt(pt))->I < 0){  // sempre hi haura un canvi de positiu a negati
	  Var_PC += pow(((pop *)ll_pop->llpt(pt))->I,2)*w_s.llpt(ptws);
      Var_res += ((pop *)ll_pop->llpt(pt))->var_k  *w_s.llpt(ptws);
      ptant = pt;
	  ll_pop->advpt(&pt);w_s.advpt(&ptws);
	  ((pop *)ll_pop->llpt(pt))->I =((pop *)ll_pop->llpt(pt))->I -sum;
  }
   /*asignacio del pop centre de la corba que es pasa a l'espai superior */

//free (xomig) ; no el borrem, conte un pop de ll_pop->
  xomig = (float *)malloc((Dim+1)*sizeof(float));
  xomig[0] =0;xomig++;
  if (!((pop *)ll_pop->llpt(pt))->I) {
     xomig = mult_esc(((pop *)ll_pop->llpt(pt))->I,
	                  ((pop *)ll_pop->llpt(ptant))->alpha);
     xomig_0 = mult_esc(((pop *)ll_pop->llpt(ptant))->I,
	                    ((pop *)ll_pop->llpt(pt))->alpha);
     xomig_1 = sum_v(xomig,xomig_0);
     delete xomig_0;
     delete xomig;
     xomig = mult_esc(((pop *)ll_pop->llpt(pt))->I
	                  *((pop *)ll_pop->llpt(ptant))->I
		  		      ,xomig_1);     // xomig que es pasa al espai superior.
     delete xomig_1;
  }    
  else memmove(xomig,((pop *)ll_pop->llpt(pt))->alpha,Dim*sizeof(float));

  ll_pop->advpt(&pt);w_s.advpt(&ptws);
  while (w_s.noend(ptws)){
	  ((pop *)ll_pop->llpt(pt))->I =((pop *)ll_pop->llpt(pt))->I -sum;
	  Var_PC += pow(((pop *)ll_pop->llpt(pt))->I,2)*w_s.llpt(ptws);
	  Var_res += ((pop *)ll_pop->llpt(pt))->var_k  *w_s.llpt(ptws);
      ll_pop->advpt(&pt);w_s.advpt(&ptws);
  }


  /* GTV:: gloval total variance */

  printf("var_PC: %f ;  var_res: %f \n", Var_PC, Var_res);
  return( Var_PC + Var_res);  // GTV

}

void  espai::obtenir_data(FILE *fd_s){

  int i;
  ll_pnt *sll_pop;
  espai *sespai;
  float *auxa;
  float *auxb;
  void *pt,*pt2;

//###   file output
//printf("llegamos al guardar en file \n");
  ll_pop->resetpt(&pt);
  while(ll_pop->noend(pt)){
	 fprintf(fd_s, "0 %f %f %f %f ",((pop *)ll_pop->llpt(pt))->I,((pop *)ll_pop->llpt(pt))->density,
		     ((pop *)ll_pop->llpt(pt))->span,((pop *)ll_pop->llpt(pt))->var_k);

	 
	 for (i =0;i<Dim;i++) 
        fprintf(fd_s,"%f ",((pop *)ll_pop->llpt(pt))->alpha[i]);
	 for (i =0;i<Dim;i++) 
        fprintf(fd_s,"%f ",((pop *)ll_pop->llpt(pt))->b_ast[i]);
     fprintf(fd_s,"\n");


     sespai = ((pop *)ll_pop->llpt(pt))->espai;
	 sll_pop = sespai->ll_pop;

     // subspai // 
	 if (PROF_REQ==2 && sll_pop){
	  sll_pop->resetpt(&pt2);
		while(sll_pop->noend(pt2)){
		 fprintf(fd_s, "1 %f %f %f %f ",((pop *)sll_pop->llpt(pt2))->I,((pop *)sll_pop->llpt(pt2))->density,
		     ((pop *)sll_pop->llpt(pt2))->span,((pop *)sll_pop->llpt(pt2))->var_k);


		 auxb = sespai->Ma->aplicar_Ma_vect(((pop *)sll_pop->llpt(pt2))->b_ast);
         auxa = sespai->Ma->aplicar_Ma_punt(((pop *)sll_pop->llpt(pt2))->alpha); 
		 for (i =0;i<Dim;i++) 
	        fprintf(fd_s,"%f ",auxa[i]);
		 for (i =0;i<Dim;i++) 
	        fprintf(fd_s,"%f ",auxb[i]);
	     fprintf(fd_s,"\n");

		 delete auxb; delete auxa;

		 sll_pop->advpt(&pt2);

	  }
	 fprintf(fd_s, "1 %f %f %f %f ",((pop *)sll_pop->llpt(pt2))->I,((pop *)sll_pop->llpt(pt2))->density ,
		     ((pop *)sll_pop->llpt(pt2))->span,((pop *)sll_pop->llpt(pt2))->var_k);

	 auxb = sespai->Ma->aplicar_Ma_vect(((pop *)sll_pop->llpt(pt2))->b_ast);
     auxa = sespai->Ma->aplicar_Ma_punt(((pop *)sll_pop->llpt(pt2))->alpha); 
	 for (i =0;i<Dim;i++) 
	    fprintf(fd_s,"%f ",auxa[i]);
	 for (i =0;i<Dim;i++) 
	    fprintf(fd_s,"%f ",auxb[i]);

	 fprintf(fd_s,"\n");
	 }
     // end subspai

	 ll_pop->advpt(&pt);
  }
 fprintf(fd_s, "0 %f %f %f %f ",((pop *)ll_pop->llpt(pt))->I,((pop *)ll_pop->llpt(pt))->density ,
		     ((pop *)ll_pop->llpt(pt))->span,((pop *)ll_pop->llpt(pt))->var_k);

 for (i =0;i<Dim;i++) 
        fprintf(fd_s,"%f ",((pop *)ll_pop->llpt(pt))->alpha[i]);
 for (i =0;i<Dim;i++) 
        fprintf(fd_s,"%f ",((pop *)ll_pop->llpt(pt))->b_ast[i]);
 fprintf(fd_s,"\n");

 fclose(fd_s);

}
//### end file output



// modificado 16/4/2002 inicializamos las variables
int espai :: PROF_REQ;
int espai :: NPARTS;
float espai :: C_H;
float espai :: C_D;

void espai::inicializar_nparts_ch_cd(int profreq,int nparts,float c_h,float c_d)
{
	if(profreq<=0) PROF_REQ=1;
	else PROF_REQ=profreq;
	if(nparts<3 || nparts>6) NPARTS=4;
	else NPARTS=nparts;
	if(c_h<0.5 || c_h>1.5) C_H = 0.75;
	else C_H=c_h;
	if(c_d<0.25 || c_d>0.5) C_D = 0.3;		// siempre menor que 0.5
	else C_D=c_d;
}


// fin




/*	Eigenvalues and eigenvectors of a real symmetric matrix
 *
 * SYNOPSIS:
 *
 * int n;
 * double A[n*(n+1)/2], EV[n*n], E[n];
 * void eigens( A, EV, E, n );
 *
 * DESCRIPTION:
 *
 * The algorithm is due to J. vonNeumann.
 *                   -     -
 * A[] is a symmetric matrix stored in lower triangular form.
 * That is, A[ row, column ] = A[ (row*row+row)/2 + column ]
 * or equivalently with row and column interchanged.  The
 * indices row and column run from 0 through n-1.
 *
 * EV[] is the output matrix of eigenvectors stored columnwise.
 * That is, the elements of each eigenvector appear in sequential
 * memory order.  The jth element of the ith eigenvector is
 * EV[ n*i+j ] = EV[i][j].
 *
 * E[] is the output matrix of eigenvalues.  The ith element
 * of E corresponds to the ith eigenvector (the ith row of EV).
 *
 * On output, the matrix A will have been diagonalized and its
 * orginal contents are destroyed.
 *
 * ACCURACY:
 *
 * The error is controlled by an internal parameter called RANGE
 * which is set to 1e-10.  After diagonalization, the
 * off-diagonal elements of A will have been reduced by
 * this factor.
 *
 * ERROR MESSAGES:
 *
 * None.
 *
 */
/*
Copyright 1973, 1991 by Stephen L. Moshier
Copyleft version.
*/

void espai::eigens(float *A,float *RR,float *E,int N )
//double A[], RR[], E[];
//int N;
{
int IND, L, LL, LM, M, MM, MQ, I, J, K, IA, LQ;
int IQ, IM, IL, NLI, NMI;
double ANORM, ANORMX, AIA, THR, ALM, QI, ALL, AMM, Xv, Y; 
double SINX, SINX2, COSX, COSX2, SINCS, AIL, AIM;
double RLI, RMI, Q, V;
//double sqrt(), fabs();
static double RANGE = 1.0e-10; /*3.0517578e-5;*/


/* Initialize identity matrix in RR[] */
//for( J=0; J<N*N; J++ )  // ja està inicialitzada.
//	RR[J] = 0.0;

MM = 0;
for( J=0; J<N; J++ )
	{
	RR[MM + J] = 1.0;
	MM += N;
	}

ANORM=0.0;
for( I=0; I<N; I++ )
	{
	for( J=0; J<N; J++ )
		{
		if( I != J )
			{
			IA = I + (J*J+J)/2;
			AIA = A[IA];
			ANORM += AIA * AIA;
			}
		}
	}
if( ANORM <= 0.0 )
	goto done;
ANORM = sqrt( ANORM + ANORM );
ANORMX = ANORM * RANGE / N;
THR = ANORM;

while( THR > ANORMX )
{
THR=THR/N;

do
{ /* while IND != 0 */
IND = 0;

for( L=0; L<N-1; L++ )
	{

for( M=L+1; M<N; M++ )
	{
	MQ=(M*M+M)/2;
	LM=L+MQ;
	ALM=A[LM];
	if( fabs(ALM) < THR )
		continue;

	IND=1;
	LQ=(L*L+L)/2;
	LL=L+LQ;
	MM=M+MQ;
	ALL=A[LL];
	AMM=A[MM];
	Xv=(ALL-AMM)/2.0;
	Y=-ALM/sqrt(ALM*ALM+Xv*Xv);
	if(Xv < 0.0)
		Y=-Y;
	SINX = Y / sqrt( 2.0 * (1.0 + sqrt( 1.0-Y*Y)) );
	SINX2=SINX*SINX;
	COSX=sqrt(1.0-SINX2);
	COSX2=COSX*COSX;
	SINCS=SINX*COSX;

/*	   ROTATE L AND M COLUMNS */
for( I=0; I<N; I++ )
	{
	IQ=(I*I+I)/2;
	if( (I != M) && (I != L) )  // no ens interesen els eigenvalues, només els eigenvectors
		{
		if(I > M)
			IM=M+IQ;
		else
			IM=I+MQ;
		if(I >= L)
			IL=L+IQ;
		else
			IL=I+LQ;
		AIL=A[IL];
		AIM=A[IM];
		Xv=AIL*COSX-AIM*SINX;
		A[IM]=AIL*SINX+AIM*COSX;
		A[IL]=Xv;
		}
	NLI = N*L + I;
	NMI = N*M + I;
	RLI = RR[ NLI ];
	RMI = RR[ NMI ];
	RR[NLI]=RLI*COSX-RMI*SINX;
	RR[NMI]=RLI*SINX+RMI*COSX;
	}

	Xv=2.0*ALM*SINCS;             // no ens interesen els eigenvalues, només els eigenvectors
	A[LL]=ALL*COSX2+AMM*SINX2-Xv;
	A[MM]=ALL*SINX2+AMM*COSX2+Xv;
	A[LM]=(ALL-AMM)*SINCS+ALM*(COSX2-SINX2);
	} /* for M=L+1 to N-1 */
	} /* for L=0 to N-2 */

	}
while( IND != 0 );

} /* while THR > ANORMX */

done:	;

/* Extract eigenvalues from the reduced matrix */
L=0;            
for( J=1; J<=N; J++ )   // no ens interesen els eigenvalues.
	{
	L=L+J;
	E[J-1]=A[L-1];
	}
}


// ops vect /////////////////

float  espai::distancia(float *pnt1,float *pnt2){
  int i;
  float sum=0.0;
  for(i =0;i<Dim;i++)
		sum += pow(pnt1[i]-pnt2[i],2);
  return sqrt(sum);
}

float *espai::mult_esc(float e,float *v){
 int i;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  v3[i] = v[i]*e;
 return v3;
}

float espai::mult_v(float *v1,float *v2){
 int i;
 float sum = 0.0;
 for(i=0;i<Dim;i++){
		sum += v1[i]* v2[i];
 }
 return sum;
}


float *espai::sum_v (float *v1,float *v2){
 int i;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  v3[i] = v1[i]+v2[i];
 return v3;
}


float *espai::dif_v (float *v2,float *v1){
 int i;
 float *v3;

 v3 = (float *)malloc(Dim* sizeof(float));
 for(i=0;i<Dim;i++)  v3[i] = v2[i]-v1[i];
 return v3;
}

float *espai::norma_v(float *v){
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

