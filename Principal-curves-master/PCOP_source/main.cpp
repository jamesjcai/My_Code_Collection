//#include "Ll_pnt.h"
#include "Fileops.h"
#include "espai.h"

#include "stdio.h"

// modificado 9/4/2002 15/4/2002 16/4/2002

//int NPARTs;
//int PROFREQ;
//float C_h;  
//float C_d;     // siempre menor que 0.5


// fin

fileops inout;
int Dim;




void main (){

espai *psp;
float **Ma;
int i;
ll_p *ll_pt;

float *punt;
int endfile;
float vtg;

float *mx;

// modificado 9/4/2002 12/4/2002  15/4/2002
char *datafile;	//se guardará el nombre del fichero donde se encuentran los datos.
char *outfile;	//se guardará el nombre del fichero donde se encuentran los datos.
datafile = (char *) malloc((26)*sizeof(char)); //le hemos puesto el tamaño de 26 char
outfile = (char *) malloc((26)*sizeof(char)); //le hemos puesto el tamaño de 26 char

int profreq,nparts;
float c_d,c_h;


// inicialización valores por defecto, despues en el fichero de setup puede que cambien
// modificado 16/4/2002
//PROFREQ =1;
//NPARTs = 4;
//C_h = 0.75;
//C_d = 0.4;		// siempre menor que 0.5

//inout.open_read_datafile_setup("test/setup.dat",&datafile);
inout.open_read_datafile_setup("setup.dat",&datafile,&outfile,&profreq,&nparts,&c_d,&c_h);
////fin

inout.open_read_datafile(datafile);

 punt = inout.donar_primer_punt(&Dim,&endfile);

 ll_pt = new ll_p(Dim);   
 
 ll_pt->add_ordX_principal(punt); //###
 
 while(!endfile){
  
 ll_pt->add_ordX_principal(inout.donar_seg_punt(Dim,&endfile)); //###
 }

 
 Ma = (float **) malloc(Dim*sizeof(float *));

 for (i=0;i<Dim;i++){
  Ma[i] = (float *)calloc(Dim,sizeof(float));
  Ma[i][i]=1;
 }

 mx = (float *) calloc(Dim,sizeof(float));
 
 psp = new espai(ll_pt,Dim,0);
 
 //modificacion 16/4/2002
 psp->inicializar_nparts_ch_cd(profreq,nparts,c_h,c_d);
 //fin

 psp->rebre_M_a(new M_a(Dim,0,Ma,mx));
 vtg = psp->obtenir_VTG(&mx);
 
 psp->obtenir_data(inout.open_write_datafile(outfile));

 delete psp;

 printf("vtg: %f \n",vtg);

 //modificacion 15/4/2002
 free(datafile); //libero la memoria ocupada por el nombre del fichero.
 free(outfile); //libero la memoria ocupada por el nombre del fichero.
 //fin 15/4/2002
}
