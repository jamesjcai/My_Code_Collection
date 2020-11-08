#include "stdio.h"
#include <stdlib.h>
#include <memory.h>
#include <string.h>


#define TRUE  1
#define FALSE 0

class fileops{
private:

FILE *fd_s;
FILE *fd_e;
char first_char; // guardamos el primer char valido despues de leer cada float. lo usamos para ver que realmente habra un siguiente float


float leer_float(int *endcoord);

public:

fileops();
~fileops(){};

//modificado 9/4/2002 12/4/2002
// entrada y lectura del setup
void open_read_datafile_setup(char *nom_fitxer,char **datafile,char **outfile,int *profreq,int *nparts,float *c_d,float *c_h);
//fin

void open_read_datafile(char *pth);
FILE *open_write_datafile(char *pth);
// entrada
float *donar_primer_punt(int *d,int *endfile);
float *donar_seg_punt(int Dim,int *endfile);
//sortida
void sortida_str(char *s);
};
