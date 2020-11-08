
#include "fileops.h"
#include <conio.h>



fileops::fileops(){
  fd_e = NULL;
  fd_s = NULL;
 
}


FILE *fileops::open_write_datafile(char *pth){
  if (fd_s) fclose(fd_s);
  if(!(fd_s = fopen(pth,"w"))){
	  printf("error: opening data file: %s \n",pth);
	  exit(-1);
  }
  return fd_s;
}
// modificado 9/4/2002
void fileops::open_read_datafile_setup(char *nom_fitxer,char **datafile,char **outfile,int *profreq,int *nparts,float *c_d,float *c_h)
{
  char command[26];
  int outfileread=FALSE;
  // 16/4/2002 inicializamos variables 
  *profreq=NULL;
  *nparts=NULL;
  *c_d=NULL;
  *c_h=NULL;

  if (fd_e) fclose(fd_e);
  //
  //modificado 12/4/2002 15/4/202
  if ((fd_e=fopen(nom_fitxer,"rt"))==NULL) {
		printf("Error en obrir el fitxer de comandes\n");
		getch();
		exit(-1);
	}
   
   while (fscanf(fd_e,"%s",command)!=EOF)
   {
	// COMANDES GENERALS 
	printf("command: %s \n",command);
//	getch();

	if (strcmp(command,"datafile")==0) 
	{
		fscanf(fd_e,"%s",*datafile);
	}
	else if (strcmp(command,"outfile")==0) 
	{
		fscanf(fd_e,"%s",*outfile);
		outfileread = TRUE;
	}
	else if (strcmp(command,"PROFREQ")==0)
	{
		fscanf(fd_e,"%d",profreq);
	}
	  else if (strcmp(command,"NPARTS")==0) 
		{
			fscanf(fd_e,"%d",nparts);
		}
		else if (strcmp(command,"C_H")==0) 
			{
				fscanf(fd_e,"%f",c_h);
			}
			else if (strcmp(command,"C_D")==0) 
				{
					fscanf(fd_e,"%f",c_d);
	  	  
				}
				// COMANDO NO RECONOCIDO 
				else 
				{
					printf("command: %s no reconeguda\n",command);
					getch();
					exit(-1);	
				}
   }
   if(!outfileread)
   		sprintf(*outfile,"%d.alpha",datafile); 
 
}
//fin

void fileops::open_read_datafile(char *pth){
  char trash[20];
  int i=-1;

  if (fd_e) fclose(fd_e);
  
  if(!(fd_e = fopen(pth,"r"))){
	  printf("error: opening data file: %s \n",pth);
	  exit(-1);
  }

  if feof(fd_e){
	  printf("error: void data file: %s \n",pth);
	  exit(-1);
  }
  do{
	 trash[++i] = fgetc(fd_e);              // saltamos tanto ' ' como '/n' iniciales
  }while ((trash[i]== ' ' || trash[i]=='\n') && !feof(fd_e));

  if feof(fd_e){
	  printf("error: void data file: %s \n",pth);
	  exit(-1);
  }
  else first_char = trash[i];

}

float *fileops::donar_primer_punt(int *d,int *endfile){
  float punt[20];
  float *d_punt;
  int endcoord;
  int i=0;

  punt[i++] = leer_float(&endcoord);
  while ( !endcoord ){
	  punt[i++] = leer_float(&endcoord);
  }
  *d = i;
  *endfile = feof(fd_e);

  d_punt = (float *) malloc((i+1)*sizeof(float));
  d_punt[0] = 1;d_punt++; // coordenada -1 pel pes.
  memmove(d_punt,punt,i*sizeof(float));
  return d_punt;
}

float *fileops::donar_seg_punt(int Dim,int *endfile){
  float *punt;
  int endcoord;
  int i=0;

  punt = (float *)calloc((Dim+1),sizeof(float));
  punt[0]=1;punt++; // coordenada -1 pel pes.
  punt[i++] = leer_float(&endcoord);
  while ( !endcoord ){
	  punt[i++] = leer_float(&endcoord);
  }
  
  *endfile = feof(fd_e);
  return punt;
}

void fileops::sortida_str(char *s){
  fwrite(&s, sizeof(s), 1,fd_s);
}

////private

float fileops::leer_float(int *endcoord){
  char ch[20];
  char trash[20];
  int i=-1;

  *endcoord = FALSE;
  memset(ch,0,20);
  ch[++i]= first_char; // guardado de la lectura anterior
  do{
	  ch[++i] = fgetc(fd_e);
  }while (ch[i]!= ' ' && ch[i]!= '\n' && !feof(fd_e));

  trash[0]=ch[i];
  i=0;
  if(trash[i]=='\n') *endcoord = TRUE; 
  while ((trash[i]== ' ' || trash[i]=='\n') && !feof(fd_e)){
	 trash[++i] = fgetc(fd_e);
	 if(trash[i]=='\n') *endcoord = TRUE; 
  }

  if feof(fd_e){
	      *endcoord = TRUE;   
	      return atof(ch); 
  }
  else{
		  first_char = trash[i];
		  return atof(ch);
  }
}



