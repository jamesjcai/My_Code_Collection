class ll_flt{

 private :

 typedef struct node{
	float info;
	node  *seg;
 };

 node *Topleft;
 node *newTopleft;
 node *Topright;

 public :

 // constructoras
 ll_flt();
 ~ll_flt();
 void add(float info);
 void addrev(float info);


 // consultoras

 void resetpt(void **pt);
 void *noend(void *pt);
 float llpt(void *pt);
 void advpt(void **pt);

 // modificadoras

 void modpt(void *pt,float inf);

 };
