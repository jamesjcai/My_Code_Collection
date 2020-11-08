class ll_pnt{

 private :

 typedef struct node{
	void *info;
	node  *seg;
 };

 node *Topleft;
 node *newTopleft;
 node *Topright;

 public :

 // constructoras
 ll_pnt();
 ~ll_pnt();
 void add(void *info);
 void addrev(void *info);
 // consultoras

 void resetpt(void **pt);
 void *noend(void *pt);
 void *llpt(void *pt);
 void advpt(void **pt);

 // modificadoras

 void modpt(void *pt,void *inf);

 };
