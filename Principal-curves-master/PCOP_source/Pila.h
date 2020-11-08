class pila{


  private:

	 typedef struct node{
		 void *pt;
		 node *seg;
	 };
	 node *top;
     node *n_top;
	 
  public:

//  constructora
	pila();
	void apilar(void *pt);
	void *desapilar();

// consultora
	int pila_buida();

};



