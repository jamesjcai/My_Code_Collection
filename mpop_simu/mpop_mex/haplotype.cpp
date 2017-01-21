#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "haplotype.h"
#include <math.h>
using namespace std;

/*
code for Haplotype class

this is a linked list, with pointers to the first and last nodes (pointer to last node speeds up adding nodes to the end). The nodes are the positions of segregating sites. 

functions:
Haplotype() : constructor
void add_node_at_end( double) : self-explanatory
Node node_before(double)
void recombine(Haplotype, double)
Node node_after(double)
void print(): simply lists all the sites in the linked list, or prints "empty"
void add_node(double): self-explanatory. if the node is always going to the end, use add_node_at_end for speed
bool has_pos(double): returns true if a position is present in the list, false otherwise
bool empty(): true if there are no segregating sites
*/

Haplotype::Haplotype(){
  start_ptr.reset();
  end_ptr.reset();
}

Haplotype::Haplotype(const Haplotype& hap) throw(){
  start_ptr.reset();
  boost::shared_ptr<Node> temp;
  temp = hap.start_ptr;
  while(temp !=NULL){
    double val = temp->pos;
    add_node_at_end(val);
    temp = temp->next;
  }
}

Haplotype::~Haplotype(){
  boost::shared_ptr<Node> temp;
  temp = start_ptr;
  while (start_ptr != NULL){
    temp = start_ptr->next;
    start_ptr.reset();
    start_ptr = temp;
  }
  end_ptr.reset();
}

const Haplotype& Haplotype::operator=(const Haplotype& hap){
  start_ptr.reset();
  boost::shared_ptr<Node> temp;
  temp = hap.start_ptr;
  while(temp !=NULL){
    double val = temp->pos;
    //cout << val << "\n";
    add_node_at_end(val);
    temp = temp->next;
  }
  return *this;
}

void Haplotype::del(){
  boost::shared_ptr<Node> temp;
  temp = start_ptr;
  while (start_ptr != NULL){
    temp = start_ptr->next;
    start_ptr.reset();
    start_ptr = temp;
  }
  end_ptr.reset();
} 

boost::shared_ptr<Haplotype::Node> Haplotype::node_before(double pos){
  // returns a pointer to the node that comes before a given position. should be useful for recombination, etc. 
  //important note about node_before and node_after: they assume no position is given which is already in the list. if such a position is given, node_before will return the node before the node in question, and node after will return the node in question. this is equivalent to the passed argument being in front of the node. 
  boost::shared_ptr<Node> temp, temp2;
  temp = start_ptr;
  temp2 = temp;
  while( temp != NULL && temp->pos < pos){
    temp2 = temp;
    temp = temp->next;
  }
  return(temp2);
  
}
 
boost::shared_ptr<Haplotype::Node> Haplotype::node_after(double pos){
  boost::shared_ptr<Node> temp, temp2;
  temp = start_ptr;
  while( temp != NULL && temp->pos < pos){
    temp = temp->next;
  }
  return(temp);
}
  
void Haplotype::add_node_at_end(double position){

  //
  // 
  // void add_node_at_end(double position)
  //
  //
  // add a node to the end of the linked list
  //
  //
  boost::shared_ptr<Node> temp(new Node);
  boost::shared_ptr<Node> temp2;
  temp2 = end_ptr;
  temp->pos = position;
  temp->next.reset();

  // 
  // if empty, make the first value
  //

  if (start_ptr == NULL){
    start_ptr = temp;
    end_ptr = temp;
  }

  //
  // otherwise, add to the end
  //

  else{
    end_ptr = temp;
    temp2->next = temp;
  }
}

boost::shared_ptr<Haplotype::Node> Haplotype::get_start_ptr(){
  boost::shared_ptr<Node> tmp;
  tmp = start_ptr;
  return(tmp);
}

boost::shared_ptr<Haplotype::Node> Haplotype::get_end_ptr(){
  boost::shared_ptr<Node> tmp;
  tmp = end_ptr;
  return(tmp);
}
void Haplotype::print(){

  //
  //
  //
  // void print()
  // 
  // print the entire haplotype, including start_ptr, all the pointers
  // in the nodes, and end_ptr to stdout
  //
  //
  //
  //

  boost::shared_ptr<Node> temp;
  temp = start_ptr;
  cout << "start ptr " << start_ptr << "\n";
  
  // 
  // if empty, print empty
  //

  if (temp == NULL){
    cout << "empty\n";
  }

  //
  // otherwise, print
  //

  else{
    while (temp !=NULL){
      cout << temp << " " << temp->pos << " " << temp->next <<"\n";
      temp  = temp->next;
    }
  }
  cout << "end_ptr " << end_ptr << "\n";
}

void Haplotype::add_node( double pos){ 

  //
  //
  //
  // void add_node(double pos). adds a node at position pos
  //
  //
  //
  
  boost::shared_ptr<Node> temp(new Node);
  boost::shared_ptr<Node> temp2, temp3;
  temp->pos = pos;
  temp2 = start_ptr;
  temp3 = temp2;
  
  //
  // Case 1: the mutation is the only segregating site on the chromosome
  // or the mutation is past the last segregating site. this is equivalent 
  // to add_node_at_end
  //
  
  if (end_ptr == NULL || end_ptr->pos < pos){
    add_node_at_end(pos);
  }
 
  //
  // Case 2: The mutation occurs before the first segregating site
  //
  
  else if(temp2->pos > pos){
    //printf("case 2\n");
    temp->next = temp2;
    start_ptr = temp;
  }

  //
  // Case 3: the mutation occurs after the first segregating site
  //

  else{
    while(temp2 != NULL && temp2->pos<pos){
      temp3 = temp2;
      temp2 = temp2->next;
    }
    
    
    // 
    //  and the mutation occurs before the last segregating site
    //
    
    if (temp2->pos >pos){
      temp->next= temp2;
      temp3-> next= temp;
    }
    
    //
    // if the haplotype already has the node, print a warning and try again
    // with a different value (if this is happening a lot
    // check your random number generator
    //

    else if (has_pos(pos)){
      cerr << "WARNING: double hit at position " << pos << ". Adding altenate site\n";
      add_node(pos+0.0001);
    }
    
    //
    // in any other case, there's a problem. print the haplotype and exit
    //

    else{
      printf("ERROR: problem with indices in mutation %f %f\n", temp3->pos, pos);
      print();
      exit(1);
    }
  }  
}




bool Haplotype::has_pos (double pos){
  //
  //
  //
  // bool has_pos(double pos)
  // 
  // test to see if the haplotype has a given position 
  //
  //
  boost::shared_ptr<Node> temp;
  temp = start_ptr;
  double test = pos - 1E-15;
  while(temp != NULL && temp->pos < test){
    temp = temp->next;
  }
  if (temp != NULL && fabs(temp->pos -pos)<1E-15){
    return(true);
  }
  else{
    return(false);
  }
}



bool Haplotype::empty(){
  //
  //
  // bool empty()
  //
  // self-explanatory
  //
  //

  if (start_ptr == NULL){
    return(true);
  }
  else{
    return(false);
  }
}


void Haplotype::recombine(Haplotype *hap2, double pos){
  //
  //
  //
  // void recombine(Haplotype *hap2, double pos)
  //
  // recombination at position pos. switch the relevant pointers in the linked list
  //
  //

  //
  // first step: position iterators on either side of the "breakpoint"
  //

  boost::shared_ptr<Node> hap1_a, hap1_b, hap2_a, hap2_b;
  hap1_a = node_before(pos);
  hap1_b = node_after(pos);
  hap2_a = hap2->node_before(pos);
  hap2_b = hap2->node_after(pos);
  //cout << "hap1a " << hap1_a << " hap1_b "<< hap1_b << " hap2a " <<hap2_a << " hap2_b " <<hap2_b <<"\n";
  
  //
  // Case 1: no segregating sites on one of the haplotypes
  //
  // (case 0 is no segregating sites on either, in which case nothing is done)
  //
  //

  //
  // case 1a: no segregating sites on haplotype 1
  //

  if (empty() && !hap2->empty()){
    
    //
    // if the recombination is after the last segregating site, do nothing
    //

    if (hap2_b != NULL){
      // 
      // test to see if the breakpoint is before the first segregating site of hap2
      //
      start_ptr = hap2_b;
      end_ptr = hap2->end_ptr;
      if (hap2_a == hap2_b){
	hap2->end_ptr.reset();
	hap2->start_ptr.reset();
      }
      else{
	hap2->end_ptr = hap2_a;
	hap2_a->next.reset();
      }
    }
  }
  
  //
  // case 1b: no segregating sites on haplotype 2

  else if (hap2->empty() && !empty()){
    
    if (hap1_b !=NULL){
      hap2->start_ptr = hap1_b;
      hap2->end_ptr = end_ptr;
      if(hap1_a == hap1_b){
	end_ptr.reset();
	start_ptr.reset();
      }
      else{
	end_ptr = hap1_a;
	hap1_a->next.reset();
      }
    }
  }
  //
  //Case 2: recombination in the middle of the chromosomes
  //

   else if (hap1_b != NULL && hap2_b !=NULL && hap1_a->pos < pos && hap2_a->pos <pos){
     boost::shared_ptr<Node> temp;
     hap1_a->next = hap2_b;
     hap2_a->next = hap1_b;
     temp = end_ptr;
     end_ptr = hap2->end_ptr;
     hap2->end_ptr = temp;
   }
  
  //
  //Case 3: recombination at the end of one of the chromosomes
  //
  //3a
   else if (hap1_b == NULL && hap2_b !=NULL && hap2_a->pos < pos){
     hap1_a->next = hap2_b;
     end_ptr = hap2->end_ptr;
     hap2_a->next.reset();
     hap2->end_ptr = hap2_a;
   }
  //
  //4b
  //
   else if (hap2_b == NULL && hap1_b !=NULL && hap1_a->pos <pos){
     hap2_a->next = hap1_b;
     hap2->end_ptr = end_ptr;
     hap1_a->next.reset();
     end_ptr = hap1_a;
   }
  
  //
  // case 4: recombination before the first segregating site on one of the chromosomes
  //

  //
  //4a before first site on haplotype 1
  //

   else if (hap1_a !=NULL && hap1_a->pos > pos){

     //
     // in the middle of haplotype 2
     //

     if (hap2_b != NULL && hap2_a->pos <pos){
       boost::shared_ptr<Node> swp;
       swp = end_ptr;
       start_ptr = hap2_b;
       end_ptr = hap2->end_ptr;
       hap2->end_ptr = swp;
       hap2_a->next = hap1_a;
     }

     //
     // at the end of haplotype 2
     //
     else if (hap2_b == NULL){
       hap2_a->next = hap1_a;
       hap2->end_ptr = end_ptr;
       start_ptr.reset();
       end_ptr.reset();
     }
   }

  //
  //4b inverse of 4a
  //

   else if (hap2_a !=NULL && hap2_a->pos >pos){
     if (hap1_b != NULL && hap1_a->pos <pos){
       boost::shared_ptr<Node> swp;
       swp = hap2->end_ptr;
       hap2->start_ptr = hap1_b;
       hap2->end_ptr = end_ptr;
       end_ptr = swp;
       hap1_a->next = hap2_a;
     }
     else if (hap1_b == NULL){
       hap1_a->next = hap2_a;
       end_ptr = hap2->end_ptr;
       hap2->start_ptr.reset();
       hap2->end_ptr.reset();
     }
   }      
}
