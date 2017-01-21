#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <vector>
#include "metapopulation.h"
#include <list>
#include <math.h>
#include "haplotype.h"
using namespace std;


Metapopulation::Metapopulation(const char in[]){
  pops.clear();
  Population inpop(in);
  pops.push_back(inpop);
  migrates.clear();
  vector<double> tmp;
  tmp.push_back(0);
  migrates.push_back(tmp);
}

Metapopulation::Metapopulation(Population pop){
  pops.clear();
  pops.push_back(pop);
  migrates.clear();
  vector<double> tmp;
  tmp.push_back(0);
  migrates.push_back(tmp);
}

Metapopulation::Metapopulation(const char genmapfile[], const char physmapfile[], const char hapfile[], int length){
  Population pop(genmapfile, physmapfile,hapfile, length);
  pops.clear();
  pops.push_back(pop);
  migrates.clear();
  vector<double> tmp;
  tmp.push_back(0);
  migrates.push_back(tmp);
}

Metapopulation::Metapopulation( vector<Population> inpops){
  pops = inpops;
  selected_site = 0;
  int num = pops.size();
  int i, j;
  for (i = 0; i<num; i++){
    vector<double> tmp;
    for (j = 0; j<num; j++){
      tmp.push_back(0);
    }
    migrates.push_back(tmp);
  }
}

const Metapopulation& Metapopulation::operator=(const Metapopulation& met){
  pops.assign(met.pops.begin(), met.pops.end());
  migrates = met.migrates;
  selected_site = met.selected_site;
  towatch = met.towatch;
  return *this;
}

void Metapopulation::set_all_mu(double m){
  vector<Population>::iterator it;
  for (it  = pops.begin(); it !=pops.end(); it++){
    it->mu = m;
  }
}
void Metapopulation::set_all_r(double r){
  vector<Population>::iterator it;
  for (it  = pops.begin(); it !=pops.end(); it++){
    it->r = r;
  }
}

bool Metapopulation::resize(gsl_rng* r, int pop, int size){
  if (pop>=pops.size()){
    cerr << "cannot resize population "<<pop<<". max is "<<(pops.size()-1);
    exit(1);
  }
  pops.at(pop).resize(r, size);
  set<int>::iterator it;
  for (it = towatch.begin(); it !=towatch.end(); it++){
    int pop = *it;
    if (!pops.at(pop).has_site(selected_site)){
      return(false);
    }
  }
  return(true);
}
  
void Metapopulation::print_mrates(){
  int len = pops.size();
  int i, j;
    for (i = 0; i<len; i++){
      for (j = 0; j<len; j++){
	cout << migrates.at(i).at(j) << " ";
      }
      cout << "\n";
    }
}

void Metapopulation::set_selection_coefs(int pop, double s, double h){
   if (pop>=pops.size()){
    cerr << "cannot set selection on population "<<pop<<". max population is "<<(pops.size()-1);
    exit(1);
   }
   pops.at(pop).set_selection_coefs(s, h);
}
   

void Metapopulation::addpop(Population toadd){
  pops.push_back(toadd);
  vector<vector<double> >::iterator it;
  for (it = migrates.begin(); it !=migrates.end(); it++){
    it->push_back(0);
  }
  vector<double> tmp;
  tmp.assign(pops.size(), 0);
  migrates.push_back(tmp);
}

int Metapopulation::split(int which){
  if (which>=pops.size()){
    cerr << "cannot split population "<<which<<". max population is "<<(pops.size()-1);
    exit(1);
  }
  int index = pops.size();
  Population tmp = pops.at(which);
  addpop(tmp);
  if (towatch.find(which) != towatch.end()){
   watchpop(index);
  }
  return(index);  
}

double Metapopulation::add_selected_site(gsl_rng *r, int pop){
  if (pop>=pops.size()){
    cerr << "ERROR: cannot add selected site to population " <<pop<<". no such population. max is "<<(pops.size()-1);
    exit(1);
  }
  pops.at(pop).add_selected_site(r);
  return(0.5);
}

double Metapopulation::add_selected_site(int pop, double freq){
  if (pop>=pops.size()){
    cerr << "ERROR: cannot add selected site to population " <<pop<<". no such population. max is "<<(pops.size()-1);
    exit(1);
  }
  double site = pops.at(pop).add_selected_site(freq);
  return(site);
}

double Metapopulation::get_mrate(int from, int to){
  if (from >=migrates.size() || from >=migrates.size()){
    cerr << "ERROR: unable to get mrate from popualtion "<<from<< " to population "<<to<<". number of population is " << pops.size()<<"\n";
    exit(1);
  }
  else{
    return(migrates.at(from).at(to));
  }
}

void Metapopulation::set_mrate(int pop1, int pop2, double rate, bool sym){
  if (pop1 >=migrates.size() || pop2 >=migrates.size()){
    cerr << "ERROR: unable to set mrate in " << pop1 <<" " <<pop2 << ". max is "<<(migrates.size()-1)<<"\n";
    exit(1);
  }
  else{
    migrates.at(pop1).at(pop2) = rate;
    if (sym){
      migrates.at(pop2).at(pop1) = rate;
    }
  }
}

Population Metapopulation::get(int i){
  return pops.at(i);
}

void Metapopulation::single_gen(gsl_rng *r){
  int len = pops.size();
  int i, j;
  for (i = 0; i<len; i++){
    int N1 = (pops.at(i).get_N())/2;
    for (j = 0; j<len; j++){
      if (i ==j){continue;}
      int N2 = (pops.at(j).get_N())/2;
      //cout << N1 << " " << N2 << "\n";
      double rate = get_mrate(i,j);
      int migs = gsl_ran_poisson(r, (rate*N1));
      int k;
      for (k = 0; k<migs; k++){
	int migrating = gsl_rng_uniform_int(r, N1);
	int replaced = gsl_rng_uniform_int(r, N2);
	migrate(i, j, migrating, replaced);
      }
    }
  }
  vector<Population>::iterator it;
  for (it = pops.begin(); it !=pops.end(); it++){
    it->single_gen(r);
  }
}

void Metapopulation::migrate(int popfrom, int popto, int indfrom, int indto){
  if (popto>(pops.size()-1) || popfrom > (pops.size()-1)){
    cerr<< "ERROR: no such population " <<  popto << " " << popfrom << ". max is "<< (pops.size()-1)<< "\n";
    exit(1);
  }
  if (pops.at(popto).get_N()<=(2*indto)+1){
    cerr << "ERROR: cannot migrate to individual "<<indto << " in population "<<popto<<".Total number of haplotypes is "<<pops.at(popto).get_N()<<"\n";
    exit(1);
  }
  if (pops.at(popfrom).get_N()<=(2*indfrom)+1){
    cerr << "ERROR: cannot migrate from individual "<<indfrom << " in population "<<popfrom<<".Total number of haplotypes is "<<pops.at(popfrom).get_N()<< "\n";
    exit(1);
  }

  pops.at(popto).assign(pops.at(popfrom).at(2*indfrom), 2*indto);
  pops.at(popto).assign(pops.at(popfrom).at((2*indfrom)+1), (2*indto)+1);
}



void Metapopulation::set_selected_site(double site){
  vector<Population>::iterator it;
  for (it = pops.begin(); it !=pops.end(); it++){
    it->set_selected_site(site);
  }
  selected_site = site;
}

bool Metapopulation::evolve(gsl_rng *r, int gens){
  int i;
  for (i = 0; i < gens; i++){
    //cout << i << "\n";
    single_gen(r);
    if (fabs(selected_site-0) > 1E-9){
      int remainder = i % 100;
      if ( i <10 || remainder ==0 || i == (gens-1)){
	set<int>::iterator it;
	for (it = towatch.begin(); it !=towatch.end(); it++){
	  int pop = *it;
	  if (!pops.at(pop).has_site(selected_site)){
	    return(false);
	  }
	}
      }
    }
  }
  return(true);
}  

void Metapopulation::watchpop(int popnum){
  towatch.insert(popnum);
}
