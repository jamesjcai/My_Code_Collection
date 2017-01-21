#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <list>
#include <cmath>
#include "haplotype.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/tuple/tuple.hpp>
using namespace std;

class Population{
 public:
  Population( int );
  Population( int, double);
  Population(const char[]);
  Population(const char[], const char[], const char[], int);
  Population(const Population&) throw();
  const Population &operator=(const Population&);
  map<double, vector<double> > get_tsnps(double, double);
  map<double, vector<double> > get_tsnps(double, double, set<double>);
  void print();
  void printms(const char[], int);
  void printms(const char[], int, set<double>);
  void add_site(int, double);
  void insert_site(int, double);
  void set_selection_coefs(double, double);
  bool evolve(gsl_rng*,int);
  double add_selected_site(gsl_rng*);
  void set_selected_site(double);
  void assign(Haplotype, int);
  double get_freq(double);
  bool has_site(double);
  void set_s(double);
  void set_h(double);
  double get_s();
  double get_h();
  double get_selected_site();
  vector<double> get_freq_spec();
  map<double, double> get_allele_freqs();
  boost::tuple<int, double, double, double> summary_stats();
  double add_selected_site(double);
  bool resize(gsl_rng*, int);
  int size();
  int get_N();
  Haplotype at(int);
  double mu, r;
  set<double> get_sites();
  void print_ldmat(char[], double);
  void print_map(const char[], int, int);
  void print_map(const char[], int, int, set<double>);
  map<double, int> get_counts();
  double read_genmap(const char[], int);
  void single_gen(gsl_rng*);
  double selected_startfreq;
  double recomb_position(gsl_rng*);
 private:
  double get_bg_rrate();
  void set_recomb_probs();
  double faywuh(double);
  double tajima_d(int, double);
  double apwd(map<double, int>);
  int seg_sites(map<double, int>);
  void mutate(gsl_rng*);
  int N;
  void print_oldpop();
  double selected_site,s,h;
  double w[3];
  vector<Haplotype> inpop;
  vector<Haplotype> oldpop;
  vector<boost::tuple<double, double, double> > hspots;
  vector<double> spotprobs;
};


   
