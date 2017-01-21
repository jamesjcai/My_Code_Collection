#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <list>
#include <cmath>
#include "population.h"
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

class Metapopulation{
 public:
  Metapopulation(const char[]);
  Metapopulation(Population);
  Metapopulation( vector<Population>);
  Metapopulation(const char genmapfile[], const char physmapfile[], const char hapfile[], int length);
  const Metapopulation &operator=(const Metapopulation&);
  vector<Population> pops;
  vector<vector<double> > migrates;
  void print_mrates();
  int split(int);
  void addpop(Population);
  bool resize(gsl_rng*, int, int);
  double add_selected_site(gsl_rng*, int);
  double add_selected_site(int, double);
  void set_all_mu(double);
  void set_all_r(double);
  double get_mrate(int, int);
  void set_mrate(int, int, double, bool);
  void migrate(int, int, int, int);
  void single_gen(gsl_rng*);
  void set_selected_site(double);
  void set_selection_coefs(int, double, double);
  Population get(int);
  bool evolve(gsl_rng*, int);
  void watchpop(int);
  double selected_site;
 private:
  set<int> towatch;
};
  
