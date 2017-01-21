#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <vector>
#include "population.h"
#include <list>
#include <math.h>
#include "haplotype.h"
using namespace std;

Population::Population( int popsize) {
  N = popsize;
  inpop = vector<Haplotype>(N);
  oldpop = vector<Haplotype>(N);
  selected_site = 0;
  selected_startfreq = 0;
  mu = 0;
  r = 0;
  s = 0;
  h = 0;
  set_selection_coefs(s,h);
}

Population::Population( int popsize, double selected_site) {
  N = popsize;
  selected_site = selected_site;
  inpop = vector<Haplotype >(N);
  oldpop = vector<Haplotype >(N);
  selected_startfreq = 0;
  mu = 0;
  r = 0;
  s = 0;
  h = 0;
  set_selection_coefs(s,h);
}

Population::Population(const Population& pop) throw(){
  N = pop.N;
  inpop = vector<Haplotype>(N);
  oldpop = vector<Haplotype >(N);
  mu = pop.mu;
  r = pop.r;
  s = pop.s;
  h = pop.h;
  selected_site = pop.selected_site;
  selected_startfreq = pop.selected_startfreq;
  set_selection_coefs(s,h);
  hspots = pop.hspots;
  spotprobs = pop.spotprobs;
  int  i;
  for(i = 0 ; i<N; i++){
    inpop.at(i).start_ptr.reset();
    inpop.at(i).end_ptr.reset();
    inpop.at(i) = pop.inpop.at(i);
  }

}

const Population& Population::operator=(const Population& pop){
  N = pop.N;
  inpop = vector<Haplotype>(N);
  oldpop = vector<Haplotype >(N);
  mu = pop.mu;
  r = pop.r;
  s = pop.s;
  h = pop.h;
  selected_site = pop.selected_site;
  selected_startfreq = pop.selected_startfreq;
  set_selection_coefs(s,h);
  hspots = pop.hspots;
  spotprobs = pop.spotprobs;
  int  i;
  for(i = 0 ; i<N; i++){
    inpop.at(i).start_ptr.reset();
    inpop.at(i).end_ptr.reset();
    inpop.at(i) = pop.inpop.at(i);
  }
  return *this;
}
  
Population::Population(const char in[]){
  // 
  //
  //initialize the population from a file in ms format
  //
  //
  int i;
  selected_site = 0;
  selected_startfreq = 0;
  mu = 0;
  r = 0;
  s = 0;
  h = 0;
  bool mpop_flag = false;
  set_selection_coefs(s,h);
  int j = 0;
  string st;
  vector<double> sites;
  ifstream infile(in);
  if (infile.fail()){
    cerr << "ERROR: cannot open file " << in << "\n";
    exit(1);
  }

  //
  // the first 6 lines contain type of file, parameters, positions
  //

  for(i = 1; i<7; i++){
    string buf;
    getline(infile,st);
    stringstream ss(st);
    vector<string> line;
    while (ss >> buf){
      line.push_back(buf);
    }
    if (i == 1){

      //
      // first line contains mpop in first position if an mpop file 
      //
      string test = line.at(0);
      //
      //if this is an mpop file (as opposed to an ms file), look for the parameters
      //

      if (test == "mpop"){ 
	mpop_flag = true;
	selected_site = atof(line.at(2).c_str());
	selected_startfreq = atof(line.at(4).c_str());
	mu = atof(line.at(6).c_str());
	r = atof(line.at(8).c_str());
	s = atof(line.at(10).c_str());
	h = atof(line.at(12).c_str());
	set_selection_coefs(s,h);
	string test("hspots:");
	vector<string>::iterator stit;
	for (stit = line.begin(); stit != line.end(); stit++){
	  string tmp = *stit;
	  if (test == tmp){
	    vector<string>::iterator stit2;
	    for (stit2 = stit+1; stit2 != line.end(); stit2++){
	      string entry2 = *stit2;
	      boost::tuple<double, double, double> hspot;
	      double st = atof(entry2.c_str());
	      stit2++;
	      double ln = atof(stit2->c_str());
	      stit2++;
	      double intense = atof(stit2->c_str());
	      hspot = boost::make_tuple(st, ln, intense);
	      //cout << st << " " << ln << " " << intense << "\n";
	      hspots.push_back(hspot);
	    }
	    break;
	  }
	}
	set_recomb_probs();
      }

      //
      //
      // no matter whether it's ms or mpop, 
      // the first line of the file in the second position
      // has the number of chromosomes in the second position
      // 
      // initialize vectors
      //

      int num = atoi(line.at(1).c_str());
      N = num;
      inpop = vector<Haplotype >(N);
      oldpop = vector<Haplotype >(N);
      if (!mpop_flag){
	//
	// read in positions and intensities of hotspots, if they exist
	//
	string v("-v");
	int len;
	string r("-r");
	vector<string>::iterator stit;
	for (stit = line.begin(); stit != line.end(); stit++){
	  string entry = *stit;
	  const char* dash = "-";
	  //cout << entry << "\n";
	  if (entry == r){
	    vector<string>::iterator tmp = stit+2;
	    len = atoi(tmp->c_str());
	  }
	  else if (entry == v){
	    vector<string>::iterator stit2;
	    for (stit2 = stit+2; stit2 != line.end(); stit2++){
	      string entry2 = *stit2;
	      if (entry2[0] == *dash){
		break;
	      }
	      boost::tuple<double, double, double> hspot;
	      double start = atof(entry2.c_str());
	      stit2++;
	      double end = atof(stit2->c_str());
	      stit2++;
	      double intense = atof(stit2->c_str());
	      double st = start/(double) len;
	      double ln = (end-start)/(double)len;
	      hspot = boost::make_tuple(st, ln, intense);
	      //cout << st << " " << ln << " " << intense << "\n";
	      hspots.push_back(hspot);
	    }
	    break;
	  }
	}
	set_recomb_probs();
      }
    }

    else if (i == 6){
      
      //
      //the sixth line of the file has the list of positions
      //
      vector<string>::iterator it;
      for (it = line.begin(); it !=line.end(); it++){
	string val = *it;
	if (it !=line.begin()){
	  double num = atof(val.c_str());
	  sites.push_back(num);
	  //cout << val <<"\n";
	}
      }
    }
    ss.flush();
  }

  //
  //
  // ensure that positions are unique. 
  //
  //


  vector<double>::iterator it1;
  set<double> myset;
  pair<set<double>::iterator, bool> test;
  for (it1 = sites.begin(); it1 != sites.end(); it1++){
    double val = *it1;
    test = myset.insert(val);
    while( test.second == false){
      //
      // don't let the position of the selected site change. if there
      // are two mutations at the selected site, exit (fix this later)
      //
   
      val = val + 0.0001;
      *it1 = val;
      test = myset.insert(val);
    }
    
    //
    // special case if there are more than one 1s. don't want positions >1
    //
    
    while( val > 1 || test.second == false){
      val = val-0.00001;
	*it1 = val;
	test = myset.insert(val);
    }	
  }


  while(getline(infile, st)){

    //
    // the rest of the file contains the haplotypes
    //
    int k = 0;
    for (k = 0; k<st.size(); k++){
      char val = st.at(k);
      if (val == '1'){
	inpop.at(j).add_node(sites.at(k));
      }
    }
    j++;
  }
  infile.close();
}

Population::Population(const char genmapfile[], const char physmapfile[], const char hapfile[], int length){
  double rrate = read_genmap(genmapfile, length);
  inpop.clear();
  r = rrate;
  s = 0;
  h = 0;
  set_selection_coefs(s, h);
  selected_site =0;
  selected_startfreq = 0;
  mu = 0;
  vector<double> sites;
  ifstream pinfile(physmapfile);
  ifstream hinfile(hapfile);
  string st;
  if (pinfile.fail()){
    cerr << "ERROR: cannot open file " << physmapfile << "\n";
    exit(1);
  }
  getline(pinfile,st);
  while(getline(pinfile, st)){
    string buf;
    stringstream ss(st);
    vector<string> line;
    while (ss >> buf){
      line.push_back(buf);
    }
    double site = atof(line.at(2).c_str());
    site = site/(double)length;
    sites.push_back(site);
  }
  //cout<< sites.size() << "\n";
  while(getline(hinfile, st)){
    string buf;
    stringstream ss(st);
    vector<string> line;
    while (ss >> buf){
      line.push_back(buf);
    }
    Haplotype tmp;
    int i;
    for(i = 2; i< line.size(); i++){
      int pos = atoi(line.at(i).c_str());
      if (pos ==1){
	tmp.add_node(sites.at((i-2)));
      }
    }
    inpop.push_back(tmp);
  }
  N = inpop.size();
  oldpop = vector<Haplotype>(N);
}


double Population::read_genmap(const char in[], int totallength){
  double minrate = 1;
  vector<int> starts, lengths;
  vector<double> rates;
  int previous;
  double previous_rate;
  double totalmaplength = 0;
  vector<double>::iterator it;
  ifstream infile(in);
  string st;
  hspots.clear();
  if (infile.fail()){
    cerr << "ERROR: cannot open file " << in << "\n";
    exit(1);
  }
  getline(infile,st);
  string buf1;
  stringstream ss1(st);
  vector<string> line1;
  while (ss1 >> buf1){
    line1.push_back(buf1);
  }
  previous = atoi(line1.at(0).c_str());
  previous_rate = atof(line1.at(1).c_str());
  while(getline(infile, st)){
    string buf;
    stringstream ss(st);
    vector<string> line;
    while (ss >> buf){
      line.push_back(buf);
    }
    int pos = atoi(line.at(0).c_str());
    double rate = atof(line.at(1).c_str());
    int length = pos-previous;
    starts.push_back(previous);
    lengths.push_back(length);
    rates.push_back(previous_rate);
    if (rate < minrate){
      minrate = rate;
    }
    previous =pos;
    previous_rate = rate;
  }
  //cout << minrate << "\n";
  starts.push_back(previous);
  lengths.push_back(totallength-previous);
  rates.push_back(previous_rate);
  int i;
  for(i = 0; i< rates.size(); i++){
    totalmaplength+=(rates.at(i)*(double)lengths.at(i));
  }
  for(i = 0; i< rates.size(); i++){
    int start = starts.at(i);
    int length = lengths.at(i);
    double rate = rates.at(i);
    if (fabs(rate-minrate)>1E-9){
      //cout << "adding\n";
      int tmpint = (int)(rate/minrate);
      double newrate = (double)tmpint;
      boost::tuple<double, double, double> tmp = boost::make_tuple((double)start/(double)totallength, (double)length/(double)totallength, newrate);
      hspots.push_back(tmp);
    }
  }
  set_recomb_probs();
  return(totalmaplength);
}
void Population::print(){
  vector<Haplotype>::iterator it;
  Haplotype temp;
  for( it = inpop.begin(); it != inpop.end(); it++ ){
    temp  = *it;
    temp.print();
    cout << "\n";
  }
}

void Population::print_oldpop(){
  vector<Haplotype>::iterator theIterator;
  Haplotype temp;
  for( theIterator = oldpop.begin(); theIterator != oldpop.end(); theIterator++ ){
    temp  = *theIterator;
    temp.print();
    cout << "\n";
  }
}

void Population::set_selected_site(double site){
  selected_site = site;
}

int Population::get_N(){
  return(N);
}

double Population::get_s(){
  return(s);
}

double Population::get_h(){
  return(h);
}

double Population::get_selected_site(){
  return(selected_site);
}

Haplotype Population::at(int ind){
  if (ind > (N-1)){
    cerr << "No such haplotype "<<ind<< ". max is "<<(N-1);
    exit(1);
  }
  return(inpop.at(ind));
}

void Population::assign(Haplotype hap, int ind){
  inpop.at(ind) = hap;
}

void Population::set_s(double news){
  s = news;
  set_selection_coefs(s,h);
}

void Population::set_h(double newh){
  h = newh;
  set_selection_coefs(s,h);
}

void Population::printms(const char out[], int seed, set<double> sites){
  ofstream outfile(out);
  outfile.precision(15);
  set<double>::iterator it;
  vector<Haplotype >::iterator it2;
  outfile << "mpop " << N << " " << selected_site << " startfreq: "<< selected_startfreq<<" mu: " << mu << " r: "<< r<< " s: "<<s<<" h: " <<h;
  if (hspots.size() >0){
    outfile << " hspots:";
    vector<boost::tuple<double, double, double> >::iterator hit;
    for (hit = hspots.begin(); hit != hspots.end(); hit++){
      boost::tuple<double, double, double> tmp = *hit;
      outfile << " " << boost::get<0>(tmp) << " " << boost::get<1>(tmp) << " " << boost::get<2>(tmp);
    }
  }
  outfile << "\n"<<seed << "\n\n//\nderived_sites: ";
  int count = sites.size();
  outfile << count << "\npositions:";
  for (it = sites.begin(); it != sites.end(); it++){
    double forout = *it;
    outfile << " " << forout;
  }
  outfile << "\n";
  for (it2 = inpop.begin(); it2 !=inpop.end(); it2++){
    Haplotype temp = *it2;
    set<double>::iterator setit;
    for(setit = sites.begin(); setit !=sites.end(); setit++){
      double set_pos = *setit;
      if (temp.has_pos(set_pos)){
	outfile << "1";
      }
      else{
	outfile <<"0";
      }
    }
    outfile << "\n";
  }
}


void Population::printms(const char out[], int seed){
  set<double> sites;
  ofstream outfile(out);
  outfile.precision(15);
  outfile << "mpop " << N << " " << selected_site << " startfreq: "<< selected_startfreq<<" mu: " << mu << " r: "<< r<< " s: "<<s<<" h: " <<h;
  if (hspots.size() >0){
    outfile << " hspots:";
    vector<boost::tuple<double, double, double> >::iterator hit;
    for (hit = hspots.begin(); hit != hspots.end(); hit++){
      boost::tuple<double, double, double> tmp = *hit;
      outfile << " " << boost::get<0>(tmp) << " " << boost::get<1>(tmp) << " " << boost::get<2>(tmp);
    }
  }
  outfile << "\n"<<seed << "\n\n//\nderived_sites: ";
  sites = get_sites();
  
  set<double>::iterator it;
  vector<Haplotype >::iterator it2;
  int count = sites.size();
  outfile << count << "\npositions:";
  
  for (it = sites.begin(); it != sites.end(); it++){
    double forout = *it;
    outfile << " " << forout;
  }
  outfile << "\n";
  
  for (it2 = inpop.begin(); it2 !=inpop.end(); it2++){
    boost::shared_ptr<Haplotype::Node> hapit;
    Haplotype temp = *it2;
    hapit = temp.start_ptr;
    set<double>::iterator setit;
    for(setit = sites.begin(); setit !=sites.end(); setit++){
      double set_pos = *setit;
      if (hapit != NULL){
	double pos = hapit->pos;
	if (fabs(set_pos-pos) < 1E-15){
	  outfile << "1";
	  hapit = hapit->next;
	}
	else{
	  outfile << "0";
	}
      }
      else{
	outfile << "0";
      }
    }
    outfile << "\n";
  }
  
  outfile.close();
}

int Population::size(){
  return(N);
}


set<double> Population::get_sites(){
  set<double> sites;
  vector<Haplotype >::iterator it;
  for(it = inpop.begin(); it != inpop.end(); it++){
    boost::shared_ptr<Haplotype::Node> it2;
    Haplotype temp = *it;
    it2 = temp.start_ptr;
    while(it2 != NULL){
      double pos = it2->pos;
      sites.insert(pos);
      it2 = it2->next;
    }
  }
  return(sites);
}
void Population::set_selection_coefs(double news, double newh){
	double max1;
	double max2;
  s = news;
  h = newh;
  w[0] = 1;
  w[1] = (h*s) +1;
  w[2] = 1+s;
  //double max1 = w[0]>?w[1];
  //double max2 = max1>?w[2];
  if (w[0]>w[1]){
	  max1=w[0];
  }else{
	  max1=w[1];
  }
  if (max1>w[2]){
	  max2=max1;
  }else{
	  max2=w[2];
  }

  w[0] = w[0]/max2;
  w[1] = w[1]/max2;
  w[2] = w[2]/max2;
  //cout << w[0] <<" " << w[1] << " " <<w[2] << "\n";
}
  
void Population::add_site(int ind, double pos){
  inpop.at(ind).add_node(pos);
}

bool Population::has_site(double pos){
  vector<Haplotype>::iterator it;
  it = inpop.begin();
  bool found = false;
  while(found == false && it != inpop.end()){
    if (it->has_pos(pos)){
      found = true;
    }
    it++;
  }
  return(found);
}

bool Population::evolve(gsl_rng *ran, int ngen){
  //
  //
  // evolves the population for ngen generations or until the selected site is lost, whichever
  // comes first
  // returns true if the selected site is still present, false otw
  //
  //
  int i;
  for (i = 0; i<ngen; i++){
    single_gen(ran);

    if (selected_site > 0){
      //
      // if there's a selected site, check to see if it's been lost. check every generation
      // for 10 generations, then every 100 afterwards
      //
      if (i < 10 && has_site(selected_site)==false){
	return(false);
      }
	
	
      else{
	int remainder = i % 100;
	if (remainder == 0 && has_site(selected_site)==false){
	  return(false);
	}
      }
    }
    
  }
  if (selected_site>0 && has_site(selected_site)==false){
    return(false);
  }
  return(true);
}

double Population::add_selected_site(gsl_rng* ran){
  //add a new selected site
  if (selected_site>0){
    cerr << "ERROR: attempt to add a selected site when one already exists. Max one selected site per population\n";
    exit(1);
  }
  int ind =  gsl_rng_uniform_int(ran, N);
  double test_selected_site = 0.5;
  while(get_freq(test_selected_site)>0){
    test_selected_site = test_selected_site+0.0001;
  }
  inpop.at(ind).add_node(test_selected_site);
  selected_site = test_selected_site;
  selected_startfreq = 1/(double)N; 
}

double Population::add_selected_site(double freq){
  //
  // adds a selected site of freqency freq. 
  // find the closest frequency SNP within 0.2 of the center of the chromosome
  //
  if (selected_site >0){
    cerr <<  "ERROR: attempt to add a selected site when one already exists. Max one selected site per population\n";
    exit(1);
  }
  
  double bestsite = 0;
  double bestfreq = 0;
  map<double, int> counts = get_counts();
  map<double, int>::iterator it;
  for (it = counts.begin(); it != counts.end(); it++){
    pair<double, int> der = *it;
    double pos = der.first;
    int count = der.second;
    double test_freq = (double) count / (double) N;
    //
    // if the frequency closer to the ideal and close to the middle of the chromosome
    // accept it (distance for close to middle is hard-coded as between 0.3 and 0.7
    //
    if (fabs(test_freq-freq)< fabs(bestfreq-freq) && fabs(pos-0.5)<0.2){
      bestfreq = test_freq;
      bestsite = pos;
    }
  }
  selected_site = bestsite;
  selected_startfreq = bestfreq;
  return(selected_site);

}

double Population::get_freq(double site){
  int count = 0;
  vector<Haplotype>::iterator it;
  for (it = inpop.begin(); it != inpop.end(); it++){
    Haplotype temp = *it;
    if (temp.has_pos(site)){
      count++;
    }
  }
  double freq = count/(double)N;
  return(freq);
}

void Population::mutate(gsl_rng *ran){  
  int muts = gsl_ran_poisson(ran, (N*mu));
  int i;
  for(i = 0; i<muts; i++){
    int ind = gsl_rng_uniform_int(ran, N);
    double pos = gsl_rng_uniform(ran);
    while(fabs(pos-selected_site)<1E-14){
      pos = gsl_rng_uniform(ran);
    }
    add_site(ind, pos);
  }
}

bool Population::resize(gsl_rng *ran, int newN){
  vector<Haplotype> temppop(newN);
  int i;
  for (i=0; i<newN; i++){
    int ind = gsl_rng_uniform_int(ran,N);
    temppop.at(i) = inpop.at(ind);
  }
  for (i = 0; i<N; i++){
    inpop.at(i).del();
    oldpop.at(i).del();
  }
  inpop.resize(newN);
  oldpop.resize(newN);
  N = newN;
  inpop = temppop;
  if(selected_site > 0 && fabs(get_freq(selected_site))< 1E-9){
    return false;
  }
  else{
    return(true);
  }
}

void Population::set_recomb_probs(){
  double total_len = 0;
  spotprobs.clear();
  double hspotlen = 0;
  vector<boost::tuple<double, double, double> >::iterator it;
  for (it = hspots.begin(); it != hspots.end(); it++){
    boost::tuple<double, double, double> tmp = *it;
    double len = boost::get<1>(tmp);
    double inten = boost::get<2>(tmp);
    total_len += (len*inten);
    hspotlen +=len;
  }
  total_len += (1-hspotlen);
  //cout << "total: " << total_len << " hspot: "<<hspotlen<< "\n"; 
  for (it = hspots.begin(); it != hspots.end(); it++){
    boost::tuple<double, double, double> tmp = *it;
    double len = boost::get<1>(tmp);
    double inten = boost::get<2>(tmp);
    double prob = (len*inten)/total_len;
    spotprobs.push_back(prob);
    //cout << prob << "\n";
  }
  
}


double Population::recomb_position(gsl_rng *ran){
  double begin = gsl_rng_uniform(ran);
  double total_prob = 0;
  bool inspot = false;
  int i;
  for (i = 0; i < spotprobs.size(); i++){
    total_prob += spotprobs.at(i);
    if (begin< total_prob){
      inspot = true;
      double start = boost::get<0>(hspots.at(i));
      double len = boost::get<1>(hspots.at(i));
      double second = gsl_rng_uniform(ran);
      double pos = start + (len*second);
      return(pos);
    }
  }
  if (!inspot){
    bool inspot2 = true;
    double pos;
    while(inspot2 == true){
      bool inspot3 = false;
      pos =gsl_rng_uniform(ran);
      vector<boost::tuple<double, double, double> >::iterator it;
      for (it = hspots.begin(); it != hspots.end(); it++){
	double hpos = boost::get<0>(*it);
	double end = hpos+ boost::get<1>(*it);
	if ((pos > hpos)  && (pos < end)){
	  inspot3 = true;
	}
      }
      inspot2 = inspot3;
    }
    return(pos);
  }    
}

void Population::single_gen(gsl_rng *ran){
  //
  // a single generation.
  // 
  // the most time consuming part of the program is the copying of haplotypes   
  // to try and reduce the number of times it's done, do as many shallow copies as possible
  // 
  // returns a vector containing 1 if the haplotype was shallow copied, 0 otw
  // shallow_copied is this vector from the previous generation
  //

  //
  // 1. mutate 
  //
  mutate(ran); 
  //
  // vector to keep track of which haplotypes can't be shallow copied
  // (either there has been recombination, or has already been sampled)
  //
  // 
  vector<int> no_shallow_cp(N);
  no_shallow_cp.assign(N,0);
  //
  //swap pointers,
  //delete haplotypes that hadn't been shallow-copied
  //reset pointers for those that had
  //
  //
  oldpop.swap(inpop);
  int clear;
  for (clear = 0; clear < N ; clear++){
    inpop.at(clear).start_ptr.reset();
    inpop.at(clear).end_ptr.reset();
  }
  int i =0;
  while (i < (N/2)){
    //
    // 1. choose two individuals to be the possible parents 
    //

    int selsite_count = 0;     
    int par1 = gsl_rng_uniform_int(ran, (N/2));
    int par2 = gsl_rng_uniform_int(ran, (N/2));
    while (par1 == par2){
     par2 = gsl_rng_uniform_int(ran, (N/2));
    }
    int par1_hap1 = 2*par1;
    int par1_hap2 = (2*par1)+1;
    int par2_hap1 = 2*par2;
    int par2_hap2 = (2*par2)+1;
    //cout << par1_hap1 << " " << par1_hap2 << " " <<par2_hap1 << " " <<par2_hap2<< "\n";
    vector<double> p_recombs1;
    vector<double> p_recombs2;
    //
    //recombine the parents, keeping track of the positions of the recombinations
    //
    if (r>0){
      int n_recombs1 = gsl_ran_poisson(ran, r);
      int n_recombs2 = gsl_ran_poisson(ran, r);
      if (n_recombs1>0){
	no_shallow_cp[par1_hap1] = 1;
	no_shallow_cp[par1_hap2] = 1;
      }
      if (n_recombs2>0){
	no_shallow_cp[par2_hap1] = 1;
	no_shallow_cp[par2_hap2] = 1;
      }
      int j;
      for (j= 0; j<n_recombs1; j++){
	double pos = recomb_position(ran);
	p_recombs1.push_back(pos);
	//cout <<  i<< " "<<par1_hap1 << " " <<par1_hap2 << " "<<pos << "\n";
	oldpop.at(par1_hap1).recombine(&oldpop.at(par1_hap2), pos);
      }
      for (j= 0; j<n_recombs2; j++){
	double pos = recomb_position(ran);
	p_recombs2.push_back(pos);
	oldpop.at(par2_hap1).recombine(&oldpop.at(par2_hap2), pos);
      }
    }
    //
    //now choose which of the two haplotypes from each parents makes the proposed child
    //
     int hap1 = 0; 
    if (gsl_rng_uniform(ran) > 0.5){
      hap1 = 1;
    }
    int hap2 = 0; 
    if (gsl_rng_uniform(ran) > 0.5){
      hap2 = 1;
    }
    //
    // count the occurrences of the selected site
    //
    if (hap1 == 0 && oldpop.at(par1_hap1).has_pos(selected_site)){
      selsite_count++;
    }
    else if (hap1 == 1 && oldpop.at(par1_hap2).has_pos(selected_site)){
      selsite_count++;
    }
    if (hap2 == 0 && oldpop.at(par2_hap1).has_pos(selected_site)){
      selsite_count++;
    }
    else if (hap2 == 1 && oldpop.at(par2_hap2).has_pos(selected_site)){
      selsite_count++;
    }
    //
    //generate a random U(0,1), accept the individual according to the acceptance probs
    //
    double test = gsl_rng_uniform(ran);
	
    if (w[selsite_count]> test){
      //cout << w[0] <<" " << w[1] << " " <<test << "\n";
      //
      //the individual is accepted, put into the next generation
      //
      if(hap1 == 0){
	if (no_shallow_cp[par1_hap1] == 1){
	  inpop.at(2*i) = oldpop.at(par1_hap1);
	}
	else{
	  (inpop.at(2*i)).start_ptr = (oldpop.at(par1_hap1)).start_ptr;
	  (inpop.at(2*i)).end_ptr = (oldpop.at(par1_hap1)).end_ptr;
	  no_shallow_cp[par1_hap1] = 1;
	  //been_shallow_cp[2*i] = 1;
       }
      }
      else{
	if (no_shallow_cp[par1_hap2] == 1){
	  inpop.at(2*i) = oldpop.at(par1_hap2);
	}
	else{
	  (inpop.at(2*i)).start_ptr = (oldpop.at(par1_hap2)).start_ptr;
	  (inpop.at(2*i)).end_ptr = (oldpop.at(par1_hap2)).end_ptr;
	  no_shallow_cp[par1_hap2] = 1;
	  //been_shallow_cp[2*i] = 1;
	}
      }
      
      
      if( hap2 == 0){
	if (no_shallow_cp[par2_hap1] == 1){
	  inpop.at((2*i)+1) = oldpop.at(par2_hap1);
	}
	else{
	  (inpop.at((2*i)+1)).start_ptr = (oldpop.at(par2_hap1)).start_ptr;
	  (inpop.at((2*i)+1)).end_ptr = (oldpop.at(par2_hap1)).end_ptr;
	  no_shallow_cp[par2_hap1] = 1;
	  //been_shallow_cp[(2*i)+1] = 1;
	}
      }
      else{
	if (no_shallow_cp[par2_hap2] == 1){
	  inpop.at((2*i)+1) = oldpop.at(par2_hap2);
	}
	else{
	  (inpop.at((2*i)+1)).start_ptr = (oldpop.at(par2_hap2)).start_ptr;
	  (inpop.at((2*i)+1)).end_ptr = (oldpop.at(par2_hap2)).end_ptr;
	  no_shallow_cp[par2_hap2] = 1;
	  // been_shallow_cp[(2*i)+1] = 1;
	}
      }
      //
      // increment
      //

      i++;
    }
    //
    //undo the recombination of the parents so that they can be resampled
    //
    vector<double>::iterator it_par1, it_par2;
    for (it_par1 = p_recombs1.begin(); it_par1 != p_recombs1.end(); it_par1++){
      double val = *it_par1;
      oldpop.at(par1_hap1).recombine(&oldpop.at(par1_hap2), val);
    }
    for (it_par2 = p_recombs2.begin(); it_par2 != p_recombs2.end(); it_par2++){
      double val = *it_par2;
      oldpop.at(par2_hap1).recombine(&oldpop.at(par2_hap2), val);
    }
  }
  // return(been_shallow_cp);
}
//
//
//
//
//
//
// Summary Statistics
//
//
//
//
//
//
boost::tuple<int, double, double, double> Population::summary_stats(){
  //
  // returns the number of segregating sites, average pairwise differences, and tajima's D
  //
  map<double, int> counts = get_counts();
  int seg = seg_sites(counts);
  double pwd = apwd(counts);
  double t_d = tajima_d(seg, pwd);
  double h = faywuh(pwd);
  return(boost::make_tuple(seg, pwd, t_d, h));
}

double Population::tajima_d(int seg, double pwd){
  //
  // returns tajima's D in the population
  //
  double a1, a2, b1, b2, c1, c2, e1, e2, var, d;
  a1 = 0;
  a2 = 0;
  int i;
  for (i = 1; i<N; i++){
    a1 = a1 + 1/(double)i;
    a2 = a2 + 1/((double)i *(double)i);
  }
  b1 = ((double)N+1)/(3*((double)N-1));
  b2 = 2*((double)N*(double)N +(double)N+3)/(9*(double)N*((double)N-1));
  c1 = b1-(1/a1);
  c2 = b2 -(((double)N+2)/(a1*(double)N)) + (a2/(a1*a1));
  e1 = c1/a1;
  e2 = c2/(a1*a1 + a2);
  var = e1*(double) seg + e2 *(double)seg *((double)seg-1);
  if (seg==0){
    d = 1000;
  }
  else{
    d = (pwd - ((double)seg/a1))/sqrt(var);
  }
  return(d);
}

double Population::faywuh(double apwd){
  map<double, int> counts = get_counts();
  int nchrom = get_N();
  int i;
  double theta_h = 0;
  for (i = 1; i<nchrom; i++){
    int si= 0;
    map<double, int>::iterator it;
    for (it = counts.begin(); it !=counts.end(); it++){
      pair<double, int> tmp = *it;
      if (tmp.second ==i){
	si++;
      }
    }
    double add = (2*(double)si*(double)i*(double)i)/((double)nchrom*(double)(nchrom-1));
    theta_h = theta_h+add;
  }
  double h = apwd = theta_h;
  return(h);
}
  
vector<double> Population::get_freq_spec(){
  map<double, int> counts = get_counts();
  vector<double> spec;
  int segsites = seg_sites(counts);
  spec.assign((N-1), 0);
  map<double, int>::iterator it;
  for (it = counts.begin(); it != counts.end(); it++){
    pair<double, int> der = *it;
    double pos = der.first;
    int count = der.second;
    if (count <N){
      int countm1 = count -1;
      spec.at(countm1) = spec.at(countm1)+ 1/(double)segsites;
    }
  }
  return(spec);
}

map<double, double> Population::get_allele_freqs(){
  map<double, int> counts = get_counts();
  map<double, double> freqs;
  map<double, int>::iterator it;
  for (it = counts.begin(); it != counts.end(); it++){
    pair<double, int> der = *it;
    double pos = der.first;
    int count = der.second;
    double freq= (double) count/ (double)N;
    freqs[pos] = freq;
  }
  return(freqs);
}

map<double, int> Population::get_counts(){
  //
  // get_counts() returns a map of positions to the number of times each position appears
  //
  map<double, int> counts;
  set<double> sites = get_sites();
  set<double>::iterator it;
  for (it = sites.begin(); it != sites.end(); it++){
    double pos = *it;
    counts[pos] = 0;
  }
  vector<Haplotype>::iterator it2;
  for (it2 = inpop.begin(); it2 != inpop.end(); it2++){
    Haplotype test = *it2;
    boost::shared_ptr<Haplotype::Node> ptr;
    ptr = test.start_ptr;
    while (ptr !=NULL){
      double pos = ptr->pos;
      counts[pos] = counts[pos]+1;
      ptr = ptr->next;
    }
  }
  return(counts);
}

int Population::seg_sites(map<double, int> counts){
  //
  //
  // seg_sites() returns the number of segregating sites in a populstion
  //
  //
  int seg = 0;
  map<double, int>::iterator it;
  for (it = counts.begin(); it != counts.end(); it++){
    pair<double, int> der = *it;
    double pos = der.first;
    int count = der.second;
    if (count < N){
      seg++;
    }
  }
  return(seg);
}
  

double Population::apwd(map<double, int> counts){
  //
  // returns average number of pairwise differences
  //
  int sum = 0;
  double pr = ((double)N*((double)N-1))/2;
  map<double, int>::iterator it;
  for (it = counts.begin(); it != counts.end(); it++){
    pair<double, int> der = *it;
    double pos = der.first;
    int count = der.second;
    int increment = count * (N-count);
    sum = sum+increment;
  }
  double pwd = (double)sum/pr;
  return(pwd);
}

//
//
//
// Tagging
//
//
//
//
//
//
map<double, vector<double> > Population::get_tsnps(double maf, double r2_thresh){
  map<double, vector<double> > tsnps;
  map<double, int> counts = get_counts();
  set<double> overmaf;
  set<int> ind_overmaf;
  map<double, int>::iterator it;
  set<double>::iterator it2, it6;
  map<double, int> pos2index;
  map<int, double> index2pos;
  vector<Haplotype>::iterator it3;
  int size;
  //
  // 1. get a vector of all the positions with MAF> maf
  //
  for (it = counts.begin(); it != counts.end(); it++){
    pair<double, int> der = *it;
    double pos = der.first;
    int count = der.second;
    double freq = count/(double)N;
    if (freq>maf && (1-freq)>maf){
      overmaf.insert(pos);
      //cout << pos << "\n";
      int tmp = overmaf.size();
      pos2index[pos] = (tmp-1);
      index2pos[tmp-1] = pos;
      ind_overmaf.insert(tmp-1);
    }
  }
  //
  // 2. generate matrix of r2 values
  //

  //
  // 2a. initialize a vector of vectors
  //
  size = overmaf.size();
  vector<vector<double> > ld_mat(size);
  int i, j;
  for (i = 0; i<(size); i++){
    vector<double> temp;
    temp.assign(size, 0);
    ld_mat.at(i) = temp;
  }
  //
  // 2b. count up how many times each pairwise haplotype occurs 
  //

  for (it3 = inpop.begin(); it3 != inpop.end(); it3++){
    boost::shared_ptr<Haplotype::Node> ptr;
    ptr = it3->start_ptr;
    while(ptr != NULL && ptr->next !=NULL){
      double pos1 = ptr->pos;
      if (overmaf.find(pos1) != overmaf.end()){
	boost::shared_ptr<Haplotype::Node> ptr2;
	ptr2 = ptr->next;
	while(ptr2 != NULL){
	  //cout << "here\n";
	  double pos2 = ptr2->pos;
	  if (overmaf.find(pos2) != overmaf.end()){
	    int index1 = pos2index[pos1];
	    int index2 = pos2index[pos2];
	    ld_mat.at(index1).at(index2) = ld_mat.at(index1).at(index2)+1;
	    ld_mat.at(index2).at(index1) = ld_mat.at(index2).at(index1)+1;
	    ptr2 = ptr2->next;
	  }
	  else{
	    ptr2 = ptr2->next;
	  }
	}
	ptr = ptr->next;
      }
      else{
	ptr = ptr->next;
      }
    }
  }
  //
  // 2c. convert counts to r2
  //
  for (it2 = overmaf.begin(); it2 != overmaf.end(); it2++){
    for (it6 = it2; it6 != overmaf.end(); it6++){
      if (it2 != it6){
	double pos1 = *it2;
	double pos2 = *it6;;
	double freq1 = (double)counts[pos1]/(double)N;
	double freq2 = (double)counts[pos2]/(double)N;
	int index1 = pos2index[pos1];
	int index2 = pos2index[pos2];
	double both = ld_mat.at(index1).at(index2);
	double d = ((double)both/(double)N) - freq1*freq2;
	double r2 = (d*d)/ (freq1*freq2*(1-freq1)*(1-freq2));
	ld_mat.at(index1).at(index2) = r2;
	ld_mat.at(index2).at(index1) = r2;
      }
    }
  }
  //
  // 3.add tsnps until all snps have a tsnp with r2 > threshold 
  //
  while(overmaf.size() >  1){
    //cout <<"entering\n";
    //
    // count up the number of snps tagged by each snp
    //
    map<double, vector<double> > overt;
    set<double>::iterator it7;
    for (it7 = overmaf.begin(); it7 != overmaf.end(); it7++){
      int sum = 0;
      double pos= *it7;
      int index = pos2index[pos];
      vector<double> tmp = ld_mat.at(index);
      set<int>::iterator tmp2;
      vector<double> over;
      for (tmp2 = ind_overmaf.begin(); tmp2 != ind_overmaf.end(); tmp2++){
	double val = tmp.at(*tmp2);
	if (val>r2_thresh){
	  double inpos = index2pos[*tmp2];
	  over.push_back(inpos);
	}
      }
      overt.insert(make_pair(pos, over));
    }
    //
    //
    // find the snp that tags the most snps. if there is none, add all the snps to the tsnps
    //
    //
    int max = 0;
    double maxval = 0;
    vector<double> maxtagged;
    map<double, vector<double> >::iterator mapit;
    for(mapit = overt.begin(); mapit !=overt.end(); mapit++){
      if (mapit->second.size()> max){
	//cout << mapit->second << "\n";
	max = mapit->second.size();
	maxtagged = mapit->second;
	maxval = mapit->first;
      }
    }
    //cout << "max: " <<  max << " maxval: " << maxval << "\n";
    if (max == 0){
      set<double>::iterator it;
      for (it = overmaf.begin(); it != overmaf.end(); it++){
	double tmppos = *it;
	tsnps.insert(make_pair(tmppos, maxtagged));
      }
      break;
    }
    else{
      //cout << "maxpos " << maxpos  << "\n";
      tsnps.insert(make_pair(maxval, maxtagged));
      int index = pos2index[maxval];
      set<double>::iterator it4;
      set<double> toremove;
      toremove.insert(maxval);
      for (it4 = overmaf.begin(); it4 != overmaf.end(); it4++){
	double pos = *it4;
	int ind = pos2index[pos];
	if(ld_mat.at(index).at(ind) > r2_thresh){
	  toremove.insert(pos);
	}
      }
      set<double>::iterator it5;
      for (it5 = toremove.begin(); it5 != toremove.end(); it5++){
	double pos = *it5;
	int index = pos2index[pos];
	set<double>::iterator tmp = overmaf.find(pos);
	set<int>::iterator tmp2 = ind_overmaf.find(index);
	ind_overmaf.erase(tmp2);
	overmaf.erase(tmp);
      }
    }
  }
  return(tsnps);
}


map<double, vector<double> > Population::get_tsnps(double maf, double r2_thresh, set<double> remove){
  map<double, vector<double> > tsnps;
  map<double, int> counts = get_counts();
  set<double> overmaf;
  set<int> ind_overmaf;
  map<double, int>::iterator it;
  set<double>::iterator it2, it6;
  map<double, int> pos2index;
  map<int, double> index2pos;
  vector<Haplotype>::iterator it3;
  int size;
  //
  // 1. get a vector of all the positions with MAF> maf
  //
  for (it = counts.begin(); it != counts.end(); it++){
    pair<double, int> der = *it;
    double pos = der.first;
    int count = der.second;
    double freq = count/(double)N;
    if (freq>maf && (1-freq)>maf){
      overmaf.insert(pos);
      //cout << pos << "\n";
      int tmp = overmaf.size();
      pos2index[pos] = (tmp-1);
      index2pos[tmp-1] = pos;
      ind_overmaf.insert(tmp-1);
    }
  }
  //
  // 2. generate matrix of r2 values
  //

  //
  // 2a. initialize a vector of vectors
  //
  size = overmaf.size();
  vector<vector<double> > ld_mat(size);
  int i, j;
  for (i = 0; i<(size); i++){
    vector<double> temp;
    temp.assign(size, 0);
    ld_mat.at(i) = temp;
  }
  //
  // 2b. count up how many times each pairwise haplotype occurs 
  //

  for (it3 = inpop.begin(); it3 != inpop.end(); it3++){
    boost::shared_ptr<Haplotype::Node> ptr;
    ptr = it3->start_ptr;
    while(ptr != NULL && ptr->next !=NULL){
      double pos1 = ptr->pos;
      if (overmaf.find(pos1) != overmaf.end()){
	boost::shared_ptr<Haplotype::Node> ptr2;
	ptr2 = ptr->next;
	while(ptr2 != NULL){
	  //cout << "here\n";
	  double pos2 = ptr2->pos;
	  if (overmaf.find(pos2) != overmaf.end()){
	    int index1 = pos2index[pos1];
	    int index2 = pos2index[pos2];
	    ld_mat.at(index1).at(index2) = ld_mat.at(index1).at(index2)+1;
	    ld_mat.at(index2).at(index1) = ld_mat.at(index2).at(index1)+1;
	    ptr2 = ptr2->next;
	  }
	  else{
	    ptr2 = ptr2->next;
	  }
	}
	ptr = ptr->next;
      }
      else{
	ptr = ptr->next;
      }
    }
  }
  //
  // 2c. convert counts to r2
  //
  for (it2 = overmaf.begin(); it2 != overmaf.end(); it2++){
    for (it6 = it2; it6 != overmaf.end(); it6++){
      if (it2 != it6){
	double pos1 = *it2;
	double pos2 = *it6;;
	double freq1 = (double)counts[pos1]/(double)N;
	double freq2 = (double)counts[pos2]/(double)N;
	int index1 = pos2index[pos1];
	int index2 = pos2index[pos2];
	double both = ld_mat.at(index1).at(index2);
	double d = ((double)both/(double)N) - freq1*freq2;
	double r2 = (d*d)/ (freq1*freq2*(1-freq1)*(1-freq2));
	ld_mat.at(index1).at(index2) = r2;
	ld_mat.at(index2).at(index1) = r2;
      }
    }
  }



  set<double>::iterator sit;
  for (sit = remove.begin(); sit != remove.end(); sit++){
    double togo = *sit;
    if(overmaf.find(togo) != overmaf.end()){
      int index = pos2index[togo];
      set<double>::iterator it4;
      set<double> toremove;
      toremove.insert(togo);
      for (it4 = overmaf.begin(); it4 != overmaf.end(); it4++){
	double pos = *it4;
	int ind = pos2index[pos];
	if(ld_mat.at(index).at(ind) > r2_thresh){
	  toremove.insert(pos);
	}
      }
      set<double>::iterator it5;
      for (it5 = toremove.begin(); it5 != toremove.end(); it5++){
	double pos = *it5;
	int index1 = pos2index[pos];
	set<double>::iterator tmp = overmaf.find(pos);
	set<int>::iterator tmp2 = ind_overmaf.find(index1);
	ind_overmaf.erase(tmp2);
	overmaf.erase(tmp);
      }
    }
  }

   //
  //
  // 3.add tsnps until all snps have a tsnp with r2 > threshold 
  //
  while(overmaf.size() >  1){

    //cout <<"entering\n";
    //
    // count up the number of snps tagged by each snp
    //
    map<double, vector<double> > overt;
    set<double>::iterator it7;
    for (it7 = overmaf.begin(); it7 != overmaf.end(); it7++){

      int sum = 0;
      double pos= *it7;
      int index = pos2index[pos];
      vector<double> tmp = ld_mat.at(index);
      set<int>::iterator tmp2;
      vector<double> over;
      //cout << ind_overmaf.size()<< "\n";
      for (tmp2 = ind_overmaf.begin(); tmp2 != ind_overmaf.end(); tmp2++){

	double val = tmp.at(*tmp2);
	//cout << "here " << *tmp2 << " "<<val << "\n";
	if (val>r2_thresh){
	  double inpos = index2pos[*tmp2];
	  over.push_back(inpos);
	}
      }
      overt.insert(make_pair(pos, over));
    }
    //
    //
    // find the snp that tags the most snps. if there is none, add all the snps to the tsnps
    //
    //
    int max = 0;
    double maxval = 0;
    vector<double> maxtagged;
    map<double, vector<double> >::iterator mapit;
    for(mapit = overt.begin(); mapit !=overt.end(); mapit++){
      if (mapit->second.size()> max){
	//cout << mapit->second << "\n";
	max = mapit->second.size();
	maxtagged = mapit->second;
	maxval = mapit->first;
      }
    }
    //cout << "max: " <<  max << " maxval: " << maxval << "\n";
    if (max == 0){
      set<double>::iterator it;
      for (it = overmaf.begin(); it != overmaf.end(); it++){
	double tmppos = *it;
	tsnps.insert(make_pair(tmppos, maxtagged));
      }
      break;
    }
    else{
      //cout << "maxpos " << maxpos  << "\n";
      tsnps.insert(make_pair(maxval, maxtagged));
      int index = pos2index[maxval];
      set<double>::iterator it4;
      set<double> toremove;
      toremove.insert(maxval);
      for (it4 = overmaf.begin(); it4 != overmaf.end(); it4++){
	double pos = *it4;
	int ind = pos2index[pos];
	if(ld_mat.at(index).at(ind) > r2_thresh){
	  toremove.insert(pos);
	}
      }
      set<double>::iterator it5;
      for (it5 = toremove.begin(); it5 != toremove.end(); it5++){
	double pos = *it5;
	int index1 = pos2index[pos];
	set<double>::iterator tmp = overmaf.find(pos);
	set<int>::iterator tmp2 = ind_overmaf.find(index1);
	ind_overmaf.erase(tmp2);
	overmaf.erase(tmp);
      }
    }
  }
  return(tsnps);
}

//
// 
// now, the same algorthm as above, but with some snps removed
//


void Population::print_ldmat(char out[], double maf){
  ofstream outfile(out);
  outfile.precision(15);
  set<double> sites;
  map<double, int> pos2index;
  map<double, int> counts = get_counts();
  map<double, int>::iterator mapit;

  int i, j;
  set<double>::iterator it2, it6;
  vector<Haplotype>::iterator it3;

  int tmp = 0;
  for (mapit = counts.begin(); mapit !=counts.end(); mapit++){
    pair<double, int> der = *mapit;
    double pos = der.first;
    int count = der.second;
    double freq = count/(double)N;
    if (freq>maf && (1-freq)>maf){
      int tmp = sites.size();
      pos2index[pos] = tmp;
      sites.insert(pos);
    }
  }
  int size = sites.size();
  vector<vector<double> > ld_mat(size);

  for (i = 0; i<(size); i++){
    vector<double> temp;
    temp.assign(size, 0);
    ld_mat.at(i) = temp;
  }
  for (it3 = inpop.begin(); it3 != inpop.end(); it3++){
    boost::shared_ptr<Haplotype::Node> ptr;
    ptr = it3->start_ptr;
    while(ptr != NULL && ptr->next !=NULL){
      double pos1 = ptr->pos;
      if (sites.find(pos1) != sites.end()){
	boost::shared_ptr<Haplotype::Node> ptr2;
	ptr2 = ptr->next;
	while(ptr2 != NULL){
	  //cout << "here\n";
	  double pos2 = ptr2->pos;
	  if (sites.find(pos2) != sites.end()){
	    int index1 = pos2index[pos1];
	    int index2 = pos2index[pos2];
	    ld_mat.at(index1).at(index2) = ld_mat.at(index1).at(index2)+1;
	    ld_mat.at(index2).at(index1) = ld_mat.at(index2).at(index1)+1;
	    ptr2 = ptr2->next;
	  }
	  else{
	    ptr2 = ptr2->next;
	  }
	}
	ptr = ptr->next;
      }
      else{
	ptr = ptr->next;
      }
    }
  }
  for (it2 = sites.begin(); it2 != sites.end(); it2++){
    for (it6 = it2; it6 != sites.end(); it6++){
      if (it2 != it6){
	double pos1 = *it2;
	double pos2 = *it6;;
	double freq1 = (double)counts[pos1]/(double)N;
	double freq2 = (double)counts[pos2]/(double)N;
	int index1 = pos2index[pos1];
	int index2 = pos2index[pos2];
	double both = ld_mat.at(index1).at(index2);
	double d = ((double)both/(double)N) - freq1*freq2;
	double r2 = (d*d)/ (freq1*freq2*(1-freq1)*(1-freq2));
	outfile << pos1 << " " << pos2 << " " << r2 << "\n";
	ld_mat.at(index1).at(index2) = r2;
	ld_mat.at(index2).at(index1) = r2;
      }
    }
  }
}

double Population::get_bg_rrate(){
  double hspotlen = 0;
  double nhspotlen = 0;
  double total_len = 0;
  double rbg;
  vector<boost::tuple<double, double, double> >::iterator it;
  for (it = hspots.begin(); it != hspots.end(); it++){
    boost::tuple<double, double, double> tmp = *it;
    double len = boost::get<1>(tmp);
    double inten = boost::get<2>(tmp);
    hspotlen +=len;
    total_len += (len*inten);
  }
  nhspotlen = 1-hspotlen;
  total_len += nhspotlen;
  rbg = r/total_len;
  return(rbg);
}


void Population::print_map(const char out[], int len, int scale){
  ofstream outfile(out);
  outfile.precision(15);
  cout.precision(15);
  double bgrate = get_bg_rrate();
  //cout << bgrate<< "\n";
  set<double> sites = get_sites();
  set<double>::iterator it;
  int i =0;
  double mappos = 0;
  double previous;
  for(it = sites.begin(); it!=sites.end(); it++){
    double pos = *it;
    bool found = false;
    if (it == sites.begin()){
      previous = pos;
    }
    else{
      vector<boost::tuple<double, double, double> >::iterator it2;
      double maplength = (pos-previous)*bgrate;
      for (it2 = hspots.begin(); it2 !=hspots.end(); it2++){
	double start = boost::get<0>(*it2);
	double len = boost::get<1>(*it2);
	double intense = boost::get<2>(*it2);
	double end = start+len;
	if (pos > start && previous < end){
	  //cout << "pos:"<<pos << " previous:"<<previous<< " start:"<<start<< " end:" << end << " intense:"<<intense << "\n"; 
	  if (start < previous && end > pos){
	    maplength += (pos-previous)*bgrate*(intense-1);
	  }
	  else if (start > previous && end < pos){
	    maplength += len*bgrate*(intense-1);
	  }
	  else if (start > previous && end > pos){
	    maplength += (pos-start)*(intense-1)*bgrate;
	  }
	  else if (start < previous && end < pos){
	    maplength += (end-previous)*bgrate*(intense-1);
	  }
	  
	}
      }
      mappos +=maplength;
    }
    outfile << "rs" << i << " " <<pos*len << " " << mappos*(double)scale << " A A\n";
    previous = pos;
    i++;
  }
  outfile.close();
}


void Population::print_map(const char out[], int len, int scale, set<double> sites){
  ofstream outfile(out);
  outfile.precision(15);
  cout.precision(15);
  double bgrate = get_bg_rrate();
  set<double>::iterator it;
  int i =0;
  double mappos = 0;
  double previous;
  for(it = sites.begin(); it!=sites.end(); it++){
    double pos = *it;
    bool found = false;
    if (it == sites.begin()){
      previous = pos;
    }
    else{
      vector<boost::tuple<double, double, double> >::iterator it2;
      double maplength = (pos-previous)*bgrate;
      for (it2 = hspots.begin(); it2 !=hspots.end(); it2++){
	double start = boost::get<0>(*it2);
	double len = boost::get<1>(*it2);
	double intense = boost::get<2>(*it2);
	double end = start+len;
	if (pos > start && previous < end){
	  //cout << "pos:"<<pos << " previous:"<<previous<< " start:"<<start<< " end:" << end << " intense:"<<intense << "\n"; 
	  if (start < previous && end > pos){
	    maplength += (pos-previous)*bgrate*(intense-1);
	  }
	  else if (start > previous && end < pos){
	    maplength += len*bgrate*(intense-1);
	  }
	  else if (start > previous && end > pos){
	    maplength += (pos-start)*(intense-1)*bgrate;
	  }
	  else if (start < previous && end < pos){
	    maplength += (end-previous)*bgrate*(intense-1);
	  }
	  
	}
      }
      mappos +=maplength;
    }
    outfile << "rs" << i << " " <<pos*len << " " << mappos*(double)scale << " A A\n";
    previous = pos;
    i++;
  }
  outfile.close();
}
