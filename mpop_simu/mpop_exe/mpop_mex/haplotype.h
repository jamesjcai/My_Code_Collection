#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

class Haplotype{
 public:
  Haplotype();
  ~Haplotype();
  Haplotype(const Haplotype&) throw();
  const Haplotype &operator=(const Haplotype&);
  void del();
  void print();
  void add_node(double);
  void add_node_at_end(double);
  void recombine(Haplotype*, double);
  bool has_pos(double);
  bool empty();
  struct Node{
    double pos;
    boost::shared_ptr<Node> next;
  };
  boost::shared_ptr<Node> get_start_ptr();
  boost::shared_ptr<Node> get_end_ptr();
  boost::shared_ptr<Node> start_ptr, end_ptr;
 private:
  boost::shared_ptr<Node> node_before(double);
  boost::shared_ptr<Node> node_after(double);
  double selected_site;
  bool has_site;
};	
