

/************************************************************************
   Unrooted REConciliation version 1.00 
   (c) Copyright 2005-2006 by Pawel Gorecki
   Written by P.Gorecki.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. 
*************************************************************************/

#ifndef _TREE__
#define _TREE__

#include <iostream>
using namespace std;


class Tree 
{
 public:
  string treename;
  double weight;
  
 Tree(double weight=1.0) : treename(""),weight(weight) {}
  virtual ostream& print(ostream&s)=0; 
  virtual int leaves()=0; 
  virtual int rooted()=0;
  virtual bool eq(Tree *s)=0;
};


#endif
