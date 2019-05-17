/*
 * tools.h
 *
 *  Created on: 18-09-2010
 *      Author: gorecki
 */

#ifndef TOOLS_H_
#define TOOLS_H_

#include <set>
using namespace std;

#include <set>
#include <vector>
#include <iomanip>
using namespace std;

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "rtree.h"
#include "urtree.h"




void readgtree(char *fn, vector<UTree*> &gtvec);
void readstree(char *fn, vector<SpeciesTree*> &stvec);

void readsmp(istream &c, vector<SpeciesTree*> &vec) ;
//void readgt(istream &c, int treenames, int weights, vector<UTree*> &);
//void readst(istream &c, int treenames, int weights, vector<SpeciesTree*> &);
void readt(istream &c, int treenames, int weights, int rooted, vector<Tree*> &);

void readfields(istream &c, int rooted, int field,vector<string> &prefix, vector<string> &suffix,vector<Tree*> &vec);

void nnist(vector<SpeciesTree*> &stvec);
void nnigt(vector<UTree*> &gtvec);
void nnit(vector<Tree*> &vec);

void tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters);

string printpstricks(UTree *t,UNode *n);

vector<string> initgenrand(char *optarg);

void setalordering();

void sortmode(char *optarg);

extern vector<string> specnames;
extern vector<int> specorder;
extern map<string,int> specnames2id;
extern int rsort; // Important in nni - used to generated unique r|u species trees
extern int gcenter; // Useful for drawing gene trees
extern int costfromroot;
extern int printleaves;
extern int stprogress; // Print progress for species trees
extern int ppgraphviz; // Print graphviz output

void readtreesinterleaved(char *fn,vector<UTree*> &gtvec, vector<SpeciesTree*> &stvec);


#endif /* TOOLS_H_ */
