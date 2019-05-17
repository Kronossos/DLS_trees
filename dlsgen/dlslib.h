#ifndef _DLSLIB_H_
#define _DLSLIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#define hasopt(X) (strchr(opt,X)!=NULL)

#define GENEDUPLICATIONS 0
#define DEEPCOALESCENCE 1
#define DEEPCOALESCENCEPATH 2
#define EXTERNALPATHLENGTH 3
#define DUPLOSS 4
#define LOSS 5
#define TRIPLET 6
#define COPHENETIC 7

// Reading trees
#define MAXBUF 10000 



#define opts_d (strchr(opts,'d'))
#define opts_o (strchr(opts,'o'))
#define opts_D (strchr(opts,'D'))
#define opts_m (strchr(opts,'m'))




typedef unsigned long long costtype;


typedef struct Node
{
  struct Node *a,*b,*p,*s;
  int depth;  
  struct Node *map;
  bool lcaspec;
  Node *scenmap;
  bool scenspec;
  char label; // label==0 - internal node
  char lfid; // for accessing clusters 
  char id;   // position in array of nodes  
  // for ruleabs
  int duplevel;
  int hasscen; // only in absorptionrule
} Node;

typedef struct Scen
{
  long int id;
  int *epi;
  int *mprofile;
  int *rule;
  int mp;
  int me;
  int level;
  int printinfo(int minmes, int minmps, int stsize, int gtsize, Node *gt, const char *opts);
} Scen;


typedef struct  LCA 
{
  Node **t; // node list
  Node ***lcat; 
  Node *root;
  int n; // leaves num
  int tfullsize; // size of t
  Node *a2n[256]; // ascii2pos
} LCA;


void treecluster(Node *s, char *buf, const char *sep="");
Node *parseTree(char *s);
void ppTree(Node *n, int e=1);
LCA* genLCA(Node *s);

#define getlca(a,b) s->lcat[a->id][b->id]
Node* map(LCA *s, Node *g);

long int scenCount(Node *gt, Node *stroot, bool fixedspec);

int treeleafnum(Node *s);
void treecluster(Node *s, char *buf, const char *sep);
Node *lca(Node *a, Node *b, LCA* l);
void insertNodes(Node *s,Node **t,int *lfp, int *itp, int depth);

typedef struct Cand
{
  Node *cand;
  Node *smap;
  int tp;
} Cand;

void ruleabsorption2(Scen *lca,Node *gt, Node *st,Node **gtnodes, LCA *s, int stsize, 
  int gtsize, Scen** scentab,int *minmetab,int *minmptab,const char *opts,int minmes, int minmps);


Scen* genscen(Scen *scen, long id, Node *gt, Node *stroot, LCA  *s, 
  bool fixedspec, int stsize, int gtsize, bool maponly=0);


void ppdls(Node *n, Node *pp=0, int e=1);
int ruletransf(Scen *a, Scen *b, int gtsize);

//-----------------------
int mpscore(int *pprofile, int tsize);
int* calcmprofile(int *epi, LCA *s, int tsize, int *pprofile);

int* calcepi(Node *gt, LCA *s, int *epi=NULL);
int* calcrule(Node *gt, int gsize, int *rule=NULL);
int mescore(int *epi, int tsize);
#endif