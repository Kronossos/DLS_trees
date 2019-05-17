
/************************************************************************
   Unrooted REConciliation version 1.00 
   (c) Copyright 2005-2006 by Pawel Gorecki
   Written by P.Gorecki.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. 
*************************************************************************/

#ifndef _ROOTED__
#define _ROOTED__

#define mmi(a)

#define T_LABS 1
#define T_TREEDIV 2
#define T_VAL 4
#define T_ALL 7
#define T_START 8


#include <iostream>
#include <map>
#include <vector>
using namespace std;

#include "tree.h"

char* getTok(char *s,int &p, int num);
int getspecies(char *s,int len=0);
string getspecname(int specid,char *name);
char* xstrndup(const char *s,int len);

extern double weight_loss;
extern double weight_dup;
extern int rsort;
extern vector<string> specnames;
extern vector<int> specorder;
extern map<string,int> specnames2id;

typedef struct DlCost 
{
    int dup;
    int loss;
    DlCost(int dl=0,int ls=0) : dup(dl), loss(ls) {}
    friend ostream& operator<<(ostream&s, DlCost p)  
	{ return s << "(" << p.dup << "," << p.loss << ")"; } 
    friend DlCost operator+(DlCost s1,DlCost s)  
	{ return DlCost(s1.dup+s.dup,s1.loss+s.loss); }
    double mut() { return weight_dup*dup+weight_loss*loss; }
} DlCost;

class RInt;

class RNode
{
 protected:
    RInt *pn;
    int depthn;
    DlCost dc;
 public:
    RNode() { pn=NULL; }
    virtual ~RNode() {}
    virtual int leaf() { return 0; }

    virtual bool eq(RNode *)=0;
    virtual ostream& print(ostream&s)  { return s; }
    int depth() { return depthn; }
    virtual void depth(int d) { depthn=d; }
    friend ostream& operator<<(ostream&s, RNode &p)  { return p.print(s); }  
    virtual RInt *p() { return pn; }
    virtual void p(RInt *p) { pn=p; }
    RNode *isParentOf(RNode *c);
    DlCost &costdet() { return dc; }
    virtual ostream& printsorted(ostream&s, int type=0)=0;
    virtual int minlab()=0;
    virtual void showcostdet(ostream&s) { s << dc << " : "; print(s); s << endl; }
    virtual DlCost subtreecost()=0; 
    virtual void nni(ostream&s,RNode *r, int unni)=0;
    virtual void pfcostdet(ostream&s) {  
	s << " dup(" << dc.dup << ")" << " loss(" << dc.loss <<")"; };
};

class RInt : public RNode
{
 protected:
    RNode *ln;
    RNode *rn;

 public:
    RInt(RNode *_l,RNode *_r) : RNode(), ln(_l), rn(_r) 
    { ln->p(this); rn->p(this); }
    virtual void depth(int d) { depthn=d; ln->depth(d+1); rn->depth(d+1); }
    ~RInt() {}
    RNode *r() { return rn; }    
    RNode *l() { return ln; }
    virtual ostream& print(ostream&s)  { return s << "(" << *ln << "," << *rn << ")"; }    
    virtual void showcostdet(ostream&s) { 
    	RNode::showcostdet(s);
    	ln->showcostdet(s);
    	rn->showcostdet(s);
    }
    virtual bool eq(RNode *n) { 
      if ((n->leaf()) || (minlab()!=n->minlab())) return false; 
      RInt *rr=(RInt*)n;
      if (l()->eq(rr->l()) && r()->eq(rr->r())) return true;
      return l()->eq(rr->r()) && r()->eq(rr->l());
    } 
    virtual DlCost subtreecost() { return ln->subtreecost() + dc + rn->subtreecost(); }
    virtual void pfcostdet(ostream&s) { 
    	s << "(";
    	ln->pfcostdet(s);
    	s << ",";
    	rn->pfcostdet(s);
    	s << ")";
    	RNode::pfcostdet(s);
    } 
    virtual int minlab() {
      int r = rn->minlab();
      int l = ln->minlab();
      if (r<l) return r;
      return l;
    }
    virtual ostream& printsorted(ostream&s, int type=0)  {
    	int fromroot = type & T_TREEDIV;
    	type&=~T_TREEDIV;
    	if (!type) s << "(";
    	if (ln->minlab()<rn->minlab())
    	{
    		ln->printsorted(s,type);
    		if (!type) s << ",";
    		if (fromroot) s << "-";
    		rn->printsorted(s,type);
    	}
    	else
    	{
	  rn->printsorted(s,type);
	  if (!type) s << ",";
	  if (fromroot) s << "-";
	  ln->printsorted(s,type);
    	}
    	if (!type) s << ")";
    	return s;
    }

    virtual void nni(ostream&s,RNode *r, int unni) {

    	ln->nni(s,r,unni);
    	rn->nni(s,r,unni);

    	int i;

    	if (unni && (this==r)) return;
    	for (i=0;i<2;i++)
    	{
    		RNode *l;
    		if (!rn->leaf())
    		{
    			RInt *n = (RInt*)rn;
    			l = ln;

    			ln=n->ln;
    			n->ln=l;
    			r->printsorted(s) << endl;

    			l=ln;
    			ln=n->rn;
    			n->rn=l;
    			r->printsorted(s) << endl;

    			// reconstruct
    			l=n->ln;
    			n->ln=ln;
    			ln=l;
    		}
    		l=ln;
    		ln=rn;
    		rn=l;
    	}
    }
};


class RLeaf : public RNode
{
 protected:
  int specid;
  char *name;
 public:
 RLeaf(int spec,char *fullname) : RNode(), specid(spec),name(fullname) 
  { }
  ~RLeaf() {}
  virtual int leaf() { return 1; }

    int species() { return specid; }
    virtual ostream& print(ostream&s)  { 
      if (name) return s << name;
      return s << (specnames[specid]); }
    virtual DlCost subtreecost() { return dc; }
    virtual ostream& printsorted(ostream&s,int type=0) { return s << specnames[specid]; }
    virtual void pfcostdet(ostream&s) { 
      s << specnames[specid] ;
      RNode::pfcostdet(s);
    }
    virtual void nni(ostream&s,RNode *r,int unni) {}
    int minlab() { return specorder[specid]; }

    virtual bool eq(RNode *r) { 
      return r->leaf() && (minlab()==r->minlab());
    }
};

class RTree;

#define F_INTERNAL 1
#define F_LEAVES 2
#define F_ALL 4

class iterator_tree 
{
    RTree *t;
    RNode *c;
    int flag;
    RNode *step();
 public:
    iterator_tree(RTree *tr, int flag_=F_ALL);
    RNode *operator()();
};

class RTree : public Tree
{
 protected:
    RNode *rootn;
    virtual RNode *parseNode(char *s, int &p, int num);
    virtual RNode *createLeaf(int species) { return new RLeaf(species,0); }
    virtual RNode *createInt(RNode *a, RNode *b) { return new RInt(a,b); } 
 public:
    RTree(RNode *_root=NULL,double weight=1.0) : Tree(weight) { rootn->depth(0); }
    RTree(char *fromstr, int num, string tn, double wg);

    virtual ~RTree() {}
    void str2tree(char *s) { int p=0; rootn=parseNode(s,p,-1); }
    RNode *root() { return rootn; } 
    virtual ostream& print(ostream&s)  { return s << *rootn; }
    friend ostream& operator<<(ostream&s, RTree &p)  { return p.print(s); }  
    virtual int rooted() { return 1; }
};



typedef map<int,RNode*> species2leaves;


class SpeciesTree : public RTree
{
 protected:
	species2leaves lmap;
	void takeLeaves(RNode *r) {
		if (r->leaf()) lmap[((RLeaf*)r)->species()]=r;
		else { takeLeaves(((RInt*)r)->l()); takeLeaves(((RInt*)r)->r()); }
	}
 public:

 SpeciesTree(char *s, int num,string tn="",double weight=1) : RTree(s,num,tn,weight) { takeLeaves(rootn); }
    SpeciesTree(RNode *_root) : RTree(_root) { mmi("SpTr-node"); takeLeaves(rootn);  rootn->depth(0); }


    virtual ~SpeciesTree() {} 
    RLeaf *getLeaf(int s) {  return (RLeaf*)lmap[s]; }
    int lsize() { return lmap.size(); }
    RNode *lca(RNode *a, RNode *b);    
    int leaves() { return lmap.size(); }
    void showcostdet(ostream&s) { rootn->showcostdet(s); } 
    DlCost totalcost() { return rootn->subtreecost(); } 
    void pfcostdet(ostream&s) { cout << "[ "; rootn->pfcostdet(s); cout << "]" << endl; }
    virtual void nni(ostream&s,int unni=0)  { rootn->nni(s,rootn,unni); }


    virtual ostream& printsorted(ostream&s, int type=0) {
      if (type) type|=T_TREEDIV;
      return rootn->printsorted(s,type);
    }

    bool eq(Tree *s) { 
      if (!s->rooted()) return false;
      return rootn->eq(((RTree*)s)->root()); 
    }

};

#endif
