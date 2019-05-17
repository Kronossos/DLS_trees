
/************************************************************************
 Unrooted REConciliation version 1.00
 (c) Copyright 2005-2006 by Pawel Gorecki
 Written by P.Gorecki.
 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/

#ifndef _UNROOTED__
#define _UNROOTED__
#include <iostream>
#include <map>
#include <set>
#include <list>
#include <cstdlib>

using namespace std;

#include "rtree.h"

#define C_MAP 1
#define C_SC 2
#define C_COST 4

#define M_MARK 1
#define M_OPT 2
#define M_START 4
#define M_OPTM 8
#define M_MINC 16
#define M_ROOTING 32
#define M_PROCESSED 512

#define GSFULL 0
#define GSPOS 1
#define GSAFTER 2
#define GSBEFORE 3

class UNode;
class ULeaf;
typedef set<UNode*> nodset;

extern int detailed_costs;
int lossprim(RNode *s, RNode *s1, RNode *s2);
void dlcostdet(RNode *s, RNode *s1, RNode *s2);
#define dupprim(s,s1,s2) (( (s==s1) || (s==s2))?1:0)

class UNode {
protected:
	UNode *pn;
	RNode *Mn;
	DlCost scn;
	DlCost costn;
	int computed;
	int ismarked;
	double branchlen;
	int smpcluster; // used for finding strong clusters like (a,(a,(a,a)))
			// edges inside such a strong cluster will be marked as a-species
			// default value -1
public:

	friend class UNode3;

	UNode(UNode *p_ = NULL) :
		pn(p_), Mn(NULL), computed(0), ismarked(0),branchlen(-1),
		smpcluster(-1)
	{
		if (pn) branchlen=p_->branchlen;
	}
	virtual ~UNode() {
	}
	void reset() {
		computed = 0;
		Mn = NULL;
		ismarked = 0;
	}

	double len() { return branchlen; }
	void setbranchlen(double len) { branchlen = len; }
	virtual int leaves()=0;
	void mark(int m = 1) {
		ismarked |= m;
	}
	int marked() {
		return ismarked;
	}
	virtual ULeaf* minleaf()=0;
	virtual int leaf()=0;
	virtual void clear()=0;
	virtual UNode *p() {
		return pn;
	}
	virtual void balance(int) {}

	void p(UNode *p_) {
		// set branch len from p
		pn = p_;
		if (pn->branchlen>=0)
		  branchlen=pn->branchlen;
		else
		  pn->branchlen=branchlen;
	}
	virtual ostream& printrootings(ostream&s, int) { return printrooting(s,0); }


	DlCost &cost(SpeciesTree *st) {
		if (!pn)
			return costn;
		if (!(computed & C_COST)) {
			RNode *s = st->lca(M(st), pn->M(st));
			costn.loss = sc(st).loss + pn->sc(st).loss + lossprim(s, M(st),
					pn->M(st));
			costn.dup = sc(st).dup + pn->sc(st).dup
					+dupprim(s,M(st),pn->M(st));
			computed |= C_COST;
		}
		return costn;
	}
	void costdet(SpeciesTree *st) {
		if (!pn)
			return; // nothing to compute (a leaf)
		RNode *s = st->lca(M(st), pn->M(st));
		dlcostdet(s, M(st), pn->M(st));
		costdetsubtree(st);
		pn->costdetsubtree(st);
	}
	virtual void costdetsubtree(SpeciesTree *st) {
	}
	virtual ostream& ppsmprooted(ostream&s, int fromroot=0)=0;
	virtual ostream& ppsmprootedext(ostream&s, RNode *sroot, int fromroot=0)=0;
	virtual RNode *smprooted()=0;
	virtual RNode *M(SpeciesTree *st)=0;
	virtual ostream& printbranchlen(ostream &s, int fromroot=0)
	{ if (branchlen>=0)
		s << ":" << ((fromroot)?0.5*branchlen:branchlen);
	return s;
	}
	// UNode
	virtual ostream& printrooting(ostream&s, int fromroot) {
	  if (pn) {
	    s << "(";
	    ppsmprooted(s,1) << ",";
	    pn->ppsmprooted(s,1) << ")";
	    return s << endl;
	  }
	  return ppsmprooted(s);
	}

	virtual ostream& printrootingext(ostream&s, RNode *sroot, int fromroot) {
		if (pn) {
		  if (sroot!=Mn && sroot!=pn->Mn) 
		    return printrooting(s,1); 
		  ppsmprootedext(s,sroot,1);
		  pn->ppsmprootedext(s,sroot,1);
		  return s;
		}		
		ppsmprooted(s) << endl;
		return s;
	}

	virtual RNode* rooted() {
		if (pn)
		  return new RInt(smprooted(), pn->smprooted());
		return smprooted();
	}

	virtual nodset* insert(nodset *n)=0;
	virtual DlCost &sc(SpeciesTree *st) {
		return scn;
	}
	virtual ostream& pf(ostream& s, DlCost &c, SpeciesTree *st) {
		if (pn) {
			s << "(";
			smppf(s, c, st) << ",";
			return pn->smppf(s, c, st) << ")";
		}
		return smppf(s, c, st);
	}
	virtual ostream& smppf(ostream &s, DlCost&, SpeciesTree*)=0;
	void pcosts(ostream &s,DlCost &c,SpeciesTree *st);
	virtual UNode* subtreecost(SpeciesTree *st)=0;
	virtual UNode* mincost(SpeciesTree *st)=0;
	virtual int printgraphviz(ostream &,int&)=0;

};

class ULeaf: public UNode {
 protected:
  int specid;
  char *name;
 public:
  
 ULeaf(int spec, char *fullname, UNode *p_=NULL) : UNode(p_), specid(spec),name(fullname) 
  { mmi("ULf"); }
  virtual ~ULeaf() {	 }
  virtual int leaf() { return 1; }
  int species() { return specid; }
  virtual void clear() {
    reset();
  }
  virtual ULeaf* minleaf() { return this; }
  
  virtual string lfspecname()
  {
    //cout << "!!" << name;
    if (name) return name;
    return getspecname(specid,name);
  }
  
  virtual ostream& ppsmprooted(ostream&s, int onlylabs=0)
  {
    s << lfspecname();
    return printbranchlen(s);
  }
  
  virtual ostream& ppsmprootedext(ostream&s, RNode *sroot,  int onlylabs=0)
  {
    s << lfspecname();
    return printbranchlen(s) << endl;
  }
  
  virtual RNode *smprooted() {
    return new RLeaf(specid,name);
  }
  virtual nodset* insert(nodset *n) {
    n->insert(this);
    return n;
  }
  virtual int leaves() { return 1; }
  
  virtual int printgraphviz(ostream & os,int &c) {
    os << c << "[shape=plaintext,label=\"" <<  lfspecname() << "\",height=0,width=0]" << ";" << endl;
    return c++;
  }
  
  virtual RNode *M(SpeciesTree *st) {
    if (!(computed & C_MAP))
      {
	Mn=st->getLeaf(specid);
	if (!Mn) {
	  cerr << "Mapping of " << specnames[specid] << " not found in the species tree." << endl;
	  exit(-1);
	}
	computed|=C_MAP;
      }
    return Mn;
  }
  
  
  virtual ostream& smppf(ostream& s,DlCost &c,SpeciesTree *st) {
		lfspecname();
		pcosts(s,c,st);
		return s;
  }
  
  virtual UNode* subtreecost(SpeciesTree *st) {
    cost(st);
    return this;
  }
  virtual UNode* mincost(SpeciesTree *st) {
    cost(st);
    if (!pn)
      return this;
    if (pn->leaf())
      return this;
    return pn->subtreecost(st);
  }
};

class UNode3: public UNode {
protected:
	UNode3 *ln;
	UNode3 *rn;
	void connect(UNode3 *a, UNode3 *b);
public:
	UNode3(UNode *p_ = NULL) :
		UNode(p_) {
	}
	~UNode3() {
	}
	virtual int leaf() {
		return 0;
	}
	UNode3 *l() {
		return ln;
	}
	UNode3 *r() {
		return rn;
	}
	virtual void clear() {
		reset();
		ln->reset();
		rn->reset();
		ln->p()->clear();
		rn->p()->clear();
	}
	void l(UNode3 *l_) {
		ln = l_;
	}
	void r(UNode3 *r_) {
		rn = r_;
	}


	virtual void balance(int lastside)
	{
		if (ln->pn->leaf() && rn->pn->leaf()) return;
		if (ln->pn->leaf())
		{
			//	    cout << "BALANCING " << endl;
			//	    ppsmprooted(cout) << endl;
			if (lastside==-1)
			{
				// cout << "BALANCING - l" << endl;
				UNode *pl = ln->pn;
				UNode *pr = rn->pn;
				ln->pn=pr;
				rn->pn=pl;
				pl->pn=rn;
				pr->pn=ln;
				ln->pn->balance(1);
			}
			else
				rn->pn->balance(-1);
			//	    cout << "RESULT: " << endl;
			//	    ppsmprooted(cout) << endl;
			return;
		}
		if (rn->p()->leaf())
		{
			//	    cout << "BALANCING" << endl;
			//	    ppsmprooted(cout) << endl;
			if (lastside==1)
			{
				//		cout << "BALANCING - r" << endl;
				UNode *pl = ln->pn;
				UNode *pr = rn->pn;
				ln->pn=pr;
				rn->pn=pl;
				pl->pn=rn;
				pr->pn=ln;
				rn->pn->balance(-1);
			}
			else
				ln->pn->balance(1);

			//	    cout << "RESULT: " << endl;
			//	    ppsmprooted(cout) << endl;
			return;
		}
		rn->pn->balance(0);
		ln->pn->balance(0);
	}

	virtual int printgraphviz(ostream & os,int &c) {
		int c1=rn->p()->printgraphviz(os,c);
		int c2=ln->p()->printgraphviz(os,c);
		os << c << "[shape=point,width=0.01];" << endl;
		os << c << " -- " << c1 << "[len=1];" << endl;
		os << c << " -- " << c2 << "[len=1];" << endl;
		return c++;

	}

	virtual int leaves() { return ln->p()->leaves()+rn->p()->leaves(); }

	virtual RNode *M(SpeciesTree *st) {
		if (!(computed & C_MAP)) {
			Mn = st->lca(ln->p()->M(st), rn->p()->M(st));
			computed |= C_MAP;
		}
		return Mn;
	}
	virtual void costdetsubtree(SpeciesTree *st) {
		dlcostdet(M(st), ln->p()->M(st), rn->p()->M(st));
		ln->p()->costdetsubtree(st);
		rn->p()->costdetsubtree(st);
	}
	virtual DlCost& sc(SpeciesTree *st) {
		if (!(computed & C_SC)) {
			scn.loss = ln->p()->sc(st).loss + rn->p()->sc(st).loss + lossprim(
					M(st), ln->p()->M(st), rn->p()->M(st));
			scn.dup = ln->p()->sc(st).dup + rn->p()->sc(st).dup
					+dupprim(M(st),ln->p()->M(st),rn->p()->M(st));
			computed |= C_SC;
		}
		return scn;
	}
	// Unode3
	virtual ostream& printrootings(ostream&s, int nonrootnonroot);
	virtual ostream& ppsmprooted(ostream&s, int fromroot=0) { 
	  s << "(";
	  ln->p()->ppsmprooted(s) << ",";
	  rn->p()->ppsmprooted(s) << ")";
	  return printbranchlen(s,fromroot);
	}
	
	virtual ostream& ppsmprootedext(ostream&s, RNode *sroot, int fromroot=0) { 	  
	  if (Mn==sroot && ( ln->p()->Mn==sroot || rn->p()->Mn==sroot)) // dupl.
	    {
	      ln->p()->ppsmprootedext(s,sroot); 
	      rn->p()->ppsmprootedext(s,sroot);
	      return s;
	    }
	  return ppsmprooted(s,fromroot) << endl;
	  
	}

	virtual nodset* insert(nodset *n) {
		n->insert(this);
		n->insert(ln);
		n->insert(rn);
		ln->p()->insert(n);
		rn->p()->insert(n);
		return n;
	}
	virtual RNode* smprooted() {
		return new RInt(ln->p()->smprooted(), rn->p()->smprooted());
	}
	virtual ostream& smppf(ostream& s, DlCost &c, SpeciesTree *st) {
		s << "( (";
		ln->p()->smppf(s, c, st) << ") ";
		ln->pcosts(s, c, st);
		s << ", ( ";
		rn->p()->smppf(s, c, st) << " ) ";
		rn->pcosts(s, c, st);
		s << " )";
		pcosts(s, c, st);
		return s;
	}

	virtual ULeaf* minleaf() {
	  ULeaf *r = rn->p()->minleaf();
	  ULeaf *l = ln->p()->minleaf();
	  if (r->species()<l->species()) return r;
	  return l;
	}
	
	
	virtual UNode* subtreecost(SpeciesTree *st) {
	  UNode *res = ln->p()->subtreecost(st);
	  UNode *res1 = rn->p()->subtreecost(st);
	  if (res->cost(st).mut() > res1->cost(st).mut())
	    res = res1;
	  if (res->cost(st).mut() > cost(st).mut())
	    return this;
	  return res;
	}
	
	virtual UNode* mincost(SpeciesTree *st) {
	  UNode *res = pn->subtreecost(st);
	  UNode *res1 = ln->p()->subtreecost(st);
	  if (res->cost(st).mut() > res1->cost(st).mut())
	    res = res1;
	  res1 = rn->p()->subtreecost(st);
	  if (res->cost(st).mut() > res1->cost(st).mut())
	    return res1;
	  return res;
	}
};

class UTree;
class iterator_utree {
protected:
	nodset::iterator nit;
	nodset* nodes;
	UNode *c;
	int flag;
public:
	iterator_utree(UTree *t, int flag_ = F_ALL);
	UNode *operator()();
};

class UTree : public Tree {
protected:

	UNode *toUNodes(RNode *t);
	UNode3* connect(UNode3 *a, UNode3 *b, UNode3 *c, UNode *u1, UNode *u2, double branchlen);
	virtual UNode *createLeaf(int spec, char *fullname) { return new ULeaf(spec,fullname); }
	virtual UNode *createNode3(UNode *u1, UNode *u2, double branchlen=-1.0) 
	{
	  return connect(new UNode3(), new UNode3(), new UNode3(), u1, u2, branchlen);
	}
	UNode *parseNode(char *s, int &p, int fromroot,int num,int extractspecies);
	void initrand(int len,double pint, double dec, vector<string> &t, int splen);



	int lsize;
public:
	UNode *start;
	UTree(char *t, int num,string tn="", double wg=1.0, int extractspecies=0);
	UTree() { mmi("UT-empty"); start=NULL; lsize=-1; }
	UTree(int len,double pint, double dec, SpeciesTree *sp);
	UTree(SpeciesTree *sp);
	UTree(UTree*, set<int> &tr);
	UTree(int len,double pint, double dec, int numlv, int uniquelv, vector<string> &t);

	virtual ~UTree() {
	}

	int leaves();

	friend class iterator_utree;
	virtual ostream& printrootings(ostream&s);
	nodset* nodes() {
		return start->insert(start->p()->insert(new nodset));
	}
	UNode *findoptimaledge(SpeciesTree *st);
	void clear() {
		start->clear();
		if (start->p())
			start->p()->clear();
	}
	void pf(ostream &s, SpeciesTree *st) {
		s << "[";
		start->pf(s, mincost(st)->cost(st), st);
		s << "]" << endl;
	}
	UNode* mincost(SpeciesTree *st) {
		return start->mincost(st);
	}
	//UNode *genRand(double pint, double dec, char **t, int s);
	UNode *genRand(double pint, double dec, vector<string> &t, int s);

	virtual ostream& print(ostream&s) {
		return cout << *start->rooted();
	}
	;

	int rooted() { return 0; }
	bool eq(Tree *s) { if (s->rooted()) return false;
	  cout << "unrooted comparison - undefined";
	  exit(-1);
	}

	// Return a species tree s.t. the root has a leaf child with min. label
	// Only for unique labellings
	virtual SpeciesTree* normalize() {
		if (!start->p()) return new SpeciesTree(start->rooted());
		ULeaf *l1 = start->minleaf();
		ULeaf *mm = start->p()->minleaf();
		if (l1->species()<mm->species()) mm=l1;
		return new SpeciesTree(mm->rooted());
	};


	ostream& printgraphviz(ostream &os)
	{
		int c=0;
		os << "graph G" << endl;
		os << "{  graph [fontsize=12];" << endl;
		os << "  edge  [fontsize=12];"  << endl;
		os << "  node  [shape=point];" << endl;
		os << "  ranksep = 0.1;" << endl;
		os << "  nodesep = .25;" << endl;
		os << "  edge [style=\"setlinewidth(2)\"];" << endl;
		int s = start->printgraphviz(os,c);
		if (start->p())
		{
			int p = start->p()->printgraphviz(os,c);
			os << s << " -- " << p << "[len=0.5];" << endl;
		}

		return os << "}" << endl;
	}
	UTree *trunc(set<int> &tr)
	{
		UTree *ut = new UTree(this,tr);
		if (ut->start==NULL) return NULL;
		return ut;
	}
	void center();

};

#endif
