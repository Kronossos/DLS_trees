
/************************************************************************
 Unrooted REConciliation version 1.00
 (c) Copyright 2005-2006 by Pawel Gorecki
 Written by P.Gorecki.
 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/

#include <ctype.h>
#include <iostream>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "tree.h"
#include "rtree.h"

double weight_loss = 1.0;
double weight_dup = 1.0;
extern int rsort;

int getspecies(char *s,int len);

int usespecnames=0;

string getspecname(int specid,char *name)
{
  if ((!usespecnames) && name)
    return string(name);
  return specnames[specid];
}

char* getTok(char *s,int &p, int num)
{
#define inctok p++
  while (isspace(s[p])) inctok;
  char *cur = s+p;
  if (isalnum(s[p]) || (s[p]=='_'))
    {
      while (isalnum(s[p]) || (s[p]=='_') || (s[p]=='.') || (s[p]=='-')) inctok;
      return cur;
    }
  if ((s[p]=='(')  || (s[p]==')') || (s[p]==',') || s[p]==':')
    {
      inctok;
      return cur;
    }
  if (num>=0) cerr << "Tree nr " << num << ". ";
  cerr << "Parse error. Expected expected LABEL, '(', ')' or ','. Found: '" << s[p]
       << "'. Current string: <<" << cur << ">>"<< endl;
  exit(-1);
}

iterator_tree::iterator_tree(RTree *tr, int flag_) {
	t = tr;
	c = tr->root();
	flag = flag_;
}


RNode *RTree::parseNode(char *s,int &p, int num)
{
	char *cur = getTok(s,p,num);
	if (cur[0]=='(')
	{
		RNode *a = parseNode(s,p,num);
		getTok(s,p,num);
		RNode *b = parseNode(s,p,num);
		getTok(s,p,num);
		RNode *r;
		if (rsort && (a->minlab() >= b->minlab())) r=createInt(b,a);
		else r= createInt(a,b);
		return r;
	}
	return createLeaf(getspecies(cur,s+p-cur));
}

RTree::RTree(char *fs, int num, string tn, double wg)
{
	int px=0;
	rootn=parseNode(fs,px,num);
	rootn->depth(0);
	treename=tn;
	weight=wg;
}

RNode *iterator_tree::operator()() {
	RNode *res = step();
	while (res) {
		if (flag & F_ALL)
			break;
		if ((flag & F_INTERNAL) && (!res->leaf()))
			break;
		if ((flag & F_LEAVES) && (res->leaf()))
			break;
		res = step();
	}
	return res;
}

RNode *iterator_tree::step() {
	RNode *res = c;
	if (!c)
		return c; // finished
	if (!c->leaf())
		c = ((RInt*) c)->l();
	else {
		while (c) {
			RNode *prev = c;
			c = c->p();
			if (!c)
				break; // last
			if ((((RInt*) c)->l()) == prev) {
				c = ((RInt*) c)->r();
				break;
			}
		}
	}
	return res;
}

RNode *SpeciesTree::lca(RNode *a, RNode *b) {
	if (b->isParentOf(a))
		return b;
	while (a) {
		if (a->isParentOf(b))
			return a;
		a = a->p();
	}
	return NULL;
}

RNode *RNode::isParentOf(RNode *c) {
	while (c) {
		if (c == this)
			return c;
		c = c->p();
	}
	return 0;
}

char* xstrndup(const char *s, int len) {
	if (len == 0)
		return strdup(s);
	char *b = new char[len + 1];
	strncpy(b, s, len);
	b[len] = 0;
	return b;
}

int getspecies(char *s,int len)
{

	char old;
	int num;
	if (len)
	{
		old=s[len];
		s[len]=0;
	}
	//	cout << "GETSP:" << s << ":" << len << endl;
	if (specnames2id.count(s)==0)
	{
		specnames2id[s]=num=specnames.size();
		specorder.push_back(num);
		specnames.push_back(s);
		//cout << "new species: #" << s << endl;
	}
	else num=specnames2id[s];
	if (len) s[len]=old;

	return num;
}
