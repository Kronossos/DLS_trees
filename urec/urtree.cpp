/************************************************************************
 Unrooted REConciliation version 1.00
 (c) Copyright 2005-2006 by Pawel Gorecki
 Written by P.Gorecki.
 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/

#include <ctype.h>
#include <iostream>
#include <cstring>
#include <stdio.h>
using namespace std;

#include "urtree.h"


int gsid=GSFULL;
char *gsdelim=NULL;
int gspos=0;

iterator_utree::iterator_utree(UTree *t, int flag_) :
					flag(flag_) {
	nodes = t->nodes();
	nit = nodes->begin();
}

UNode *iterator_utree::operator()() {
	UNode *un;
	while (1) {
		if (nit == nodes->end())
			return NULL;
		un = *nit;
		if (flag & F_ALL)
			break;
		if ((flag & F_INTERNAL) && !un->leaf())
			break;
		if ((flag & F_LEAVES) && un->leaf())
			break;
		nit++;
	}
	nit++;
	return un;
}

ostream& UNode3::printrootings(ostream&s, int nonroot) {
	if (!nonroot) {
		pn->printrootings(s, 1);
		rn->p()->printrootings(s, 1);
		ln->p()->printrootings(s, 1);
		return s;
	}
	printrooting(s, nonroot); // single rooting
	ln->p()->printrootings(s, nonroot);
	rn->p()->printrootings(s, nonroot);
	return s;
}

ostream& UTree::printrootings(ostream&s) {
  UNode *u=start;
  if (start->leaf()) u=start->p();
  return u->printrootings(s, 0);
}

double parseBranchLen(char *s,int &p, int num)
{
	double branchlen=-10;
	getTok(s,p, num); // eat :
	char *cur = getTok(s,p,num);
	if (sscanf(cur,"%lf",&branchlen)!=1)
	{
		cerr << "Parse error: number expected after :";
		exit(-10);
	}
	return branchlen;

}

UNode *UTree::parseNode(char *s, int &p, int fromroot, int num, int extractspecies) {
	char *cur = getTok(s, p, num);
	double branchlen=-1;
	if (cur[0] == '(') {
	  UNode *a = parseNode(s, p, 0, num, extractspecies );
		getTok(s, p,num);
		UNode *b = parseNode(s, p, 0, num, extractspecies);
		if (fromroot) {
			if (s[p] == ',') {
				getTok(s,p,num);
				a = createNode3(a, b);
				b = parseNode(s, p, 0, num, extractspecies);
				getTok(s, p, num);  // branch length from the root is ignored
			}
			// join a<->b
			a->p(b);
			b->p(a);
			return start = b; // branch length from the root is ignored
		}
		getTok(s,p, num); // eat )
		if (s[p]==':') branchlen=parseBranchLen(s,p,num);
		return createNode3(a, b, branchlen);
	}

	//start = createLeaf(cur, s + p - cur);
	int spid=-1;
	char *fullname=NULL;
	int len=s+p-cur;

	if (extractspecies)
	  {
	    fullname=strndup(cur,len);
	    //cout << "###" << fullname << "###" << extractspecies << " " << gsid << " " << gspos << endl;
	    switch (gsid) { 
	    case GSPOS:
	      if (gspos>0) spid=getspecies(cur,gspos);
	      else spid=getspecies(cur+len+gspos,gspos);
	      break;
	    case GSAFTER:
	      {
		char *st=strstr(fullname,gsdelim);
		if (!st) 
		  {
		    cerr << "No species name found in " << fullname << endl;
		    exit(-1);
		  }
		st=st+strlen(gsdelim);
		//cout << "st=" << st << endl;
		char *ste=strstr(st,gsdelim);
		if (!ste)		
		  spid=getspecies(st,strlen(st));
		else
		  spid=getspecies(st,ste-st);
	      }
	    case GSBEFORE:
	      {
		char *st=strstr(fullname,gsdelim);
		if (!st) 
		  {
		    cerr << "No species name found in " << fullname << endl;
		    exit(-1);
		  }
		spid=getspecies(cur,st-fullname);
	      }
	      break;
	    }	    

	  }
	
	if (spid<0) spid=getspecies(cur,len);
	start=createLeaf(spid,fullname);

	if (s[p]==':')
	  start->setbranchlen(parseBranchLen(s,p, num));
	return start;
}

UNode3* UTree::connect(UNode3 *a, UNode3 *b, UNode3 *c, UNode *u1, UNode *u2, double branchlen) {
	a->l(b);
	b->l(c);
	c->l(a);
	a->r(c);
	b->r(a);
	c->r(b);

	u1->p(a);
	u2->p(b);
	a->p(u1);
	b->p(u2);
	c->setbranchlen(branchlen);
	return c;
}

UNode *UTree::findoptimaledge(SpeciesTree *st) {
#define shw(k) 
	UNode *cur = start;
	shw("start");
	cur->mark(M_START | M_MARK);
	if (!cur->p())
		return cur;
	if (cur->leaf() && cur->p()->leaf())
		return cur;
	if (cur->leaf())
		cur = cur->p();

	shw("init");
	// cur - internal
	int i;
	int found = 0;
	RNode *MG = st->lca(cur->M(st), cur->p()->M(st));
	if (MG->leaf())
		return cur; // |L(G)|=1
	for (i = 0; i < 3; i++, cur = ((UNode3*) cur)->l())
		if (cur->M(st) != MG) {
			found = 1;
			break;
		}
	cur->mark();
	shw("ins");
	if (found) {
	  while (!cur->p()->leaf()) {
	    shw("wh");
	    UNode3 *cur3p = (UNode3*) cur->p();
	    if (cur3p->l()->M(st) != MG)
	      cur = cur3p->l();
	    else if (cur3p->r()->M(st) != MG)
	      cur = cur3p->r();
	    else {
	      cur = cur3p;
	      break;
	    }
	    cur->mark();
	  }
	  if (cur->M(st) != MG)
	    return cur;
	}
	cur->mark();
	for (i = 0; i < 3; i++, cur = ((UNode3*) cur)->l())
		if (cur->p()->M(st) == MG)
			return cur;
	return cur;
}

int lossprim(RNode *s, RNode *s1, RNode *s2) {
	if ((s != s1) && (s != s2))
		return s1->depth() + s2->depth() - 2 * s->depth() - 2;
	if (s != s1)
		return s1->depth() - s->depth();
	return s2->depth() - s->depth();
}

void dlcostdetintermediates(RNode *child, RNode *cur, RNode *last,
		int skiplast = 1) {
	while (1) {
		if ((cur == last) && (skiplast))
			return;
		if (child == ((RInt*) cur)->l())
			((RInt*) cur)->r()->costdet().loss++;
		else
			((RInt*) cur)->l()->costdet().loss++;
		if (cur == last)
			return;
		cur = cur->p();
		child = child->p();
	}
}

void dlcostdet(RNode *s, RNode *s1, RNode *s2) {
	//loss
	if ((s != s1) && (s != s2)) {
		dlcostdetintermediates(s1, s1->p(), s);
		dlcostdetintermediates(s2, s2->p(), s);
	} else {
		if (s != s1)
			dlcostdetintermediates(s1, s1->p(), s, 0);
		else if (s != s2)
			dlcostdetintermediates(s2, s2->p(), s, 0);
		s->costdet().dup++;
	}
}

//UNode *UTree::genRand(double pint, double dec, char **t, int s) {
//	if ((1.0 * rand() / RAND_MAX) < pint) {
//		UNode *a = genRand(pint * dec, dec, t, s);
//		UNode *b = genRand(pint * dec, dec, t, s);
//		return createNode3(a, b);
//	}
//	return createLeaf(strdup(t[rand() % s]));
//}


UNode *UTree::genRand(double pint, double dec, vector<string> &t, int s)
{
	if ((1.0*rand()/RAND_MAX)<pint)
	{
		UNode *a = genRand(pint*dec,dec,t,s);
		UNode *b = genRand(pint*dec,dec,t,s);
		return createNode3(a,b);
	}
	char *nm=strdup(t[rand()%s].c_str());
	return createLeaf( getspecies( nm) , nm);
}



void UTree::initrand(int len,double pint, double dec, vector<string> &t, int splen)
{
    UNode *cur = genRand(pint,dec,t,splen);
    UNode *cur2 = genRand(pint,dec,t,splen);
    for (int i=0; i<len-2; i++)
    {
	cur=createNode3(cur,cur2);
	cur2=genRand(pint,dec,t,splen);
    }
    cur->p(cur2);
    cur2->p(cur);
    start=cur;
}


UTree::UTree(int len,double pint, double dec, int numlv, int uniquelv, vector<string> &src)
{

	mmi("UT-vec");
	lsize=-1;

	// take all strings from src
	int splen=src.size();

	if ((numlv<0) && (!uniquelv))
	{
		initrand(len,pint,dec,src,splen);
	}
	else
	{

		int lf = numlv;
		if (numlv>0)
			if (uniquelv)
				if (splen<numlv) lf=splen; // no more
				else;
			else;
		else
		{ if (uniquelv)
			{
			if (splen>1) lf=1+(rand()%(splen));
			else lf=1;
			}
		}
		// only 2 parameters: lf - number of leaves to generate
		// uniquelv - are they unique?
		// Generating tree

		UNode *tb[lf];
		char lfusage[splen];

		// Usage of unique leaves
		int i;
		if (uniquelv)
			for (i=0;i<splen;i++) lfusage[i]=0;

		for (i=0;i<lf;i++)
		{
			// get a leaf (can be more efficient...)
			int pos=-1;
			while (1)
			{
				pos=rand()%splen;
				if (!uniquelv) break;
				if (!lfusage[pos])
				{
					// ok found
					lfusage[pos]=1;
					break;
				}
			
			}
			char *nm= strdup(src[pos].c_str());
			tb[i]=createLeaf( getspecies( nm),nm);
		}


		if (lf==1)
		{
			start = tb[0];
			return;
		}

		for (i=0;i<lf-2;i++)
		{
			// 	    cout << "lf-i=" << lf -i;
			// 	    cout << "  i=" << i << endl;
			int p = rand()%(lf-i);
			int q;
			do q = rand()%(lf-i);
			while (p==q);
			// 	    cout << " " << p << " " << q << endl;
			// 	    cout << "joining:" ;
			// 	    tb[p]->ppsmprooted(cout);
			// 	    cout << " with: ";
			// 	    tb[q]->ppsmprooted(cout);
			// 	    cout << endl;

			// join the trees
			tb[p] = createNode3(tb[p],tb[q]);

			int j;
			for (j=q+1;j<lf-i;j++) tb[j-1]=tb[j];

		}

		tb[0]->p(tb[1]);
		tb[1]->p(tb[0]);
		start=tb[0];

		// 	cout << " " << lf << " " << uniquelv << endl;
		// 	start->ppsmprooted(cout); cout << endl;
		// 	print(cout); cout << endl;

		//	cout << "=============================================" << endl;

	}
}


//
//UTree::UTree(int len, double pint, double dec, int numlv, int uniquelv,
//		char *src) {
//	int splen = strlen(src);
//	if ((numlv < 0) && (!uniquelv)) {
//		int i = 0;
//		char *t[splen];
//		char buf[2];
//		buf[1] = 0;
//		for (i = 0; i < splen; i++) {
//			buf[0] = src[i];
//			t[i] = strdup(buf);
//		}
//		initrand(len, pint, dec, t, splen);
//	} else {
//
//		int lf = numlv;
//		if (numlv > 0)
//			if (uniquelv)
//				if (splen < numlv)
//					lf = splen; // no more
//				else
//					;
//			else
//				;
//		else if (uniquelv)
//			{ if (splen > 1)
//				lf = 1 + (rand() % (splen));
//			else
//				lf = 1; }
//
//		// only 2 parameters: lf - number of leaves to generate
//		// uniquelv - are they unique?
//
//
//		// Generating tree
//
//		UNode *tb[lf];
//		char lfusage[splen];
//
//		// Usage of unique leaves
//		int i;
//		if (uniquelv)
//			for (i = 0; i < splen; i++)
//				lfusage[i] = 0;
//
//		for (i = 0; i < lf; i++) {
//			// get a leaf (can be more efficient...)
//			int pos = -1;
//			while (1) {
//				pos = rand() % splen;
//				if (!uniquelv)
//					break;
//				if (!lfusage[pos]) {
//					// ok found
//					lfusage[pos] = 1;
//					break;
//				}
//			}
//			char bf[2];
//			bf[1] = 0;
//			bf[0] = src[pos];
//			tb[i] = createLeaf(strdup(bf));
//		}
//
//		if (lf == 1) {
//			start = tb[0];
//			return;
//		}
//
//		for (i = 0; i < lf - 2; i++) {
//			// 	    cout << "lf-i=" << lf -i;
//			// 	    cout << "  i=" << i << endl;
//			int p = rand() % (lf - i);
//			int q;
//			do
//				q = rand() % (lf - i);
//			while (p == q);
//
//			// join the trees
//			tb[p] = createNode3(tb[p], tb[q]);
//
//			int j;
//			for (j = q + 1; j < lf - i; j++)
//				tb[j - 1] = tb[j];
//
//		}
//
//		tb[0]->p(tb[1]);
//		tb[1]->p(tb[0]);
//		start = tb[0];
//
//	}
//
//}

//UTree::UTree(int len, double pint, double dec, SpeciesTree *sp) {
//	int i = 0;
//	char *t[sp->lsize() + 1];
//	iterator_tree it(sp, F_LEAVES);
//	RLeaf *r;
//	while ((r = (RLeaf*) it()) != 0)
//		t[i++] = r->label();
//	t[sp->lsize()] = 0;
//	initrand(len, pint, dec, (char**) t, sp->lsize());
//}



int UTree::leaves()
{
  if (lsize>0) return lsize;
  if (!start->p()) return 1;
  lsize=start->leaves()+start->p()->leaves();
  return lsize;
}


UTree::UTree(char *t, int num,string tn, double wg, int extractspecies)
{
  //  mmi("UT-str");
  while (isspace(t[0])) t++;
  if ((strlen(t)>0)  && (((t[0]>='0') && (t[0]<='9')) || (t[0]=='.')))
    {
      if (sscanf(t,"%lf",&weight)!=1)
	{
	  cerr << "Weight expected in <<" << t << ">>" << endl;
	  exit(-1);
	}
      while (isdigit(t[0]) || isspace(t[0]) || t[0]=='.') t++;
    }
  else weight=1.0;
  int p=0;
  parseNode(t,p,1,num,extractspecies);
  lsize=-1;
  treename=tn;
  weight=wg;
}

UTree::UTree(int len,double pint, double dec, SpeciesTree *sp)
{
	mmi("UT-rand-st");
	int i=0;
	lsize=-1;
	vector<string> t;
	iterator_tree it(sp,F_LEAVES);
	RLeaf *r;
	while ((r=(RLeaf*)it())!=0) t[i++]=specnames[r->species()];
	initrand(len,pint,dec,t,sp->leaves());
}

void UNode::pcosts(ostream &s,DlCost &c,SpeciesTree *st)
{
  if (st)
    {
      s << " totalc({" << cost(st).dup << "," << cost(st).loss << "})"
	<< " treec({" << sc(st).dup << "," << sc(st).loss << "}) ";
      if (sc(st).dup==c.dup && sc(st).loss==c.loss)
	{
	  s << " minc(1) ";
	}
    }
  else
    {
		//s << " totalc({" << costuu.dup << "," << costuu.loss << "}) "
		//		<< " uusupport(" << (1.0*uusupport) << ") ";
		//if (pn)
		//{
		//	s << " uudiv(\"";
		//	rooted(1)->printsorted(s,T_LABS|T_TREEDIV);
		//	s <<"\") ";
		//}
	}


	if (ismarked & M_MARK) s << " mark(1)";
	if (ismarked & M_OPT) s << " markopt(1)";
	if (ismarked & M_START) s << " markstart(1)";
	if (ismarked & M_OPTM) s << " markoptm(1)";
	if (ismarked & M_MINC) s << " minc(1)";
	if (ismarked & M_ROOTING) s << " rooting(1)";
	if (st)
	{
		if (M(st)->p()) s << " destn(\"" << *M(st) << "\") ";
		else s << " destn(\"\") ";
	}
}

void UTree::center()
{
  //  cout << "centering\n";
  if (!start->p()) return;
  int all=leaves();

  UNode *cur = start;
  int lv = start->leaves();
  int lvp = all-lv;

  if (abs(lv-lvp) < 2) return;

  if (lv < lvp) {
    cur=start->p();
    lvp=lv;
    lv=all-lvp;
  }

//   cout << "lvp=" << lvp << endl;
//   cout << "lv=" << lv << endl;
//   cout << "====" << endl;

  while (1)
    {
      if (cur->leaf()) break;
      int lv_sl = ((UNode3*)cur)->l()->p()->leaves();
      int lv_sp = ((UNode3*)cur)->r()->p()->leaves();

//       cout << "lvp=" << lvp << endl;
//       cout << "lv_sl=" << lv_sl << endl;
//       cout << "lv_sp=" << lv_sp << endl;

      if (lv_sl > (lv_sp+lvp))
	{
	  cur=((UNode3*)cur)->l()->p();
	  lvp=lvp+lv_sp;
	  //	  cout << "left" << endl;
	}
      else
	if (lv_sp > (lv_sl+lvp))
	  {
	    cur=((UNode3*)cur)->r()->p();
	    lvp=lvp+lv_sl;
	    //	  cout << "rgh" << endl;
	  }
	else break;
      //      cout << "====" << endl;
    }
  start=cur;
  start->balance(0);
  start->p()->balance(0);
}
