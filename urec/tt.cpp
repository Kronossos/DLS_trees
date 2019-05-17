/************************************************************************
 Unrooted REConciliation version 1.00
 (c) Copyright 2005-2006 by Pawel Gorecki
 Written by P.Gorecki.
 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/
#include <set>
#include <vector>
using namespace std;
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include "tools.h"

#define OPT_PRINTGENE 1
#define OPT_PRINTSPECIES 2
#define OPT_BYCOST 8
#define OPT_WEAK 32
#define OPT_NNI 4
#define OPT_SPR 64
#define OPT_PRINTWEAKEDGES 128
#define OPT_STDERR 256

#define POSUMMARY 1 
#define POTAB 2
#define PODESCRIPTION 4
#define POFULL 8
#define POWEAKEDGES1 16
#define POWEAKEDGES2 32
#define POWEAKEDGES3 64

void spr(UTree *gt, RTree *sp, UNode *t, UNode3 *e, DlCost cost);

double MAXVAL=10e100;
extern int gsid;
extern char *gsdelim;
extern int gspos;
extern int usespecnames;
int genopt = 0;

double ttcompute(UTree *gt, SpeciesTree *st, int printopt, int ops)
{
	gt->clear();
	UNode *un = gt->findoptimaledge(st);
	DlCost &c = un->cost(st);
	double m = c.mut();

	if (printopt & POTAB)
	  {
	    cout << " s-cost " << m << " ";
	    un->ppsmprooted(cout);
	    cout << " | ";
	    un->p()->ppsmprooted(cout);
	    cout << endl;
	  }
	if (printopt & POFULL)
	  {
	    cout << m << " " << c.dup  << " " << c.loss << " " << ops << " (";
	    un->ppsmprooted(cout);
	    cout << ",";
	    un->p()->ppsmprooted(cout);
	    cout << ")" << endl;
	  }
	return m;
}

void _ttcycle( UNode3 *src ,UNode3 *dst, UTree *gt, SpeciesTree *st, int show)
{  
  UNode3 *nz= src->l();
  UNode3 *a = nz->l(); 
  UNode3 *b = dst->l(); 
  UNode3 *c = b->l(); 

#define LR(x,y) x->l(y); y->r(x)

  LR(a,b);
  LR(nz,c);
  LR(c,src);
  LR(dst,a);
  LR(b,dst);
}


UNode* _ttforward( UNode3 *src ,UNode3 *dst, UTree *gt, SpeciesTree *st, int show)
		{
	//  src--------src-p SUBTREE
	//  move SUBTREE to dst --------- dst-p

	UNode *x = src->l()->p();
	UNode *y = src->r()->p();
	double centerlen=x->len(); // to be reconstructed

	x->p(y);
	y->p(x);

	UNode *v = dst->p();
	src->l()->p(dst);
	dst->p(src->l());

	v->setbranchlen(centerlen);
	src->r()->p(v);
	v->p(src->r());
	return y;
 }

void _ttbackward( UNode3 *src ,UNode *y, UTree *gt, SpeciesTree *st, int show)
{
	UNode *x = y->p();
	UNode *dst = src->l()->p();
	UNode *v = src->r()->p();

	// center edge
	x->setbranchlen(v->len());
	src->l()->p(x);
	x->p(src->l());

	// len from y
	src->r()->p(y);
	y->p(src->r());

	// len from dst
	v->p(dst);
	dst->p(v);

}

double _ttbfsingle( UNode3 *src ,UNode3 *dst, UTree *gt, SpeciesTree *st, int show, int ops)
{
	UNode *y = _ttforward(src,dst,gt,st,show);
	double res=ttcompute(gt,st,show,ops);
	_ttbackward(src,y,gt,st,show);
	return res;
}

double _ttbfgen( UNode3 *src ,UNode *dst, UTree *gt, SpeciesTree *st, int spr, int show)
{
	if (dst->leaf()) {
		cout << "dst is a LEAF: ";
		cout << *dst->smprooted() << endl;
		return MAXVAL;
	}

	UNode3 *u = (UNode3*)dst;

	double res=_ttbfsingle(src,u->l(),gt,st, show, 1);
	double r=_ttbfsingle(src,u->r(),gt,st,show,1);
	if (res>r) res=r;

	if (spr)
	{
		r=_ttbfgen(src,u->l()->p(),gt,st,spr,show);
		if (res>r) res=r;
		_ttbfgen(src,u->r()->p(),gt,st,spr,show);
		if (res>r) res=r;
	}
	return res;
}

double _ttbf(UNode3 *src, UTree *gt, SpeciesTree *st, int spr, int show)
{
	// src ---- src-pn: CURRENT edge
	// initialize operation
	// take subtree: SUBTREE src-l and attach to all subtrees starting from
	// src->p, src-p-l, src->p-r, ....

	return _ttbfgen(src->l(),src->p(),gt,st,spr, show);
	//return res;
}

double ttbfspr(UTree *gt, SpeciesTree *st,
		 int show,double weak, int maxweak)
{
	cout << "SPR to be implemented" << endl;
	return 0.0;
}


vector<UNode3*> findweakedges(UTree *gt,SpeciesTree *st,double weak)	{
	vector<UNode3*> weakedges;
	UNode *ur;
	iterator_utree itu(gt);

	while ((ur = itu()) != 0) {
		if ((!ur->leaf()) && (!ur->p()->leaf()))
		{
			UNode3* u3= (UNode3*)ur;

			if (((weak<0) || (u3->len()<weak)))
			{
				int i,fnd=0;
				for (i=0;(int)i<(int)weakedges.size();i++)
					if (weakedges[i]->p()==u3)
					{
						fnd=1;
						break;
					}
				if (!fnd)
				{
					UNode3* pn=(UNode3*)u3->p();
					RNode *a=u3->M(st);
					RNode *b=pn->M(st);

//					cout << "CLUSTER CHECK: "
//							<< *u3->smprooted()
//							<< " | " << *pn->smprooted() <<  " ";

					int smp=0;
					// filter smpcluster edges
					if (a->leaf())
					{
						if (pn->l()->p()->M(st)==a) smp=1;
						else if (pn->r()->p()->M(st)==a) smp=1;
					}
					if ((!smp) && b->leaf())
					{
						if (u3->l()->p()->M(st)==b) smp=1;
						else if (u3->r()->p()->M(st)==b) smp=1;
					}
					if (!smp)
						//cout << "SMPCLUSTER!" << endl;
					//else
					//{
						//cout << "COMPLEX" << endl;
						weakedges.push_back(u3);

				}
			}
		}
	}
	return weakedges;
}

double ttbfnni(UTree *gt, SpeciesTree *st,
		int show, double weak, int maxweak,
	       int &err,int depth)
{
	iterator_utree itu(gt);
	double ermin=-1;

	// improve(!)
	vector<UNode3*> weakedges = findweakedges(gt,st,weak);
	if ((weak>=0) && (int)weakedges.size()>maxweak) {
		err=1;
		return 0;
	}

	int i;
	for (i=0;i<(int)weakedges.size(); i++)
	{
		UNode3* u3=weakedges[i];
		if (show)
		{
			cout << endl << "#subtree: ";
			cout << u3->len() << " " << u3->p()->len() << " ";
			u3->ppsmprooted(cout);
			cout << " | ";
			u3->p()->ppsmprooted(cout);
			cout << endl;
		}

		double r=_ttbfgen(u3->l(),u3->p(),gt,st,0, show);
		if ((ermin<0) || (ermin>r)) ermin=r;
	}
	return ermin;
}

int cnt=1;

void _tt2(int spr, UTree *gt, SpeciesTree *st,
	    int printopt, int depth, int ops, ofstream &f,
	    vector<UNode3*> &weakedges,
	    double &optcost, int &opterr)
{
   
#define indent() for (j=0;j<ops;j++) cout << "  ";

  int j;
  
  if (printopt & POTAB)
    {
      indent();
      cout << "solution " << ops << " ";
    }
  

  double res=ttcompute(gt,st,printopt,ops);
  if (!ops || (optcost>res) || (optcost==res && opterr>ops))
    {
      optcost=res;
      opterr=ops;
    }
  
#ifdef _DEBUG_
  if (printopt & POTAB) 
    {     
      f << "\\begin{figure}[h]" << endl;
      f << printpstricks(gt,NULL);
      f << "\\caption{solution ops " << ops << " cost= " << res;
      f << "}" << endl << "\\end{figure}" << endl << endl;
    }
#endif
  
  if (!depth) return; // last
    
  int i;


  for (i=0;i<(int)weakedges.size(); i++)
    {
      UNode3* src=weakedges[i];
      UNode* d=src->p();

      if (d->leaf()) {
	indent();
	cout << "dst is a LEAF: ";
	cout << *d->smprooted() << endl;
	continue;
      }
      
      UNode3 *dst= (UNode3*)d;


      if (!ops && genopt & OPT_STDERR) cerr << ".";

      
      if (printopt & POTAB) 
	{
	  indent();
	  cout << endl << "#subtree ops=" << ops << "   ";
	  cout << src->len() << " " << dst->len() << " ";
	  src->ppsmprooted(cout);
	  cout << " | ";
	  dst->ppsmprooted(cout);
	  cout << endl;

#ifdef _DEBUG_
	  f << "\\clearpage" << endl << endl << endl;
	  f << "\\begin{figure}[htbp!]" << endl;
	  f << printpstricks(gt,src);
	  f << "\\caption{";
	  f << endl << "subtree ops=" << ops << "   ";
	  f << src->len() << " " << dst->len() << " ";
	  src->ppsmprooted(f);
	  f << " | ";
	  dst->ppsmprooted(f);
	  f << "}" << "\\end{figure}" << endl << endl;
#endif
	}

#ifdef _DEBUG
      int cser=cnt;
#endif
      cnt++;


      // 1st NNI
      if (printopt & POTAB)
	{
#ifdef _DEBUG_
	  f << "\\begin{figure}[htbp!]" << endl;
	  f << printpstricks(gt,src);
	  f << "\\caption{";
	  f << "cnt=" << cser << " ";
	  f << endl << "A-1NNI before 1cycle" << ops << "   ";
	  f << "}" << "\\end{figure}" << endl << endl;
#endif
	}

      _ttcycle(src,dst,gt,st,printopt);

      if (printopt==1)
	{
#ifdef _DEBUG_
	  f << "\\begin{figure}[htbp!]" << endl;
	  f << printpstricks(gt,src);
	  f << "\\caption{";
	  f << "cnt=" << cser << " ";
	  f << endl << "B-1NNI after 1cycle" << ops << "   ";
	  f << "}" << "\\end{figure}" << endl << endl;
#endif
	}

      _tt2(spr, gt, st, printopt, depth-1,ops+1,f, weakedges,optcost,opterr);

      if (printopt & POTAB)
	{
#ifdef _DEBUG_
	  f << "\\begin{figure}[htbp!]" << endl;
	  f << printpstricks(gt,src);
	  f << "\\caption{";
	  f << "cnt=" << cser << " ";
	  f << endl << "C-1NNI  after tt/1cycle and before 2cycle" << ops << "   ";
	  f << "}" << "\\end{figure}" << endl << endl;
#endif
	}

      _ttcycle(src,dst,gt,st,printopt);

      if (printopt & POTAB)
	{
#ifdef _DEBUG_
	  f << "\\begin{figure}[htbp!]" << endl;
	  f << printpstricks(gt,src);
	  f << "\\caption{";
	  f << "cnt=" << cser << " ";
	  f << endl << "D-NNI after 2cycle" << ops << "   ";
	  f << "}" << "\\end{figure}" << endl << endl;
#endif
	}

      _tt2(spr, gt, st, printopt, depth-1,ops+1,f, weakedges,optcost,opterr);

      if (printopt & POTAB)
	{
#ifdef _DEBUG_
	  f << "\\begin{figure}[htbp!]" << endl;
	  f << printpstricks(gt,src);
	  f << "\\caption{";
	  f << "cnt=" << cser << " ";
	  f << endl << "E-NNI after tt/2cycle - before 3cycle" << ops << "   ";
	  f << "}" << "\\end{figure}" << endl << endl;
#endif
	}

      _ttcycle(src,dst,gt,st,printopt);

      if (printopt & POTAB)
	{
#ifdef _DEBUG_
	  f << "\\begin{figure}[htbp!]" << endl;
	  f << printpstricks(gt,src);
	  f << "\\caption{";
	  f << "cnt=" << cser << " ";
	  f << endl << "F-NNI after after 3cycle" << ops << "   ";
	  f << "}" << "\\end{figure}" << endl << endl;
#endif
	}
    }

}

double tt2(int spr, UTree *gt, SpeciesTree *st,
	   int &numerrors, int &removedweaktrees,
	   int printopt, double weak, int maxweak,int depth, ofstream &f)
{

  vector<UNode3*> weakedges = findweakedges(gt,st,weak);


  if (genopt & OPT_STDERR)
    cerr << " "<< weakedges.size() << " ";
  
  if ((weak>=0) && (int)weakedges.size()>maxweak) {
    removedweaktrees++;
    return 0;
  }    

  double optcost;
  int opterr;

  if (printopt & PODESCRIPTION && printopt & POFULL)
    {
      cout << "MutationCost Duplications Losses NNI-distance(errors) Rooted-opt-gene-tree" << endl;
    }
  
  _tt2(spr, gt, st, printopt, depth,0,f, weakedges,optcost,opterr); 

  if (printopt & PODESCRIPTION && printopt & POFULL)
    {
      cout << "Please note that it can be more optimal rooted gene trees with the same cost" << endl;
    }

  numerrors+=opterr;

  return optcost;
  
  
}

// Old version for single NNI
double tt(int spr, UTree *gt, SpeciesTree *st,
		double &costexact, double &costtt, int &numerrors, int &removedweaktrees,
	  int printopt, double weak, int maxweak,int depth)
{
	//if (printopt) cout <<"#initial tree:"<<endl;
	int err=0; // set if too many weak

	//cout << " " << weak << " " << maxweak << endl;
	if (spr)
	  costtt=ttbfspr(gt,st,printopt,weak,maxweak); // not implemented yet
	else
	  costtt=ttbfnni(gt,st,printopt,weak,maxweak,err,depth);
	
	if (err)
	  {
	    removedweaktrees++;
	    return 0;
	  }
	if (printopt) cout << "exact ";
	
	costexact=ttcompute(gt,st,printopt,0);
	if ((costtt>=0) && (costexact>costtt))
	{
		numerrors+=1;
		return costtt;
	}
	return costexact;

}

int usage(int argc, char **argv) {
  cout << 
    " Error corrections and rootings evaluation. k-NNI algorithms.\n"
    " Usage: " << argv[0] << " [options]\n"
    " -g gene tree \n"
    " -s species tree\n"
    " -G filename - defines a set of gene trees\n"
    " -S filename - defines a set of species trees\n"
    " -p - print a gene tree\n"
    " -P - print a species tree\n"
    " -D dupweight  - set weight of gene duplications\n"
    " -L lossweight - set weight of gene losses\n"
    " -k DEPTH - how many errors to correct (NNI operations, default 1)\n"
    //	" -n run NNI\n"
    " -b min cost GS vs ST\n"
    " -R [dst]+ - printing options for error correction where"
    "     s - summary only (default)\n"
    "     d - show descriptions\n"
    "     f - complete information for all NNI variants\n"
    "     u - print species names instead of full gene ids in gene trees\n"
    "     w - for each gene tree print the number of gene edges\n"
    "     W - print the gene trees filtered by -m and -w parameters\n"
    "     v - print only the numbers of gene trees filtered by -m and -w parameters\n"
    
    //	  "     t - (debug) tab separated for all NNI variants\n"
    
    " -r TREE|SPECIESLIST initialize species order\n"
    " -w RNUM nni analysis with weak edges (branch lengths required, omega parameter)\n"
    " -m max num of weak edges (valid with -w only, mu parameter)\n"
    " -v info on stderr after every gene tree\n"
    //    " -W e|f|i weak edges info (e-full info, f-trees filtered by -m, i-tree number only)\n"
    
    " -l pNUM|p-NUM|aDELIM|bDELIM - mapping gene identifiers to species names\n"
    "   pNUM - species name in the first NUM characters of gene ids\n"
    "   p-NUM - species name in the last NUM characters of gene ids\n"
    "   aDELIM - species name after delimiter DELIM\n"
    "   bDELIM - species name before delimiter DELIM\n"

    "\n Examples: \n"
    "> Full error correction (all edges are candidates)\n"
    "tt -G genetrees.txt -S speciestree.txt -b\n\n"
    
    "> Compute min cost with weak edges and gene tree filtering\n"
    "tt -G genetrees.txt -S speciestree.txt -m4 -w0.1 -b\n\n"	  
    
    "> Compute min cost with no weak edges \n"
    "tt -G genetrees.txt -S speciestree.txt -w0 -b\n\n"	  
    
    "> Print counts of weak edges in each gene tree\n"
    "tt -G genetree.txt -S speciestree.txt -m4 -w0.1 -Re\n"
  
    "> Print detailed NNI variants for the gene-species tree\n"
    "tt  -G genetree.txt  -S speciestree.txt -b -k1 -Rdf\n\n"

    "> Print detailed variants (-Rf) for correction of at most 2 errors (-k2) and sort the results\n"
    "tt  -G genetree.txt  -S speciestree.txt -b -k2 -Rf | sort -k1 -n \n\n"   

    "> Gene ids have the following from SPECIESNAME_GENEDATA\n"
    "tt -Rf -la_ -g '((cat_142:2,(mouse_123:10,mouse_1444:12):13):4,rat_3333asdaasdasd_sd:1):3'  -s '(cat,(rat,mouse))'  -b \n\n"
      ;

  exit(-1);
}

int main(int argc, char **argv) {
	int opt;
	int printopt = 0;

	double weak=-1;
	int maxweak=2;
	int depth=1;
	if (argc < 2) usage(argc, argv);
	vector<SpeciesTree*> stvec;
	vector<UTree*> gtvec;

	srand(time(0));

	gsid=GSFULL;
	while ((opt = getopt(argc, argv, "s:S:g:G:pPL:D:bR:r:w:m:t:W:k:vl:"))
			!= -1)
		switch (opt) {
		case 'g':
		  gtvec.push_back(new UTree(optarg,-1,"",1,1));
		  break;
		case 's':
		  stvec.push_back(new SpeciesTree(optarg,-1));
		  break;
		case 'S':
		  readstree(optarg, stvec);
		  break;
		case 'l':
		  switch (optarg[0])
		    {
		    case 'p': 
		      gspos=atoi(optarg+1);
		      gsid=GSPOS;
		      break;
		    case 'a':
		      gsid=GSAFTER;
		      gsdelim=strdup(optarg+1);
		      break;
		    case 'b':
		      gsid=GSBEFORE;
		      gsdelim=strdup(optarg+1);
		      break;
		    default:
		      cerr << "Invalid code for -l option expected aDELIM,bDELIM or p[-]number" << endl;
		      exit(-1);
		    }
		  break;
		  
		case 'G':
			readgtree(optarg, gtvec);
			break;
		case 'p':
			genopt |= OPT_PRINTGENE;
			break;
		case 'P':
			genopt |= OPT_PRINTSPECIES;
			break;

		case 'b':
			genopt |= OPT_BYCOST;
			break;

		case 'R':
		  {
		    int i;
		    for (i=0;i<int(strlen(optarg));i++)
		      {
			switch (optarg[i])
			  {
			  case 's': printopt|=POSUMMARY; break;
			  case 't': printopt|=POTAB; break;
			  case 'd': printopt|=PODESCRIPTION; break;
			  case 'f': printopt|=POFULL; break;
			  case 'w': printopt|=POWEAKEDGES1;	 break;
			  case 'W': printopt|=POWEAKEDGES2;	break;

			  case 'v': printopt|=POWEAKEDGES3;	break;
			  case 'u': usespecnames=1; break;
			  default:
			    cerr << optarg[i] << " unknown suboption in -R" << endl;
			    exit(-1);
			  }
		      }
		  }
		  break;

		case 'v':
		  genopt |= OPT_STDERR;
		  break;

		case 'r':
		{
			initgenrand(optarg);
			break;
		}

		case 'L':
			if (sscanf(optarg, "%lf", &weight_loss) != 1) {
				cerr << "Number expected in -L" << endl;
				exit(-1);
			}
			break;
		case 'D':
			if (sscanf(optarg, "%lf", &weight_dup) != 1) {
				cerr << "Number expected in -D" << endl;
				exit(-1);
			}
			break;

		case 't':
			if (optarg[0]=='s') genopt |= OPT_SPR;
			else genopt|=OPT_NNI;
			break;

		case 'w':
			genopt |= OPT_WEAK;
			if (sscanf(optarg, "%lf", &weak) != 1) {
				cerr << "Number expected in -w" << endl;
				exit(-1);
			}
			break;

		case 'k':
			if (sscanf(optarg, "%d", &depth) != 1) {
				cerr << "Number expected in -k" << endl;
				exit(-1);
			}
			break;

		case 'm':
			if (sscanf(optarg, "%d", &maxweak) != 1) {
				cerr << "Number expected in -m" << endl;
				exit(-1);
			}
			break;



		default:
			cerr << "Unknown option: " << ((char) opt) << endl;
			exit(-1);
		}

	vector<SpeciesTree*>::iterator stpos;
	vector<UTree*>::iterator gtpos;


	if (genopt & OPT_PRINTGENE)
		for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)
			(*gtpos)->print(cout) << endl;


	if (printopt & (POWEAKEDGES1|POWEAKEDGES2|POWEAKEDGES3))	  
	  for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)
	    {
	      int id=0;

	      if (printopt & PODESCRIPTION) 
		{
		  cout << "Species tree: ";
		  (*stpos)->print(cout) << endl;

		  if (printopt & POWEAKEDGES1)
		    cout << "NumberOfWeakEdges "; 

		  if (printopt & (POWEAKEDGES2|POWEAKEDGES3))
		    {		      
		      if (printopt & POWEAKEDGES3) cout << "TreeNum ";		      
		      if (printopt & POWEAKEDGES2) cout << "GeneTree ";
		    }
		  else
		    if (printopt & POWEAKEDGES1)
		      if (printopt & POWEAKEDGES2) cout << "GeneTree ";

		  cout << endl;
		}
	      
	      for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)
		{	

		  int k=findweakedges(*gtpos,*stpos,weak).size(); 

		  if (printopt & (POWEAKEDGES2|POWEAKEDGES3))
		    {		      
		      if (k<=maxweak)
			{
			  if (printopt & POWEAKEDGES1)
			    cout << k << " "; 

			  if (printopt & POWEAKEDGES3)
			    cout << id << " ";

			  if (printopt & POWEAKEDGES2)
			    (*gtpos)->start->printrooting(cout,1); 
			  else cout << endl;
			  
			}
		      
		    } 
		  else		  
		    if (printopt & POWEAKEDGES1)
		      {		      
			cout << k << " ";
			(*gtpos)->start->printrooting(cout,1); 
		      }
		  
		  id++;
		}		    
	      exit(0);
	    }
	
	if (genopt & OPT_PRINTSPECIES) {
	  for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)
	    (*stpos)->print(cout) << endl;
	}
	

	if (genopt & OPT_BYCOST) {

	  ofstream f;
#ifdef _DEBUG_
	  if (printopt & POTAB) 
	    {
	      f.open ("tree.tex");  
	      for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos) {
		f << printpstricks(*gtpos,NULL);
	      }
	    }
#endif
	  
	  int cnts=0;
	  for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos) {
	    SpeciesTree *s = *stpos;
	    
	    int errnum=0;
	    double total=0;
	    int removedweaktrees=0;
	    int cnt=0;

	    cnts++;
	    if (!printopt) printopt=POSUMMARY;

	    for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos) {
	      UTree *g = *gtpos;
	      g->clear();
	      cnt++;
	      if (genopt & OPT_STDERR)
		{
		  cerr << cnts << "/" << stvec.size() << "-";
		  cerr << cnt << "/" << gtvec.size();
		}
	      total+=tt2(genopt & OPT_SPR, g,
			 *stpos,
			 errnum,
			 removedweaktrees,
			 printopt,weak,maxweak,depth,f);			
	    }

	    if (genopt & OPT_STDERR)
	      {
		cerr << endl;
	      }
	    
	    if (printopt & PODESCRIPTION && printopt & POSUMMARY)
	      cout << "OptimalCost NumberOfErrors(NNI-dist) RejectedTrees SpeciesTree" << endl;
	    if (printopt & POSUMMARY)
	      cout << total << " " << errnum << " " << removedweaktrees << " "<< *s << endl;	  
	  }

	  if (printopt & POTAB)
	    f.close();
	  
	  
	  
	}
	
	else
	  {
		int j = 0;
		for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos) {
			j++;
			UTree *g = *gtpos;
			double costexact;
			double costtt;
			int errnum=0;
			int removedweaktrees=0;
			for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)
			  tt(genopt & OPT_SPR, g,*stpos,costexact,costtt,
			     errnum,
			     removedweaktrees,
			     printopt,
			     weak,maxweak,depth);

		}
	}

}
