/************************************************************************
 Unrooted REConciliation version 1.00
 (c) Copyright 2005-2006 by Pawel Gorecki
 Written by P.Gorecki.
 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/
#include <string.h>
#include "tools.h"

#define OPT_RECDETAILS 1
#define OPT_RECINFO 2
#define OPT_PRINTROOTINGS 4
#define OPT_RECMINROOTING 8
#define OPT_RECTREECOSTDETAILS 16
#define OPT_PRINTGENE 32
#define OPT_PRINTSPECIES 64
#define OPT_RECMINCOST 128
#define OPT_SUMMARYTOTAL 256
#define OPT_SUMMARYDISTRIBUTIONS 512
#define OPT_SUMMARYDLTOTAL (1<<12)
#define OPT_TREEDISTRIBUTIONS (1<<13)
#define OPT_VOTING (1<<14)
#define OPT_BYCOST (1<<15)
#define OPT_RANDUNIQUE (1<<16)
#define OPT_PRINTBRANCHLENGHTS (1<<17)
#define OPT_RECMINROOTINGEXT (1<<18)
#define OPT_CLUSTERS (1<<19)
#define OPT_RECDETAILSEXT (1<<20)
#define OPT_RECEPIS (1<<21)

#define OPT_NNI (1<<25)
#define OPT_RSORT (1<<26)

extern int gspos;
extern int gsid;
extern char* gsdelim;

int adjacent(UNode *a, UNode *b)
{
  if (a==b) return 1;
  if (a->leaf()) return 0;
  UNode3 *ap=(UNode3*)a;
  return ap->r()==b || ap->r()->r()==b;
}

int usage(int argc, char **argv) {
  cout << " Unrooted REConciliation v2.00. (C) P.Gorecki 2005-2010" << endl;

  cout << " Usage: " << argv[0] << " [options]" << endl;
  cout << " -g gene tree " << endl;
  cout << " -s species tree" << endl;
  cout << " -G filename - defines a set of gene trees" << endl;
  cout << " -S filename - defines a set of species trees" << endl;
  cout << " -T filename - read file in interleaved format: gene_tree <EOLN> species_tree <EOLN> ..." << endl;

  cout << " -R - show rootings for every gene tree" << endl;
  cout << " -p - print a gene tree" << endl;
  cout << " -P - print a species tree" << endl;
  cout << " -D dupweight  - set weight of gene duplications" << endl;
  cout << " -L lossweight - set weight of gene losses" << endl;
  cout << " -r leaves - random unrooted gene trees" << endl;
  cout << "   -n num - length for random gene tree" << endl;
  cout << "   -i rnum - prob. of internal node" << endl;
  cout << "   -e rnum - decreased prob. of internal node" << endl;
  cout << "   -l num - number of random gene trees" << endl;
  cout << "   -u - unique leaves (a species tree)" << endl;
  cout << "   -E num - number of leaves (0 - use all species)" << endl;
  cout << " -b - computing costs" << endl;
  cout << " -y - show branch lengths (preserved during unrooted -> rooted convertion)";
  cout << " -m pNUM|p-NUM|aDELIM|bDELIM - mapping gene identifiers to species names" << endl 
       << "    pNUM - species name in the first NUM characters of gene ids" << endl
       << "    p-NUM - species name in the last NUM characters of gene ids" << endl
       << "    aDELIM - species name after delimiter DELIM" << endl
       << "    bDELIM - species name before delimiter DELIM" <<endl; 
  cout
    << " For every reconciliation of an unrooted gene tree with a species tree (details of costs):"
    << endl;
  cout << "   -o - show an optimal cost" << endl;
  cout << "   -O - show an optimal rooting" << endl;
  cout << "   -I - show all non-root-duplication rooted subtrees" << endl;
  cout << "   -a1 - show attributes and mappings" << endl;
  cout << "   -a2 - show cluster & split representation" << endl;
  cout << "   -a3 - show detailed attributes" << endl;
  cout << "   -a4 - show all rootings with costs" << endl;
  cout << "   -a5 - show candidates for episode rooting (up to 5 edges)" << endl;
  cout << "   -a51 - show the rooting on empty edge or nothing" << endl;
  cout << "   -a52 - show the rooting on double edge or max 2 candidates for empty edge gene trees" << endl;
  cout << "   -a53 - show all rootings with types: D-double, E-empty, L-leftempty, R-rightempty " << endl;
  cout << "   -a6 - show dl-plateau size" << endl;
  cout << "   -a7 - show d-plateau size" << endl;
  cout << "   -a8 - show all dl-optimal rootings" << endl;
  cout << "   -a81 - empty edge not-incident to a leaf " << endl;
  cout
    << " For every species tree, i.e., summary of costs when reconciling a species tree with a set of gene trees):"
    << endl;
  cout << "   -c - print total mutation cost" << endl;
  cout << "   -C - print total dl-cost (dup,loss)" << endl;
  cout << "   -d - print detailed total cost (distributions)" << endl;
  cout << "   -x - print species tree with detailed total costs (nested parenthesis notation with attributes)"
       << endl;
  
  cout << "   -N - print nni variants for rooted (use -s/-S) or unrooted (-g/G) species trees" << endl;
  
  cout << " -Mg - graphviz output, example: " << endl;
  cout << "          urec -Mg -G gtree.txt -p | dot -Tps > gtree.ps"  << endl;
  
  
  exit(-1);
}


int plateauborder(UNode *u, RNode *mg)
{
  if (u->leaf()) return 1;
  UNode3 *u3=(UNode3*)u;
  return u3->l()->p()->M(0)!=mg && u3->r()->p()->M(0)!=mg;
 
}

const char * nodetype(UNode *u)
{
  if (u->leaf()) return "L";
  UNode3 *u3=(UNode3*)u;
  if (u3->l()->p()->M(0)==u3->M(0) || u3->r()->p()->M(0)==u3->M(0))
    return "D";
  return "S";     
  
}

UNode *check2epis(UNode* un, SpeciesTree *s, int optdup)
{ 
  if (un->leaf()) return NULL;
  
  UNode *r=((UNode3*)un)->r()->p();
  DlCost dlu=r->cost(s);
  if (dlu.dup!=optdup) return NULL;
  // in D-plateau
  UNode *l=((UNode3*)un)->l()->p();		      
  //cout << *r->smprooted() << "----" << *l->smprooted() << endl;
  //cout << *r->M(s) << "----" << *l->M(s) << "     " << *un->M(s)   << endl;
  if (un->M(s)==l->M(s)) return l; 
  return r; 
}

void pinfo(UNode *u, RNode *mg,SpeciesTree *s)
{
  if (plateauborder(u,mg)) cout << "Bor ";
  else cout << "Plt ";
  cout << nodetype(u);
  cout << " " << *u->smprooted() << " " << *u->M(s) << endl;
}

int main(int argc, char **argv) {
  int opt;
  int rt_len = 2;
  int rt_numlv = -1;
  int loop = 10;
  double rt_pint = 0.5;
  double rt_dec = 0.75;
  
  if (argc < 2)
    usage(argc, argv);
  vector<SpeciesTree*> stvec;
  vector<UTree*> gtvec;
  srand(time(0));
  
  int genopt = 0;
  int optclusters=0;
  while ((opt = getopt(argc, argv, "T:M:bvg:s:pPE:ua:r:Rl:i:e:n:IOoG:XcCdxL:D:S:yNXm:"))
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
    case 'G':
      readgtree(optarg, gtvec);
      break;
    case 'T':
      readtreesinterleaved(optarg, gtvec, stvec);
      break;



      
    case 'p':
      genopt |= OPT_PRINTGENE;
      break;
    case 'v':
      genopt |= OPT_VOTING;
      break;
    case 'P':
      genopt |= OPT_PRINTSPECIES;
      break;
      
    case 'N':
      genopt|=OPT_NNI;
      break;
    case 'l':
      if (sscanf(optarg, "%d", &loop) != 1) {
	cerr << "Number expected in -l" << endl;
	exit(-1);
      }
      break;
    case 'r':
      {
	vector<string>  ss=initgenrand(optarg);
	int i;
	if (rt_numlv==0)
	  rt_numlv=specnames.size();
	for (i=0; i<loop; i++)
	  gtvec.push_back(new UTree(rt_len,rt_pint,rt_dec,rt_numlv,(genopt&OPT_RANDUNIQUE), ss));
	break;
      }
    case 'n':
      if (sscanf(optarg, "%d", &rt_len) != 1) {
	cerr << "Number expected in -l" << endl;
	exit(-1);
      }
      break;
      
    case 'i':
      if (sscanf(optarg, "%lf", &rt_pint) != 1) {
	cerr << "Number expected in -i" << endl;
	exit(-1);
      }
      break;
    case 'e':
      if (sscanf(optarg, "%lf", &rt_dec) != 1) {
	cerr << "Number expected in -e" << endl;
	exit(-1);
      }
      break;
    case 'E':
      if (sscanf(optarg, "%d", &rt_numlv) != 1) {
	cerr << "Number expected in -E" << endl;
	exit(-1);
      }
      break;

    case 'u':
      genopt |= OPT_RANDUNIQUE;
      break;
    case 'b':
      genopt |= OPT_BYCOST;
      break;
    case 'a':
      {
	int atype=1;
      if (sscanf(optarg, "%d", &atype) != 1) {
	cerr << "Number expected in -a" << endl;
	exit(-1);
      }
      if (atype==1) genopt |= OPT_RECINFO;
      else if (atype==3) genopt |= OPT_RECINFO | OPT_RECDETAILS;
      else if (atype<9 || atype==81 || atype==51 || atype==52 || atype==53) optclusters=atype;
      else {
	cerr << "Expected 1,2,...,8" << endl;
	exit(-1);
      }
      break;
      }

    case 'm':
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
	  cerr << "Invalid code for -m option expected aDELIM,bDELIM or p[-]number" << endl;
	  exit(-1);
	}
      break;

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

    case 'R':
      genopt |= OPT_PRINTROOTINGS;
      break;
      
    case 'o':
      genopt |= OPT_RECMINCOST;
      break;
      
    case 'O':
      genopt |= OPT_RECMINROOTING;
      break;
    case 'I':
      genopt |= OPT_RECMINROOTINGEXT;
      break;
      
    case 'X':
      genopt |= OPT_RECTREECOSTDETAILS;
      break;
    case 'y':
      genopt |= OPT_PRINTBRANCHLENGHTS;
      
    case 'c':
      genopt |= OPT_SUMMARYTOTAL;
      break;
      
    case 'C':
      genopt |= OPT_SUMMARYDLTOTAL;
      break;
      
    case 'd':
      genopt |= OPT_SUMMARYDISTRIBUTIONS;
      break;
      
    case 'x':
      genopt |= OPT_TREEDISTRIBUTIONS;
      break;
      
      
    case 'M':
      // sort nodes of rooted trees during reading
      sortmode(optarg);
      break;
      
      
    default:
      cerr << "Unknown option: " << ((char) opt) << endl;
      exit(-1);
    }
  
  vector<SpeciesTree*>::iterator stpos;
  vector<UTree*>::iterator gtpos;
  
  if (genopt & OPT_PRINTGENE)
    {
      for (gtpos=gtvec.begin(); gtpos !=gtvec.end(); ++gtpos)
	{
	  if (gcenter) (*gtpos)->center();
	  if (rsort) (*gtpos)->normalize()->printsorted(cout,0) << endl;
	  else
	    if (ppgraphviz)
	      (*gtpos)->printgraphviz(cout);
	    else (*gtpos)->print(cout) << endl;
	}
    }
  

  if (genopt & OPT_NNI)
    {
      // NNI operations for the set of gene/species trees
      nnist(stvec);
      nnigt(gtvec);
    }
  
  if (genopt & OPT_PRINTSPECIES) {
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)
      (*stpos)->print(cout) << endl;
  }
  
  if (genopt & OPT_PRINTROOTINGS) {
    for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)
      (*gtpos)->printrootings(cout);
  }
  
  if (genopt & OPT_VOTING) {
    
    int trnum = stvec.size();
    double mincnts[trnum];
    int i;
    
    for (i = 0; i < trnum; i++)
      mincnts[i] = 0;
    
    int j = 0;
    for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos) {
      j++;
      double min = 0;
      int minc = 0;
      UTree *g = *gtpos;
      i = 0;
      for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos) {
	g->clear();
	UNode *un = g->findoptimaledge(*stpos);
	double m = (un->cost(*stpos)).mut();
	if (i == 0) {
	  min = m;
	  minc = 1;
	} else if (min > m) {
	  min = m;
	  minc = 1;
	} else if (min == m)
	  minc++;
	i++;
      }
      i = 0;
      for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos) {
	g->clear();
	UNode *un = g->findoptimaledge(*stpos);
	if ((un->cost(*stpos)).mut() == min)
	  mincnts[i] += 1.0 / minc;
	i++;
      }
      
    }
    i = 0;
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)
      cout << **stpos << " " << mincnts[i++] << endl;
  }
  
  if (genopt & OPT_BYCOST) {
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos) {
      SpeciesTree *s = *stpos;
      if (genopt & OPT_RECINFO)
	cout << " SPECIES TREE: " << endl << *s << endl;
      
      DlCost total;
      
      for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos) {
	UTree *g = *gtpos;
	g->clear();
	
	
	if (genopt & OPT_RECINFO) {
	  cout << " GENE TREE: " << endl;
	  iterator_utree itu(g);
	  UNode *ur;
	  while ((ur = itu()) != 0) {
	    if (ur->leaf())
	      cout << "** leaf " << ((ULeaf*) ur)->species();
	    else {
	      cout << "** int  ";
	      if (genopt & OPT_RECDETAILS)
		cout << "  " << *ur->smprooted() << endl;
	    }
	    if (genopt & OPT_RECDETAILS)
	      cout << "  p=" << *ur->p()->smprooted() << endl;
	    cout << "\t sc=" << ur->sc(s);
	    cout << "\t cost=" << ur->cost(s) << "\t ";
	    cout << *ur->smprooted() << " ==> " << *ur->M(s)
		 << endl;
	  }
	}
	UNode *un = NULL;


	if (optclusters)
	  {

	  iterator_utree itu(g);
	  UNode *ur,*urp;
	  un = g->findoptimaledge(s);
	  UNode *unp = un->p();
	  DlCost dl=un->cost(s);
	  int dlplateausize=0;
	  int dplateausize=0;
	  while ((ur = itu()) != 0) {
	    DlCost dlu=ur->cost(s);
	    if (dlu.dup==dl.dup && dlu.loss==dl.loss) dlplateausize++;
	    if (dlu.dup==dl.dup) dplateausize++;
	  }
	  dlplateausize/=2;
	  dplateausize/=2;
	  RNode *MG = s->lca(un->M(s), un->p()->M(s)); // M(G)

#define PROOTINGS(ur) cout << "(" << *ur->smprooted() << "," << *ur->p()->smprooted() << ") "
#define PROOTINGSCOST(ur,dlu) PROOTINGS(ur) <<  dlu.dup << " "<< dlu.loss << " " <<  dlu.dup+dlu.loss<<endl;

	  if (optclusters==6)
	    cout << dlplateausize << endl;
	  else if (optclusters==7)
	    cout << dplateausize << endl;
	  else if (optclusters==81)
	    {
	      if (dlplateausize==1 && un->M(s)!=MG && !un->leaf() && !unp->leaf())
		cout << "Internal empty edge" << endl;
	    }
	  else if (optclusters==51)
	    {
	      if (dlplateausize==1 && un->M(s)!=MG) 
		PROOTINGS(un)<<endl; 
	    }
	  else if (optclusters==52)
	    {
	      // EPIS gen
	      if (dlplateausize==1 && un->M(s)!=MG) 
		{ // EMPTY EDGE 
		  // Find other connected
		  ur=check2epis(un,s,dl.dup);
		  if (ur) PROOTINGS(ur)<<endl;

		  ur=check2epis(unp,s,dl.dup);
		  if (ur) PROOTINGS(ur)<<endl;

		}
	      else {
		PROOTINGS(un)<<endl; 
	      }
	    }

	  else if (optclusters==53)
	    {
	      // EPIS gen
	      if (dlplateausize==1 && un->M(s)!=MG) 
		{ // EMPTY EDGE
		  cout << "E ";
		  PROOTINGS(un)<<endl; 
		  // Find other connected
		  
		  ur=check2epis(un,s,dl.dup);
		  if (ur) { cout << "L "; PROOTINGS(ur)<<endl; }

		  ur=check2epis(unp,s,dl.dup);
		  if (ur) { cout << "R "; PROOTINGS(ur)<<endl; }

		}
	      else {
		cout << "D "; PROOTINGS(un)<< endl; 
	      }

	    }
	  else if (optclusters==5)
	    {
	      // EPIS gen
	      if (dlplateausize==1 && un->M(s)!=MG) 
		{ // EMPTY EDGE 
		  // Find other connected
		  iterator_utree itu2(g);   	  
		  while ((ur = itu2()) != 0) {
		    DlCost dlu=ur->cost(s);
		    if (ur->marked() & M_PROCESSED) continue;
		    if (dlu.dup!=dl.dup) continue;
		    urp=ur->p();
		    if (adjacent(ur,un) || adjacent(urp,un) || adjacent(ur,unp) || adjacent(urp,unp))
		      {
			PROOTINGS(ur) << endl;
			//PROOTINGSCOST(ur,dlu);
		      }
		    ur->mark(M_PROCESSED);
		    ur->p()->mark(M_PROCESSED);
		  }
		}
	      else {
		//PROOTINGSCOST(un,dl); 
			PROOTINGS(un)<<endl; 
	      }

	    }
	  else
	  if (optclusters==4)
	    {
	      cout << "GENE TREE: " << endl;
	      cout << "Optimal cost: " << dl << endl;
	      cout << "DL-plateau size: " << dlplateausize << endl;
	      cout << "D-plateau size: " << dplateausize << endl;	  	
	      
	      // By edges

	      iterator_utree itu2(g);   	  
	      while ((ur = itu2()) != 0) {
		DlCost dlu=ur->cost(s);
		if (ur->marked() & M_PROCESSED) continue;
		PROOTINGSCOST(ur,dlu);
		ur->mark(M_PROCESSED);
		ur->p()->mark(M_PROCESSED);
		//pinfo(ur,MG,s);
		//pinfo(ur->p(),MG,s);
	      }
	    }
	  else if (optclusters==8)
	    {
	      // opt. dl-rootings
	      iterator_utree itu2(g);   	  
	      while ((ur = itu2()) != 0) {
		DlCost dlu=ur->cost(s);
		if (ur->marked() & M_PROCESSED) continue;
		if (dlu.dup==dl.dup && dlu.loss==dl.loss)
		  {			
		    ur->mark(M_PROCESSED);
		    ur->p()->mark(M_PROCESSED);
		    PROOTINGS(ur)<<endl;
		  }
	      }
	    }	  
	  else if (optclusters==2)
	    {
	      cout << "GENE TREE: " << endl;
	      cout << "Optimal cost: " << dl << endl;
	      cout << "Plateau size: " << dlplateausize << endl;

	      if (dlplateausize==1 && un->M(s)!=MG) cout << "Root S";
	      else cout << "Root D";
	      cout << " (" << *un->smprooted() << "," << *un->p()->smprooted() << ")" 
		   << " " << *MG << endl;
	  
	      iterator_utree itu2(g);   	  
	      while ((ur = itu2()) != 0) {
		DlCost dlu=ur->cost(s);
		if (ur->marked() & M_PROCESSED) continue;
		if (dlu.dup==dl.dup && dlu.loss==dl.loss)
		  {			
		    ur->mark(M_PROCESSED);
		    ur->p()->mark(M_PROCESSED);
		    pinfo(ur,MG,s);
		    pinfo(ur->p(),MG,s);
		  }
	      }
	      
	      iterator_utree itu3(g);   	  
	      while ((ur = itu3()) != 0) {
		DlCost dlu=ur->cost(s);
		if (ur->marked() & M_PROCESSED) continue;
		if (dlu.dup!=dl.dup || dlu.loss!=dl.loss)
		  {			
		    if (ur->M(s)!=MG)
		      {
			cout << "Ext " << nodetype(ur) << " ";
			cout << *ur->smprooted() << " " << *ur->M(s) << endl;
		      }
		  }	    
	      }
	    }	  
	}	
	
	if (genopt & (OPT_RECMINROOTING | OPT_RECMINROOTINGEXT | OPT_RECMINCOST
		      | OPT_RECTREECOSTDETAILS | OPT_SUMMARYTOTAL
		      |OPT_SUMMARYDLTOTAL | OPT_SUMMARYDISTRIBUTIONS
		      |OPT_TREEDISTRIBUTIONS))
	  un = g->findoptimaledge(s);
	
	if (genopt & OPT_RECMINROOTING)
	  {
	    un->printrooting(cout,1);
	  }

	if (genopt & OPT_RECMINROOTINGEXT)
	  {
	    un->printrootingext(cout,s->root(),1);
	  }
		
	if (genopt & OPT_RECMINCOST)
	  cout << un->cost(s) << endl;
	
	if (genopt & OPT_RECTREECOSTDETAILS) {
	  if (un->p())
	    un->p()->mark(2 | 8);
	  un->mark(2);
	  g->pf(cout, s);
	  un->pcosts(cout,un->cost(s),s);
	}
	
	if ((genopt & OPT_SUMMARYTOTAL)
	    || (genopt & OPT_SUMMARYDLTOTAL)) {
	  DlCost s1 = un->cost(s);
	  total.loss += s1.loss;
	  total.dup += s1.dup;
	}
	
	if ((genopt & OPT_SUMMARYDISTRIBUTIONS) || (genopt
						    & OPT_TREEDISTRIBUTIONS))
	  un->costdet(s);
      } // gt-loop
      
      if (genopt & (OPT_SUMMARYTOTAL | OPT_SUMMARYDLTOTAL
		    | OPT_SUMMARYDISTRIBUTIONS))
	cout << *s << "\t";
      
      if (genopt & OPT_SUMMARYTOTAL)
	cout << total.mut() << "\t";
      if (genopt & OPT_SUMMARYDLTOTAL)
	cout << total << "\t";
      
      if (genopt & (OPT_SUMMARYTOTAL | OPT_SUMMARYDLTOTAL
		    | OPT_SUMMARYDISTRIBUTIONS))
	cout << endl;
      
      if (genopt & OPT_SUMMARYDISTRIBUTIONS)
	s->showcostdet(cout);
      if (genopt & OPT_TREEDISTRIBUTIONS)
	s->pfcostdet(cout);
      
    } // st-loop
  } // (OPT_BYCOST)
  
}
