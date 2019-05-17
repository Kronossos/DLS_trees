/************************************************************************
 Unrooted REConciliation version 1.00
 (c) Copyright 2005-2006 by Pawel Gorecki
 Written by P.Gorecki.
 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/
#include <string.h>
#include "tools.h"

#include <iostream>
using namespace std;

#define OPT_NNI 1
#define OPT_PRINTSPNAMES 2
#define OPT_PRINTALNAMES 8
#define OPT_PRINTROOTINGS 4
//#define OPT_NNI (1<<25)
//#define OPT_RSORT (1<<26)


int usage(int argc, char **argv) {
  cout << 
    " TreeMan - tree manipulations\n"
    " Read from a file/string or stdin and transform a tree\n"
    " Usage: " << argv[0] << " [options] files\n"
    " -t u|r - type of input trees\n"
    " -s leaves|tree - defines species names\n"
    " -o CMD - where CMD is one of the following\n" 
    "    r - show rootings for every unrooted tree\n"
    "    n - show nni's for every rooted/unrooted tree\n"
    "    E - extended notation (e.g., nni's with comments)\n"	
    "    c - show center variant of unrooted (check it)\n"
    "    s - print alphanum sorted trees (unrooted - can be incorrect)\n"
    "    A - rename a,b,c,... labels"
    "    g - ppgraphviz output (unrooted only)\n"
    "    S - print species names\n"	
    "    u - print unique rooted trees\n"
    "    x - simple reading\n"
    "    e - read pairs, print 1 if rooted trees are equal\n"
    "        otherwise 0\n"
    " -c defs_withtreenames - replace trees by names\n"
    " -f NUM - column number where tree is defined (use with -c)\n"
    " -d NUM - compute difference between first NUM trees and the remaining ones\n"
    "        print its treename otherwise print whole tree\n"
    " -x - files only; dont convert to trees in case of err file\n\n"
    " Examples:"
    "	  ./treeman -f2 -tu -og ../remrun/g3bl.known.txt | dot -Tps  > g.ps"

    "";


	exit(-1);
}


int treenames=0;
int weights=0;

int diff=-1;


int main(int argc, char **argv) {
  int opt;
  
  vector<Tree*> vec;
  vector<string> prefix;
  vector<string> suffix;

  vector<UTree*> gtvec;
  srand(time(0));
  
  int genopt = 0;
  int rooted=1;

  int eqtest=0;
  char *convert=NULL;
  int fileonly=0; // if error reading a file treat the string as tree definition
  int smp=0;
  int field=0;
  int unique=0;
  int extended=0;

  if (argc<2) usage(argc,argv);
  
  while ((opt = getopt(argc, argv, "t:o:s:h:c:xf:d:"))!=-1)
    {
    switch (opt) 
      {
      case 't':
	rooted=(optarg[0]!='u');
	break;
	
      case 'o':
	int k;
	for(k=0;k<(int)strlen(optarg);k++)
	  {
	    switch(optarg[k])
	      {
	      case 'n': genopt |= OPT_NNI; break;
	      case 'r': genopt |= OPT_PRINTROOTINGS; break;
	      case 'c': gcenter=1; break;
	      case 's': rsort=1; break;
	      case 'u': unique=1; break;
	      case 'E': extended=1; break;
	      case 'x': smp=1; break;
	      case 'g': ppgraphviz=1; break;
	      case 'S': genopt |= OPT_PRINTSPNAMES; break;
	      case 't': treenames=1; break;  // first column
	      case 'w': weights=1; break; // first or second (if treename is defined)
	      case 'e': eqtest=1; break; // check if trees are equal
	      default:
		cerr << "Undefined variant in -o option";
		exit(-1);
	      }
	  }
	break;

      case 'd':
	diff=atoi(optarg);
	break;
	
      case 'c':
	convert=strdup(optarg);
	break;

      case 'f':
	field=atoi(optarg);
	break;

      case 'x':
	fileonly=1;
	break;
	      
      case 's':
	initgenrand(optarg);
	break;
	      
      default:
	cerr << "Unknown option: " << ((char) opt) << endl;
	exit(-1);
      }
    }
	
  vector<Tree*>::iterator stpos, stpos2;
  
  int pr=0;
  int rd=0;
  vector<Tree*> st2; 

  if (convert)
    {
      ifstream c(convert);
      readt(c, 1, 0, rooted, st2);
      c.close();
      
      if (st2.size()>1)
	{
	  unsigned int i,j;
	  
	  for (i=0;i<st2.size()-1;i++)
	    for (j=i+1;j<st2.size();j++)
	      if (st2[i]->eq(st2[j]))
		cerr << "Warning: candidates " << st2[i]->treename << " " << st2[j]->treename << " are equal" << endl;
	}

    }
  
  // read input files
  int index;
  for (index = optind; index < argc; index++)
    {
      ifstream c(argv[index]);
      if (field)
	readfields(c,rooted,field,prefix,suffix,vec);	  
      else readt(c,treenames,weights,rooted,vec);
      c.close();
      rd=1;
    }
  if (!rd) 
    {
      if (field) 
	readfields(cin,rooted,field,prefix,suffix,vec);
      else readt(cin,treenames,weights,rooted,vec);
    }
  rd=1;

  // convert
  if (convert)
    {
      int i,j;
      for (i=0;i<(int)vec.size();i++)
	{
	  if (field) cout << prefix[i] << " ";
	  int fnd=0;
	  for (j=0;j<(int)st2.size();j++)
	    if (vec[i]->eq(st2[j]))
	      {
		cout << st2[j]->treename << " ";
		fnd=1;
		break;
	      }	    
	  if (!fnd) 
	    (vec[i])->print(cout) << " "; 
	  if (field) cout << suffix[i] << " ";
	  cout << endl;
	}
      pr=1;      
    }
  
  
  // gen unique
  if (!pr && unique)
    {
      int i,j;
      for (i=0;i<(int)vec.size();i++)
	{
	  int fnd=0;
	  for (j=i+1;j<(int)vec.size();j++)
	    if (vec[i]->eq(vec[j]))
	      {
		fnd=1;
		break;
	      }
	  if (!fnd)
	    {
	      if (field) cout << prefix[i] << " ";
	      (vec[i])->print(cout);
	      if (field) cout << suffix[i] << " ";
	      cout << endl;	      
	    }

	}      
      pr=1;
    }

  // comp diff.
  if (!pr && diff>0)
    {
      int i,j;
      if (diff>(int)vec.size()) diff=vec.size();
      for (i=0;i<diff;i++)
	{
	  int fnd=0;
	  for (j=diff;j<(int)vec.size();j++)
	    if (vec[i]->eq(vec[j]))
	      fnd=1;
	  if (!fnd)
	    {
	      if (field) cout << prefix[i] << " ";
	      (vec[i])->print(cout);
	      if (field) cout << suffix[i] << " ";
	      cout << endl;	      
	    }
	}
      pr=1;
    }

  if (!pr && convert)
    {
      for (stpos = vec.begin(); stpos != vec.end(); ++stpos)
	{
	  int fnd=0;
	  for (stpos2 = st2.begin(); stpos2 != st2.end(); ++stpos2)	
	    if ((*stpos)->eq(*stpos2))
	      {
		cout << (*stpos2)->treename << endl;
		fnd=1;
		break;
	      }	    	    
	  if (!fnd) (*stpos)->print(cout) << endl;
	}
      
      pr=1;
    }



  if (genopt & OPT_PRINTSPNAMES)
    {
      int i;
      for (i=0; i<(int)specnames.size(); i++)
	cout << specnames[i] << endl;
      pr=1;
    }
   
  if (genopt & OPT_NNI)
    {
      // NNI operations for the set of gene/species trees      
      int i;
      for (i=0;i<(int)vec.size();i++)
	{	  
	  if (extended)
	    {
	      if (field) cout << "#" << prefix[i] << endl;
	      cout << "#";
	      (vec[i])->print(cout) << endl;
	      if (field) cout << "#"<< suffix[i] << endl;
	    }
	  if ((vec[i])->rooted())
	    ((SpeciesTree*)(vec[i]))->nni(cout);
	  else
	    ((UTree*)(vec[i]))->normalize()->nni(cout,1);
	}

      pr=1;
    }
  
  
  if (genopt & OPT_PRINTROOTINGS) {
    
    for (stpos = vec.begin(); stpos != vec.end(); ++stpos)
      if ((*stpos)->rooted())
	(*stpos)->print(cout) << endl;
      else
	((UTree*)(*stpos))->printrootings(cout);
    pr=1;
  }

  if (eqtest)
    {
      if (vec.size()<2)
	{
	  cerr  << "Pairs of trees expected..." << endl;
	  exit(-1);
	}
      if (vec[0]->eq(vec[1])) cout << 1 << endl;
      else cout << 0 << endl;
      pr=1;
    }

  if (rsort) setalordering();
  
  if (!pr && (gcenter || rsort || ppgraphviz))
    {
      int i;
      for (i=0;i<(int)vec.size();i++) 
	{
	  
	  if (field && !ppgraphviz) cout << prefix[i] << " ";

	  if (!vec[i]->rooted())
	    {
	      UTree *t=(UTree*)vec[i];
	      
	      if (gcenter) t->center();
	      if (rsort) t->normalize()->printsorted(cout,0);
	      else if (ppgraphviz)
		t->printgraphviz(cout);
	      else t->print(cout);
	    }
	  else
	    {
	      SpeciesTree *s=(SpeciesTree*)vec[i];
	      if (rsort)
		s->printsorted(cout);
	      else
		s->print(cout);
	    }

	  if (field && !ppgraphviz) cout << suffix[i] << " ";
	  cout << endl;	      

	}
      pr=1;
    }
  
}

