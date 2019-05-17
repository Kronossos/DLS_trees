/*
 * tools.cpp
 *
 *  Created on: 18-09-2010
 *      Author: gorecki
 */

#include <string.h>

#include <stdlib.h>
#include <iostream>
#include <sstream>
using namespace std;

#include "tools.h"


vector<string> specnames;
vector<int> specorder;
map<string,int> specnames2id;

int rsort=0;   // Important in nni - used to generated unique r|u species trees
int gcenter=0; // Useful for drawing gene trees
int costfromroot=0;
int printleaves=0;
int stprogress=0; // Print progress for species trees
int ppgraphviz=0; // Print graphviz output

void nnit(vector<Tree*> &vec)
{
  vector<Tree*>::iterator stpos;
  for (stpos=vec.begin(); stpos !=vec.end(); ++stpos)
      if ((*stpos)->rooted())
	((SpeciesTree*)(*stpos))->nni(cout);
      else
	((UTree*)(*stpos))->normalize()->nni(cout,1);
}

void nnist(vector<SpeciesTree*> &stvec)
{
	vector<SpeciesTree*> dest;
	vector<SpeciesTree*>::iterator stpos;
	for (stpos=stvec.begin(); stpos !=stvec.end(); ++stpos)
	  {
	    (*stpos)->nni(cout);
	  }
}

void nnigt(vector<UTree*> &gtvec)
{
	vector<UTree*> dest;
	vector<UTree*>::iterator gtpos;
	for (gtpos=gtvec.begin(); gtpos !=gtvec.end(); ++gtpos)
	{
		SpeciesTree *st = (*gtpos)->normalize();
		st->nni(cout,1);
	}
}


//#define BUFSIZE 10000
//void readgtree(char *fn, set<UTree*> &gtset) {
//	FILE *f;
//	f = fopen(fn, "r");
//	if (!f) {
//		cerr << "Cannot open file " << fn << endl;
//		exit(-1);
//	}
//	while (1) {
//		char buf[BUFSIZE];
//		if (!fgets(buf, BUFSIZE, f))
//			break;
//		gtset.insert(new UTree(buf));
//	}
//	fclose(f);
//}


#define BUFSIZE 10000
void readgtree(char *fn, vector<UTree*> &gtvec)
{
	FILE *f;
	f= fopen(fn,"r");
	if (!f)
	{
	  cerr << "Cannot open file " << fn << endl;
	  exit(-1);
	}
	int num=0;
	while (1)
	  {
	    char buf[BUFSIZE];
	    if (!fgets(buf,BUFSIZE,f)) break;
	    num++;
	    UTree *ut = new UTree(buf,num,"",1.0,1);
	    if (ut->leaves()!=1) gtvec.push_back(ut);
	    else cerr << "Warning: a gene tree with one leaf - omitting" << endl;
	}
	fclose(f);
}
//
//void readstree(char *fn, set<SpeciesTree*> &stset) {
//	FILE *f;
//	f = fopen(fn, "r");
//	if (!f) {
//		cerr << "Cannot open file " << fn << endl;
//		exit(-1);
//	}
//	while (1) {
//		char buf[BUFSIZE];
//		if (!fgets(buf, BUFSIZE, f))
//			break;
//		stset.insert(new SpeciesTree(buf));
//	}
//	fclose(f);
//}

void readsmp(istream &c, vector<SpeciesTree*> &vec) 
{ 
  int num=0;
  while(!c.eof()) {
    string t;
    if (t.size()==0) return;
    num++;
    vec.push_back(new SpeciesTree(strdup(t.c_str()),num));
  }
}

Tree *gettree(string t, int rooted, int num, string tn, float wg)
{
  char *s=strdup((t).c_str());
  if (rooted)
    return new SpeciesTree(s,num,tn,wg);
  return new UTree(s,num,tn,wg);
}


void readfields(istream &c, int rooted, int field,vector<string> &prefix, vector<string> &suffix, vector<Tree*> &vec)
{
  string s;
  while (getline(c,s)) 
    {
      vector<string> tokens;
      vector<string>::iterator t;
      tokenize(s, tokens," ");	    
      int cnt=1;
      string pref="";
      string suff="";
      for (t = tokens.begin(); t != tokens.end(); ++t) 
	{
	  if (cnt<field) pref=pref+" "+*t;
	  if (cnt>field) suff=suff+" "+*t;
	  if (cnt==field)
	    vec.push_back(gettree(*t,rooted,vec.size(),"",1.0));
	  cnt++;
	}
      prefix.push_back(pref);
      suffix.push_back(suff);
    }
}
  


void readt(istream &c, int treenames, int weights, int rooted, vector<Tree*> &vec) 
{ 
  int num=0;
  while(!c.eof()) {
    string l,tn;
    double wg;
    if (treenames) c >> tn; 
    if (weights) c >> wg;
    if (!getline(c,l)) break;       
    num++;
    vec.push_back(gettree(l,rooted,num,tn,wg));
  }
}



void readstree(char *fn,vector<SpeciesTree*> &stvec)
{
	FILE *f;
	f= fopen(fn,"r");
	if (!f)
	{
		cerr << "Cannot open file " << fn << endl;
		exit(-1);
	}
	int num=0;
	while (1)
	{
		char buf[BUFSIZE];
		if (!fgets(buf,BUFSIZE,f)) break;
		num++;
		stvec.push_back(new SpeciesTree(buf,num));
	}
	fclose(f);
}

void readtreesinterleaved(char *fn,vector<UTree*> &gtvec, vector<SpeciesTree*> &stvec)
{
	FILE *f;
	f= fopen(fn,"r");
	if (!f)
	{
		cerr << "Cannot open file " << fn << endl;
		exit(-1);
	}
	int num=0;
	while (1)
	{
		char buf[BUFSIZE];
		if (!fgets(buf,BUFSIZE,f)) break;
		if (num%2)
		  stvec.push_back(new SpeciesTree(buf,num/2));
		else
		  {
		    UTree *ut = new UTree(buf,num/2,"",1.0,1);
		    if (ut->leaves()!=1) gtvec.push_back(ut);
		    else cerr << "Warning: a gene tree with one leaf - omitting" << endl;
		  }
		num++;
	}
	fclose(f);
}

vector<string> initgenrand(char *optarg) {
  
	int m=0;
	int i;
	vector<string> ss;
	int num;
	if (sscanf(optarg,"%d",&num)==1)
	  {
	    char buf[10];
	    int md ='z'-'a'+1;
	    for (i=0;i<num;i++)
	      {
		char ll='a'+(i%md);
		int r=i/md;
		if (!r) sprintf(buf,"%c",ll);
		else sprintf(buf,"%c%d",ll,r);
		getspecies(strdup(buf));
		ss.push_back(strdup(buf));
	      }
	    return ss;
	    
	  }
	for (i=0;i<(int)strlen(optarg);i++)
		if (strchr(",() ",optarg[i]))
		{
			m=1;
			break;
		}

	if (m)
	{
		char *f=strtok(optarg,"(), ");
		getspecies(f);
		ss.push_back(f);
		while ((f=strtok(NULL,"(), ")))
		{
			getspecies(f);
			ss.push_back(f);
		}
	}
	else
	{
		char buf[2];
		buf[1]=0;
		for (i=0;i<(int)strlen(optarg);i++)
		{
			buf[0]=optarg[i];
			ss.push_back(buf);
			getspecies(strdup(string(buf).c_str()));
		}
	}
	return ss;

}

void sortmode(char *optarg)
{
	int j;
	for (j=0;j<(int)strlen(optarg);j++)
	{
		if (optarg[j]=='c') gcenter=1;
		if (optarg[j]=='v') gcenter=1;
		if (optarg[j]=='s') rsort=1;
		if (optarg[j]=='g') ppgraphviz=1;
		if (optarg[j]=='l') printleaves=1;
		if (optarg[j]=='S') {
		  int i;
		  for (i=0; i<(int)specnames.size(); i++)
		    cout << i << ". " << specnames[i] << endl;
		}
	}
}

void setalordering()
{
  //  specorder[num]->num
  // ugly 
  int i,j;
  for (i=0;i<(int)specnames.size();i++)
    {
      int cnt=0;
      for (j=0;j<(int)specnames.size();j++)
	if (specnames[i]>specnames[j]) cnt++;
      specorder[i]=cnt;
      //      cout << i << " " << cnt << endl;
    }
  
}


int _addnode(ostringstream &f,UNode *n,double x,double y,double &a,
double &b, double &c, double &d)
{
  static int nodeid=1;
  nodeid++;
  if (n->leaf())
    f << "\\rput[c](" << x << "," << y << "){\\ovalnode{A" 
      << nodeid << "}{" << specnames[((ULeaf*)n)->species()] << "}}" << endl;
  else
    f << "\\cnode*(" << x << "," << y << "){2pt}{A" << nodeid << "}" << endl;
  a=min(a,x);
  b=min(b,y);
  c=max(c,x);
  d=max(d,y); 
  return nodeid;
}

void _addedge(ostringstream &f,UNode *n,int s, int e, int withlen=0)
{
  if (withlen)
    {
    f << "\\ncarc[nodesep=1pt,ncurv=0.9]{->}{A" << s << "}{A" << e << "}" << endl;
    f << "\\aput[1pt]{:U}{" << n->len() << "}" << endl;
    //\aput{:U}{Hypotenuse}
    }
  else
    f << "\\ncarc[nodesep=0pt]{A" << s << "}{A" << e << "}" << endl;

}

#define SMALLDELTA 0.2
#define DECRFACTOR 0.8
#define EDGELEN 1.0
#define YFACTOR 0.5

void _gentree(ostringstream &f,UNode *n, double x, double y, int nid, double coef, 
	      double depth,double &a,double &b, double &c, double &d)
{
  if (n->leaf()) return; 
  UNode3 *u = (UNode3*)n;
  
  double xpos=x+SMALLDELTA*coef;
  double ypost=y+SMALLDELTA*YFACTOR;
  double yposb=y-SMALLDELTA*YFACTOR;

  int tn=_addnode(f,u->l(),xpos,ypost,a,b,c,d);
  int bn=_addnode(f,u->r(),xpos,yposb,a,b,c,d);
  _addedge(f,NULL,nid,tn);
  _addedge(f,NULL,nid,bn);
  _addedge(f,NULL,bn,tn);

  double xpos2=xpos+EDGELEN*coef;
  double ypos2t=ypost+depth*EDGELEN*YFACTOR;
  double ypos2b=yposb-depth*EDGELEN*YFACTOR;

  UNode* tpar=u->l()->p();
  UNode* bpar=u->r()->p();

  int tnp=_addnode(f,tpar,xpos2,ypos2t,a,b,c,d);	     
  int bnp=_addnode(f,bpar,xpos2,ypos2b,a,b,c,d);

  _addedge(f,u->l(),tn,tnp,1);
  _addedge(f,tpar,tnp,tn,1);

  _addedge(f,u->r(),bn,bnp,1);
  _addedge(f,bpar,bnp,bn,1);

  double decr=depth; 
  if (!tpar->leaf() && !bpar->leaf()) decr*=DECRFACTOR;

  _gentree(f,tpar,xpos2,ypos2t,tnp,coef,decr,a,b,c,d);
  _gentree(f,bpar,xpos2,ypos2b,bnp,coef,decr,a,b,c,d);
  
}

/*\cnode*(0,0){3pt}{A}
\cnode*(4,2){3pt}{B} \ncline[nodesep=3pt]{A}{B}
\mput*{1} */

string printpstricks(UTree *t,UNode *n)
{
  ostringstream f;
  
  if (!n) n=t->start;

  double a=0,b=0,c=0,d=0;

  //  f << "\\begin{pspicture}(-3,-4)(3,4)" << endl;  // @TODO dimensions
  int sn= _addnode(f,n,0,0,a,b,c,d);
  int en= _addnode(f,n->p(),EDGELEN,0,a,b,c,d);

  _addedge(f,n,sn,en,1);
  _addedge(f,n->p(),en,sn,1);

  _gentree(f,n,     0,0,sn,-1.0,2,a,b,c,d);
  _gentree(f,n->p(),EDGELEN,0,en,1.0,2,a,b,c,d);
  
  f << "\\end{pspicture}" << endl;

  ostringstream g;

  g << "\\begin{pspicture}(" << a << "," << b << ")(" << c << "," << d << ")" << endl;
  return g.str() + f.str();

}

bool eqSpeciesTrees(SpeciesTree *t1, SpeciesTree *t2)
{
  return t1->eq(t2);
}



void tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters = " ")
{
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}
