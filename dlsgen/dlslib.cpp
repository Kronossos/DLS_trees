



#include "dlslib.h"

costtype **costmatrix=NULL;
int usecomma=0;
int costcomp=GENEDUPLICATIONS;
int swapping=0;
int opt2=1;
int minsearch=0;
int alloweqcosts=0;
int stnum;
int gtnum;

int treereprsize(int n) { return 4*(n-1)+2; }

//#define _TRIPLETDEBUG_

//#define _SWAPHEUR_



char sortNode(Node *n)
{
  if (n->label) return n->label;  
  char a=sortNode(n->a);
  char b=sortNode(n->b); 
  if (a>b)
  {
    Node *z=n->a;
    n->a=n->b;
    n->b=z;
    return b;
  }
  return a;
}

// Trees are sorted
int eqNodes(Node *n,Node *m)
{
  if (n->label!=m->label) return 0; 
  // 0=0 internal or lf=lf
  if (n->label) return 1; // leaves
  // internal
  return eqNodes(n->a,m->a) && eqNodes(n->b,m->b);
}





void deleteTree(Node *n)
{
  if (!n->label)    
    {
      deleteTree(n->a);
      deleteTree(n->b);      
    }
  free(n);
}


void setdepth(Node *n, int st)
{
  n->depth=st;
  if (!n->label)  
    {
      setdepth(n->a,st+1);
      setdepth(n->b,st+1);
    }
}


void ppTree(Node *n, int e)
{
  if (n->label)
    printf("%c",n->label);
  else
    {
      printf("(");
      ppTree(n->a,0);
      printf(",");
      ppTree(n->b,0);
      printf(")");
    }
  if (e) printf("\n");
}


void ppscen(Node *n, int e=1)
{
  if (n->label)
    printf("%c",n->label);
  else
    {
      printf("(");
      ppscen(n->a,0);      
      printf(",");
      ppscen(n->b,0);
      printf(")");
      int cnt=0;
      Node *stm=n->map;
      while (stm!=n->scenmap) {
        stm=stm->p;
        cnt+=1;
      }
      if (n->lcaspec) printf("~");
      if (!n->scenspec) printf("+%d",cnt);
    }
  if (e) printf("\n");
}


Node *sibling(Node *n)
{
  if (!n->p) return NULL;
  if (n->p->a==n) return n->p->b;
  return n->p->a;  
}

void ppdls(Node *n, Node *pp, int e)
{
  
    Node *s=n->scenmap;
    Node *s1=s,*s2=s;      
    if (pp)
    {
      s=n->scenmap;
      while (s!=pp)
      {
        printf("(");    
        s=s->p;
      }
    }

    if (n->label) printf("%c",n->label);
    else  
    {
      s=n->scenmap;      
      if (n->scenspec)
      {
        s1=n->a->scenmap;
        // ppTree(s1);
        // ppTree(s);

        while (s1->p!=s) s1=s1->p;
        s2=n->b->scenmap;
        while (s2->p!=s) s2=s2->p;        
      }      
      printf("(");
      ppdls(n->a,s1,0);      
      printf(",");
      ppdls(n->b,s2,0);               
      printf(")");    
      if (n->scenspec) printf("~"); // auto anyway
      else printf("+");
    }
      // add losses to the pp
  if (pp)
  {
    char buf[1000];
    s=n->scenmap;
    while (s!=pp)
    {
      printf(",-");
      treecluster(sibling(s),buf," ");
      printf("%s)",buf);
      s=s->p;
    }
  }
    
  if (e) printf("\n");
}



int scenlevel(Node *n)
{
  if (n->label) return 0;

  int l=scenlevel(n->a)+scenlevel(n->b);
  Node *stm=n->map;
  while (stm!=n->scenmap) {
    stm=stm->p;
    l+=1;
  }
  if (n->lcaspec && !n->scenspec) l+=1;   // spec->dup 
  return l;
}





int treeleafnum(Node *s)
{
   if (s->label) return 1;
   return treeleafnum(s->a)+treeleafnum(s->b);
}

void treecluster(Node *s, char *buf, const char *sep)
{
  if (s->label) {
      sprintf(buf,"%c%s",s->label,sep);
      return;
  }
  treecluster(s->a,buf,sep);
  treecluster(s->b,buf+strlen(buf),sep);
}




Node *lca(Node *a, Node *b, LCA* l)
{
    if (a==b) return a;
    if (l->lcat[a->id][b->id]) return l->lcat[a->id][b->id];
    if (a->depth>b->depth) return l->lcat[a->id][b->id]=lca(a->p,b,l);
    if (a->depth<b->depth) return l->lcat[a->id][b->id]=lca(a,b->p,l);
    return l->lcat[a->id][b->id]=lca(a->p,b->p,l);
}

void insertNodes(Node *s,Node **t,int *lfp, int *itp, int depth)
{
  s->depth=depth;
  if (s->label)
  {
    t[*lfp]=s;
    s->lfid=s->id=*lfp;
    *lfp+=1;
    return;
  }
  insertNodes(s->a,t,lfp,itp,depth+1);
  insertNodes(s->b,t,lfp,itp,depth+1);
  
  t[*itp]=s;
  s->id=*itp;
  *itp+=1;  
  if (s->p)
  {
    if (s==s->p->a) s->lfid=s->a->lfid;
    else s->lfid=s->b->lfid;
  } 
}



LCA* genLCA(Node *s)
{
  LCA *l=(LCA*)malloc(sizeof(LCA));
  l->n=treeleafnum(s);
  l->root=s;
  l->tfullsize=2*l->n-1;
  l->t=(Node**)malloc(l->tfullsize*sizeof(Node*));
  int lfp=0,itp=l->n;
  insertNodes(s,l->t,&lfp,&itp,0);
 
  // prep LCA
  l->lcat=(Node***)malloc(l->tfullsize*sizeof(Node**));
  int i,j;
  for(i=0;i<l->tfullsize;i++)
  {
    l->lcat[i]=(Node**)malloc(l->tfullsize*sizeof(Node*));
    for(j=0;j<l->tfullsize;j++) l->lcat[i][j]=NULL;    
  }

  // fill LCA
  for(i=0;i<l->tfullsize;i++)
  {      
    l->lcat[i][i]=l->t[i];  
    for(j=0;j<l->tfullsize;j++)
      if (i!=j)
        l->lcat[i][j]=lca(l->t[i],l->t[j],l);
  }
  

  // a2n
  for(i=0;i<l->n;i++) l->a2n[l->t[i]->label]=l->t[i];

  return l;
}

Node *gtree(Node *a, Node *b)
{
  Node *n = (Node*)malloc(sizeof(Node));
  n->a=a;
  n->b=b;
  n->label=0;
  n->a->p=n;
  n->b->p=n;
  n->a->s=n->b;
  n->b->s=n->a;  
  n->p=NULL;
  return n;
}

Node *gleaf(char a)
{
  Node *n = (Node*)malloc(sizeof(Node));
  n->label=a;  
  return n;
}

Node *_parseTree(char *s,int *p)
{
  if (s[*p]=='(')
    {
      Node *a; 
      Node *b;             
      *p+=1;
      a = _parseTree(s,p);    
      if (s[*p]!=',') 
      {
        fprintf(stderr,"Error - comma expected. Found <%s> in <%s> at pos %d.",s+*p,s,*p);      
        exit(-1);
      }
      *p+=1;
      b = _parseTree(s,p);
      if (s[*p]!=')') 
      {
        fprintf(stderr,"Error - ) expected. Found <%s> in <%s> at pos %d.",s+*p,s,*p);      
        exit(-1);
      }
      *p+=1;
      return gtree(a,b);
    }
  *p+=1;
  return gleaf(s[*p-1]);
}



Node *parseTree(char *s)
{
  int p=0;
  Node *r=_parseTree(s,&p);
  setdepth(r,0);
  return r;
}





int eqshapes(Node *a,Node *b)
{
  if (a->label && b->label) return 1;
  if (a->label || b->label) return 0;
  return (eqshapes(a->a,b->a) && eqshapes(a->b,b->b)) || (eqshapes(a->b,b->a) && eqshapes(a->a,b->b));
}



#define getlca(a,b) s->lcat[a->id][b->id]
Node* map(LCA *s, Node *g)
{
  if (g->label) 
    g->map=s->a2n[g->label];    
  else
  {
    Node *a=map(s,g->a);
    Node *b=map(s,g->b);  
    g->map=getlca(a,b); 
    g->lcaspec=(g->map!=a && g->map!=b);  
  }
  return g->map;
}







void printcluster(Node *s, int n)
{
  char  buf[n+1];
  treecluster(s,buf,"");
  printf("%s",buf);
}

int subseteq(char *gc, char *sc)
{
  while (*gc)
    if (!strchr(sc,*gc++)) return 0;
  return 1;
}



void printclusterclusterstats(LCA *s, Node *gt, char *opt)
{
  LCA *g=genLCA(gt);
  int i,j;
  int n=treeleafnum(gt);
  char gc[n+1],sc[n+1];
  int inc[s->tfullsize][g->tfullsize];
  printf("S="); 
  ppTree(s->root);
  printf("\nG\\S");
  for (i=0;i<s->tfullsize;i++) { 
    printf("\t");
    printcluster(s->t[i],n);
  }
  printf("\n");
  for (j=0;j<g->tfullsize;j++) { 
    treecluster(g->t[j],gc);
    printf("%s",gc);
    int rowsum=0;
    for (i=0;i<s->tfullsize;i++)
      {
        treecluster(s->t[i],sc);
        int res=subseteq(gc,sc);
        printf("\t%d",res);
        inc[i][j]=res;
        rowsum+=res;
      }
    printf("\t%d\n",rowsum);
  }
  int total=0;
  for (i=0;i<s->tfullsize;i++)
  {
    int colsum=0;
    for (j=0;j<g->tfullsize;j++) 
      colsum+=inc[i][j];
    printf("\t%d",colsum);
    total+=colsum;
  }
  printf("\t%d\n",total);  
}


Node *randStTree(int n)
{
  
  Node* t[n];
  for (int i=0; i<n; i++) t[i]=gleaf(i+'a');  
    
  while (n>1) {
    int a=rand()%n;
    int b=rand()%n;
    
    if (a==b) continue;
    Node *nx=gtree(t[a],t[b]);    
    if (a<b) { t[a]=nx; }
    else { t[b]=nx; b=a; }
    for (int j=b;j<n-1;j++) t[j]=t[j+1];
    n--; 
  }
  return t[0];
}



void getInternals(Node *st, Node **in, int *p)
{
  if (st->label) return;
  in[(*p)++]=st;
  getInternals(st->a,in,p);
  getInternals(st->b,in,p);
}





long int _scenCount(Node *n, Node *parmap, bool parspec, bool fixedspec)
{
    if (n->label) return 1;

    bool nspec=n->lcaspec;
    Node *nmap=n->map;
    long cnt=0;
    while (1)
    {
      cnt+=_scenCount(n->a,nmap,nspec,fixedspec)*_scenCount(n->b,nmap,nspec,fixedspec);

      if (!parspec && !nspec && nmap==parmap) break;

      if (nspec)
        { 
          if (fixedspec) break; 
          else nspec=0;
        } 
      else
        nmap=nmap->p;
      if (parspec && nmap==parmap) break;       
    }
    // printf(" %ld parspec=%d parmap=",cnt,parspec);
    // ppTree(parmap,0);
    // printf(" n=");
    // ppTree(n);
    
    return cnt;
}

long int scenCount(Node *gt, Node *stroot, bool fixedspec)
{
  return _scenCount(gt,stroot,0,fixedspec);
}


long int _genscen(long id, Node *n, Node *parmap, bool parspec, bool fixedspec)
{
    if (n->label) { 
      n->scenmap=n->map;
      return 1;
    }
    // ppTree(n,0);
    // printf(" %ld\n",id);

    bool nspec=n->lcaspec;
    Node *nmap=n->map;
    long cnt=0;
    while (1)
    {
      long ccnta=_scenCount(n->a,nmap,nspec,fixedspec);
      long ccntb=_scenCount(n->b,nmap,nspec,fixedspec);
      // printf("FIX id=%ld ac=%ld bc=%ld \n",id,ccnta,ccntb);
        
      if (ccnta*ccntb>id)
      {
        n->scenmap=nmap;
        n->scenspec=nspec;
        _genscen(id%ccnta,n->a,nmap,nspec,fixedspec);
        _genscen(id/ccnta,n->b,nmap,nspec,fixedspec);              
        return 1;
      }
      id=id-ccnta*ccntb;

      if (!parspec && !nspec && nmap==parmap) break;

      if (nspec)
        { 
          if (fixedspec) break; 
          else nspec=0;
        } 
      else
        nmap=nmap->p;
      if (parspec && nmap==parmap) break;       
    }
    // printf(" %ld parspec=%d parmap=",cnt,parspec);
    // ppTree(parmap,0);
    // printf(" n=");
    // ppTree(n);
    
    return cnt;
}



void travme(Node *gt,int *epi, int dupup)
{   
    Node *parmap=NULL;
    if (gt->p) parmap=gt->p->scenmap;
    if (gt->label || gt->scenspec)
    {
      // leaf or speciation
      // store dupup
      if (dupup && dupup>epi[parmap->id]) epi[parmap->id]=dupup;
      dupup=0;        
    }
    else
    {
      // internal and dupl.
      // store if map are different
      if (parmap!=gt->scenmap)
        { 
          if (parmap && dupup>epi[parmap->id]) epi[parmap->id]=dupup;
          dupup=0;
        }      
      dupup+=1;    
    }

    if (gt->label) return;
    travme(gt->a,epi,dupup);
    travme(gt->b,epi,dupup);
}


void travrule(Node *n,int *rule)
{   
    int l=0;
    Node *stm=n->map;
    while (stm!=n->scenmap) {
      stm=stm->p;
      l+=1;
    }
    if (n->lcaspec && !n->scenspec) { l+=1; l=-l; }   // spec->dup 
    rule[n->id]=l;
    //printf("%d ",n->id);
    if (n->label) return;
    travrule(n->a,rule);
    travrule(n->b,rule);    
}



int* calcepi(Node *gt, LCA *s, int *epi)
{
  if (!epi) epi=(int*)malloc(sizeof(int)*(s->tfullsize)); 
  int i;
  for (i=0;i<s->tfullsize;i++) epi[i]=0;
  travme(gt,epi,0);
  return epi;  
}


int* calcrule(Node *gt, int gsize, int *rule)
{
  if (!rule) rule=(int*)malloc(sizeof(int)*(2*gsize+1)); 
  int i;
  for (i=0;i<2*gsize+1;i++) rule[i]=0;
  travrule(gt,rule);
  return rule;  
}

int mescore(int *epi, int tsize)
{
  int mes=0;
  int i;
  for (i=0;i<2*tsize-1;i++) mes+=epi[i];
  return mes;   
}


int* calcmprofile(int *epi, LCA *s, int tsize, int *pprofile)
{  
  int i;
  if (!pprofile) pprofile=(int*)malloc(sizeof(int)*tsize); 
  for (i=0;i<tsize;i++) 
  {
    Node *sn=s->t[i];
    int pps=0;
    while (sn)
    {
      pps+=epi[sn->id];
      sn=sn->p;
    }
    pprofile[i]=pps;
  }
  return pprofile;   
}

int mpscore(int *pprofile, int tsize)
{
  int mps=0;
  int i;
  for (i=0;i<tsize;i++) if (pprofile[i]>mps) mps=pprofile[i];
  return mps;
}

void ppvect(int *t,int s)
{
  int i;
  printf("[");
  for (i=0;i<s;i++)
  {
      if (i>0) printf(",");
      printf("%d",t[i]);
  }
  printf("]");
}
// set args "(((((a,b),(c,(d,e))),(f,(g,h))),i),j)" "(((((a,(b,c)),((d,e),(f,g))),(h,i)),((j,(k,l)),(((m,n),o),p))),(q,r))" i

    // int *epi=scentab[i].epi=calcepi(gt,s);
    // int *mprofile=scentab[i].mprofile=calcmprofile(epi,s,s->tsize,NULL);
    // mes=scentab[i].me=mescore(epi, s->tsize);
    // mps=scentab[i].mp=mpscore(mprofile, s->tsize); 
    // scentab[i].mprofile=calcrule(gt,gtsize*2+1);

Scen* genscen(Scen *scen, long id, Node *gt, Node *stroot, LCA  *s, 
  bool fixedspec, int stsize, int gtsize, bool maponly)
{
  _genscen(id,gt,stroot,0,fixedspec);
  if (maponly) return NULL;

  if (!scen) scen=(Scen*)malloc(sizeof(Scen));
  scen->epi=calcepi(gt,s);
  scen->me=mescore(scen->epi,stsize);
  scen->mprofile=calcmprofile(scen->epi,s,stsize,NULL);  
  scen->mp=mpscore(scen->mprofile, stsize); 
  scen->rule=calcrule(gt,gtsize);
  scen->id=id;
  scen->level=scenlevel(gt);
  return scen;

}


int Scen::printinfo(int minmes, int minmps,  int stsize, int gtsize, Node *gt, const char *opts)
{
    int printlabel=1;
    if (strchr(opts,'M'))
    {
        if (this->me!=minmes && this->mp!=minmps) return 0; // skip
    }
    if (strchr(opts,'v')) printlabel=0;
    int i=0;
    for (i=0;i<strlen(opts);i++)
    {
        char o=opts[i];
        switch(o)
        {
          case 'i': 
            if (printlabel) printf("id=");
            printf("%ld ",this->id);        
            break;

          case 'M':        
            if (printlabel) printf("m=");            
            if (this->me==minmes) printf("E");
            if (this->mp==minmps) printf("P");
            printf(" ");
            break;
          case 'o':
            ppscen(gt,0); printf(" "); 
            break;                      
          case 'l': 
             if (printlabel) printf("level=");
             printf("%d ",this->level); 
             break;
          case 'e': 
              if (printlabel) printf("me=");
               printf("%d ",this->me);           

             break;
          case 'p': 
             if (printlabel) printf("mp=");
             printf("%d ",this->mp); 
             break;
          case 'P': 
          if (printlabel) printf("P=");
             ppvect(this->mprofile, stsize); printf(" ");  
             break;
          case 'E': 
             if (printlabel) printf("E=");
             ppvect(this->epi, stsize*2-1); printf(" "); 
             break;
          case 'R': 
             if (printlabel) printf("R=");
             ppvect(this->rule, gtsize*2-1); printf(" ");                          
             break;
        }
      }
    return 1;
}


// a lower (closer to lca), b upper
int ruletransf(Scen *a, Scen *b, int gtsize) 
{
  //printf("Test %ld %ld %d %d\n",a->id,b->id,a->level,b->level);
  if (a->level+1!=b->level) return 0;
  int i;
  int ids=2*gtsize-2;
  int diff=0;
  for (i=0;i<gtsize*2-1;i++)
  {
    if (a->rule[i]==b->rule[i]) ids--;
    else diff=i;
  }
  //printf("CCC %d %d",diff,ids);
  if (ids) return 0;
  if (abs(a->rule[diff])+1!=abs(b->rule[diff])) return 0; // no edge
  if (a->rule[diff]==0 and b->rule[diff]==-1) return -1; // spec conv.
  return 1; // dup conv.
}



