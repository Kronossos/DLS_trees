
#include "dlslib.h"

#include <stdio.h>  
#include <unistd.h> 

int levellimit=100000;


#define GENALL 1
#define GENFIXEDSPEC 2
#define GENMP 4
#define GENME 8

int usage(int argc, char **argv) {
  printf(
  " DLS generator v1.0. (C) P.Gorecki 2019\n"
  " Usage: %s [options] GeneTree SpeciesTree\n",argv[0]);
  printf(
    " -g ascC - gen scenarios\n"
    "    a - all, default\n"
    "    s - only from PG model; no TMOVE-1 rules\n"    
    "    c - print scenario count\n"    
    "    C - print scenario count and exit\n"    
    " -v INT - verbose level (def. 0)\n"
    " -d aep - gen diagram graph\n"
    "    a - all scenarios\n"
    "    e - minme scenarios\n"
    "    p - minmp scenarios \n"
    "    ep - minmp+minme scenarios \n"    
    " -l NUM - analyse scenario with level<=NUM\n"    
    " -p  - scenario printing options\n"
    "    M - skip non ME|MP; print E if ME, P if MP\n"
    "    v - skip field names with="
    "    l - level\n"
    "    e - me\n"
    "    p - mp\n"
    "    P - mp profile; per each leaf\n"
    "    E - episode vector; prefix order\n"
    "    R - rule vector;how many rules to apply; negative -=TCLOST; prefix order\n"
    "    o - print compressed scenario\n"
    "    d - print scenario &t notation\n"
    "    D - \n"
    "\n"
    "Examples:\n"    
    "Diagram for all scenarios"
    "dlsgen2 \"(a,(b,c))\" \"(a,(b,c))\" -ga -da -pilaE -l4\n\n"
    "Skip variable info:\n"
    "dlsgen2 \"(a,(b,c))\" \"(a,(b,c))\" -ga -da -pilaEv -l4\n\n"
    "Diagram for PG-scenarios:\n"
    "dlsgen2 \"(a,(b,c))\" \"(a,(b,c))\" -gs -da -pilaEv -l4\n\n"
    );
  exit(0);
}


int main(int argc, char **argv)
{ 
  int verbose=0;
  if (argc < 2) usage(argc, argv);


  // if (argc<3)
  //   {
  //     fprintf(stderr,
  //       "Dls tree generator\n"
  //       "Usage: %s gtree stree OPT \n" 
  //       "OPT: \n" 
  //       " i - print number of scenarios\n" 
  //       " s - fixed speciations=only TMOVE-1 rule\n"                

  //       " eE - print me score/profile\n"
  //       " pP - print mp score/profile\n"        
  //       " l - print level\n"
  //       " G - gse output\n"
  //       " R - print rule profile\n"
  //       " d - &t output\n"
  //       " D - as above but diagram based ouput (single line per each tree)\n"
  //       " M - comput minme, minmp\n"                
  //       " m - comput minmp\n"
  //       " g - rule-edges\n"       
  //       " x - show only minme/minmp in the graph\n"
  //       " c - compress?\n"
  //       " A - heur: absorption rule\n"
  //     ,argv[0]);      
  //     exit(-1);
  //   }

  int optgenscen=GENALL; // ALL
  int optcnt=0; 
  int optcountsonly=0;
  char *opts=strdup("");
  char opt;
  int optdiagram=0;
  while ((opt = getopt(argc, argv, "d:g:iv:l:p:"))!= -1)
    switch (opt) {
    case 'g':
      if (strchr(optarg,'a')) optgenscen=GENALL;
      if (strchr(optarg,'s')) optgenscen=GENFIXEDSPEC;
      if (strchr(optarg,'c')) optcnt=1;
      if (strchr(optarg,'C')) optcnt=optcountsonly=1;
      break;

    case 'v':
      if (sscanf(optarg, "%d", &verbose) != 1) {
        fprintf(stderr,"Number expected in -v\n");
        exit(-1);
      }
      break;
    case 'd': //old d          
      if (strchr(optarg,'e')) optdiagram|=GENME;
      if (strchr(optarg,'p')) optdiagram|=GENMP;
      if (strchr(optarg,'a')) optdiagram=GENALL;
      break;
    case 'l': 
      if (sscanf(optarg, "%d", &levellimit) != 1) {
        fprintf(stderr,"Number expected in -l\n");
        exit(-1);
      }
      break;
    case 'p':      
      opts=strdup(optarg);      
      break;


    default:
      fprintf(stderr,"Unknown option: %c\n",((char) opt));
      exit(-1);
  } //switch

  if (argc-optind!=2)
  {
    printf("Expected two arguments: gene tree and species tree\n");
    exit(-1);
  }

  Node *gt = parseTree(argv[optind]);
  Node *st = parseTree(argv[optind+1]);

  if (verbose>0) { 
    printf("Gene and species trees: ");
    printf("\"");
    ppTree(gt,0);
    printf("\" \""); 
    ppTree(st,0);
    printf("\"\n");
  }
  
  LCA *s=genLCA(st);
  map(s,gt);
  
  
  //bool fixedspec=strchr(argv[3],'s');
  long scencnt=scenCount(gt,st,optgenscen==GENFIXEDSPEC);
  int gtsize=treeleafnum(gt);
  int stsize=treeleafnum(st);

  if (verbose>0) { 
    if (optgenscen==GENALL) printf("Model: all scenarios\n");    
    if (optgenscen==GENFIXEDSPEC) printf("Model: PG model\n");    
  }

  if (optcnt) printf("Scenarios count: %ld\n",scencnt);    
  
  Node **gtnodes=(Node**)malloc((2*gtsize+1)*sizeof(Node*));
  int lfp=0,itp=gtsize;
  insertNodes(gt,gtnodes,&lfp,&itp,0); // set id nodes

  long i=0;
  int mes,mps;  


  // Gen counts only
  if (optcountsonly) exit(0);


  // Big tab for all scenarios
  Scen **scentab=(Scen**)malloc(scencnt*sizeof(Scen*));
  if (!scentab)
  {
    perror("Not enough memory\n");
    exit(-1);
  }  
  // printf("%ld\n",scencnt*sizeof(Scen));
  int minmps=gtsize*4;
  int minmes=gtsize*4;  
  int maxlevel=0;
  int minmpcnt=0;
  int minmecnt=0;

  for (i=0;i<scencnt;i++)
  {
    if (verbose>2)
      if (!(i%10000)) printf("%ld/%ld\n",i,scencnt);    
    scentab[i]=genscen(NULL,i,gt,st,s,optgenscen==GENFIXEDSPEC,stsize,gtsize);     
    if (scentab[i]->level>maxlevel) maxlevel=scentab[i]->level;    
    if (scentab[i]->me<minmes) { minmes=scentab[i]->me; minmecnt=0; }
    if (scentab[i]->mp<minmps) { minmps=scentab[i]->mp; minmpcnt=0; }                   
    if (scentab[i]->me==minmes) minmecnt++;
    if (scentab[i]->mp==minmps) minmpcnt++;               
  }

  if (verbose)
  {
    printf("Levels: %d\n",maxlevel+1);
    printf("ME: %d\n",minmes);
    printf("MP: %d\n",minmps);
    printf("Min ME scenarios cnt: %d\n",minmecnt);
    printf("Min MP scenarios cnt: %d\n",minmpcnt);
  }

  int *scentabsorted=NULL;
  
  int levelcount[maxlevel+2];
  
  if (optdiagram)
  {    
    // sort by levels
    maxlevel+=1;    
    int i;
    for (i=0;i<=maxlevel;i++) levelcount[i]=0;
    for (i=0;i<scencnt;i++) levelcount[scentab[i]->level+1]++;
    for (i=1;i<=maxlevel;i++) levelcount[i]+=levelcount[i-1];
    scentabsorted=(int*)malloc(scencnt*sizeof(int));
    for (i=0;i<scencnt;i++) 
      scentabsorted[levelcount[scentab[i]->level]++]=i;

    for (i=0;i<=maxlevel;i++) levelcount[i]=0;
    for (i=0;i<scencnt;i++) levelcount[scentab[i]->level+1]++;    
    
    for (i=1;i<=maxlevel;i++) levelcount[i]+=levelcount[i-1];
    maxlevel--;
  }


  int mmlevelcount[maxlevel+2];
  for (i=0;i<=maxlevel;i++) mmlevelcount[i]=0;
    
  int minmingraph=(optdiagram&(GENMP|GENME));
        
  for (i=0;i<scencnt;i++)   
    {  
      long int id=i;
      if (scentabsorted) id=scentabsorted[id]; // sorted diagram
      if (scentab[id]->level>levellimit) break;
     
      if (minmingraph)
      {       
        if (( minmingraph&GENMP && minmps==scentab[id]->mp) || 
            ( minmingraph&GENME && minmes==scentab[id]->me))              
        {
          genscen(NULL,id,gt,st,s,optgenscen==GENFIXEDSPEC,stsize,gtsize,1);     
          scentab[id]->printinfo(minmes,minmps,stsize, gtsize, gt, opts);     
          mmlevelcount[scentab[id]->level]++;    
          if (opts_d) printf("&t "); ppdls(gt,0,0);               
          printf("\n");      
        }
       
        continue; // skip rest
      }
      else
      {

        if (opts_o || opts_d || opts_D)
          genscen(NULL,id,gt,st,s,optgenscen==GENFIXEDSPEC,stsize,gtsize,1);        

        if (!scentab[id]->printinfo(minmes,minmps,stsize, gtsize, gt, opts)) continue;
              
        if (optdiagram)
        {
          int clevel=scentab[id]->level;
          char buf[100000];          
          char *p=buf;
                              
          int j,k;            
          for (k=0;k<2;k++)
          {
              if (!strchr(opts,'v')) 
              if (k) printf("Rules=");
              else printf("Conn="); 
              printf("[");
              int first=1;
              if (clevel<maxlevel && clevel<=levellimit)
              for (j=levelcount[clevel+1];j<levelcount[clevel+2];j++)
              {
                long int id2=scentabsorted[j];              
                int rtr=ruletransf(scentab[id],scentab[id2],gtsize);
                if (rtr)
                {
                  if (!first) printf(","); else first=0;
                  if (k) 
                  {                    
                    if (rtr<0)  //printf("-1,"); // spec
                    printf("0");
                    else printf("1"); 
                  }
                  else
                  {                  
                    printf("%ld",id2);
                  }
                }
              }              
              printf("] ");
          } //for k            
          
          
        }        
      }

      if (opts_d) { printf("&t "); ppdls(gt,0,0);  }             
      printf("\n");
    }

    if (minmingraph)
    {

      int mi=0;
      for (i=0;i<=maxlevel;i++) 
        if (mmlevelcount[i]>0) mi=i;
      printf("MAXLEVEL=%d\n",mmlevelcount[mi]);
    }
        
    if (opts_m)
      printf("minme=%d minmp=%d\n",minmes,minmps);

  

  // if (strchr(opts,'A'))
  // {

  //   int minmetab[minmecnt+1],ime=0;
  //   int minmptab[minmpcnt+1],imp=0;
  //   for (i=0;i<scencnt;i++) 
  //   {
  //     if (scentab[i]->me==minmes) minmetab[ime++]=i;      
  //     if (scentab[i]->mp==minmps) minmptab[imp++]=i;
  //   }    
  //   minmetab[ime]=-1;
  //   minmptab[imp]=-1;
  //   printf("minme=%d minmp=%d\n",minmes,minmps);
  //   printf("minmecnt=%d minmpcnt=%d\n",minmecnt,minmpcnt);
  //   Scen* lca=genscen(NULL,0,gt,st,s,optgenscen==GENFIXEDSPEC,stsize,gtsize);      // LCA
  //   ruleabsorption2(lca,gt,st,gtnodes,s,stsize,gtsize,scentab,minmetab,minmptab,opts,minmes,minmps);

    
  // }

  free(scentab);
  return 0;
}

