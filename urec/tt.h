
/************************************************************************
   Unrooted REConciliation version 1.00 
   (c) Copyright 2005-2006 by Pawel Gorecki
   Written by P.Gorecki.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. 
*************************************************************************/

#ifndef _SPR__
#define _SPR__

#include <iostream>
#include <map>
using namespace std;


void spr(UTree *gt, RTree *sp, UNode *t, UNode3 *e, DLCost cost);
/*
 *  SPR on gt. Result gt'
 *  t - the root of the pruned subtree
 *  <e,e->pn> - the edge where the subtree will be regrafted
 *  ming - ?
 * 
 *  Assumptions:
 *  gt - has correct mappings
 *  cost - is the minimal rec. cost
 *  
 */
