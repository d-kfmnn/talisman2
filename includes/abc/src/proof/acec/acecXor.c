/**CFile****************************************************************

  FileName    [acecXor.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [CEC for arithmetic circuits.]

  Synopsis    [Detection of XOR trees.]

  Author      [Alan Mishchenko]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - June 20, 2005.]

  Revision    [$Id: acecXor.c,v 1.00 2005/06/20 00:00:00 alanmi Exp $]

***********************************************************************/

#include "acecInt.h"
#include "misc/vec/vecWec.h"
#include "misc/extra/extra.h"
#include "misc/util/utilTruth.h"
#include "misc/vec/vecInt.h"

#ifdef ABC_USE_PTHREADS
#include <pthread.h>
#endif

ABC_NAMESPACE_IMPL_START

////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

// Global verbose level for controlling debug output
// 0 = silent, 1 = basic info, 2+ = detailed debug
static int g_acecVerbose = 0;

// Set verbose level for ACEC module
void Acec_SetVerbose(int level) { g_acecVerbose = level; }
int  Acec_GetVerbose(void)      { return g_acecVerbose; }

////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Acec_CheckXors( Gia_Man_t * p, Vec_Int_t * vXors )
{
    Gia_Obj_t * pFan0, * pFan1;
    Vec_Int_t * vCount2 = Vec_IntAlloc( Gia_ManObjNum(p) );
    int i, Entry, Count = 0;
    for ( i = 0; 4*i < Vec_IntSize(vXors); i++ )
        if ( Vec_IntEntry(vXors, 4*i+3) == 0 )
            Vec_IntAddToEntry( vCount2, Vec_IntEntry(vXors, 4*i), 1 );
    Vec_IntForEachEntry( vCount2, Entry, i )
        if ( Entry > 1 )
            printf( "*** Obj %d has %d two-input XOR cuts.\n", i, Entry ), Count++;
        else if ( Entry == 1 && Gia_ObjRecognizeExor(Gia_ManObj(p, i), &pFan0, &pFan1) )
            printf( "*** Obj %d cannot be recognized as XOR.\n", i );
    if ( Count == 0 )
        printf( "*** There no multiple two-input XOR cuts.\n" );
    Vec_IntFree( vCount2 );
}

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Vec_Int_t * Acec_OrderTreeRoots( Gia_Man_t * p, Vec_Int_t * vAdds, Vec_Int_t * vXorRoots, Vec_Int_t * vRanks )
{
    Vec_Int_t * vOrder = Vec_IntAlloc( Vec_IntSize(vXorRoots) );
    Vec_Int_t * vMove = Vec_IntStartFull( Vec_IntSize(vXorRoots) );
    int i, k, Entry, This;
    // iterate through adders and for each try mark the next one
    for ( i = 0; 6*i < Vec_IntSize(vAdds); i++ )
    {
        int Node = Vec_IntEntry(vAdds, 6*i+4);
        if ( Vec_IntEntry(vRanks, Node) == -1 )
            continue;
        for ( k = 0; k < 3; k++ )
        {
            int Fanin = Vec_IntEntry(vAdds, 6*i+k);
            if ( Vec_IntEntry(vRanks, Fanin) == -1 )
                continue;
            //printf( "%4d:  %2d -> %2d\n", Node, Vec_IntEntry(vRanks, Node), Vec_IntEntry(vRanks, Fanin) );
            Vec_IntWriteEntry( vMove, Vec_IntEntry(vRanks, Node), Vec_IntEntry(vRanks, Fanin) );
        }
    }
//Vec_IntPrint( vMove );
    // find reodering
    Vec_IntForEachEntry( vMove, Entry, i )
        if ( Entry == -1 && Vec_IntFind(vMove, i) >= 0 )
            break;
    assert( i < Vec_IntSize(vMove) );
    while ( 1 )
    {
        Vec_IntPush( vOrder, Vec_IntEntry(vXorRoots, i) );
        Entry = i;
        Vec_IntForEachEntry( vMove, This, i )
            if ( This == Entry )
                break;
        if ( i == Vec_IntSize(vMove) )
            break;
    }
    Vec_IntFree( vMove );
//Vec_IntPrint( vOrder );
    return vOrder;
}
                       
/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
// marks XOR outputs
Vec_Bit_t * Acec_MapXorOuts( Gia_Man_t * p, Vec_Int_t * vXors )
{
    Vec_Bit_t * vMap = Vec_BitStart( Gia_ManObjNum(p) ); int i;
    for ( i = 0; 4*i < Vec_IntSize(vXors); i++ )
        Vec_BitWriteEntry( vMap, Vec_IntEntry(vXors, 4*i), 1 );
    return vMap;
}
// marks XOR outputs participating in trees
Vec_Bit_t * Acec_MapXorOuts2( Gia_Man_t * p, Vec_Int_t * vXors, Vec_Int_t * vRanks )
{
    Vec_Bit_t * vMap = Vec_BitStart( Gia_ManObjNum(p) ); int i;
    for ( i = 0; 4*i < Vec_IntSize(vXors); i++ )
        if ( Vec_IntEntry(vRanks, Vec_IntEntry(vXors, 4*i)) != -1 )
            Vec_BitWriteEntry( vMap, Vec_IntEntry(vXors, 4*i), 1 );
    return vMap;
}
// marks MAJ outputs
Vec_Bit_t * Acec_MapMajOuts( Gia_Man_t * p, Vec_Int_t * vAdds )
{
    Vec_Bit_t * vMap = Vec_BitStart( Gia_ManObjNum(p) ); int i;
    for ( i = 0; 6*i < Vec_IntSize(vAdds); i++ )
        Vec_BitWriteEntry( vMap, Vec_IntEntry(vAdds, 6*i+4), 1 );
    return vMap;
}
// marks MAJ outputs participating in trees
Vec_Int_t * Acec_MapMajOuts2( Gia_Man_t * p, Vec_Int_t * vAdds, Vec_Int_t * vRanks )
{
    Vec_Int_t * vMap = Vec_IntStartFull( Gia_ManObjNum(p) ); int i;
    for ( i = 0; 6*i < Vec_IntSize(vAdds); i++ )
        if ( Vec_IntEntry(vRanks, Vec_IntEntry(vAdds, 6*i+4)) != -1 )
            Vec_IntWriteEntry( vMap, Vec_IntEntry(vAdds, 6*i+4), i );
    return vMap;
}
// marks nodes appearing as fanins to XORs
Vec_Bit_t * Acec_MapXorIns( Gia_Man_t * p, Vec_Int_t * vXors )
{
    Vec_Bit_t * vMap = Vec_BitStart( Gia_ManObjNum(p) ); int i;
    for ( i = 0; 4*i < Vec_IntSize(vXors); i++ )
    {
        Vec_BitWriteEntry( vMap, Vec_IntEntry(vXors, 4*i+1), 1 );
        Vec_BitWriteEntry( vMap, Vec_IntEntry(vXors, 4*i+2), 1 );
        Vec_BitWriteEntry( vMap, Vec_IntEntry(vXors, 4*i+3), 1 );
    }
    return vMap;
}
// collects XOR roots (XOR nodes not appearing as fanins of other XORs)
Vec_Int_t * Acec_FindXorRoots( Gia_Man_t * p, Vec_Int_t * vXors )
{
    Vec_Bit_t * vMapXorIns = Acec_MapXorIns( p, vXors );
    Vec_Int_t * vXorRoots = Vec_IntAlloc( 100 );  int i;
    for ( i = 0; 4*i < Vec_IntSize(vXors); i++ )
        if ( !Vec_BitEntry(vMapXorIns, Vec_IntEntry(vXors, 4*i)) )
            Vec_IntPushUniqueOrder( vXorRoots, Vec_IntEntry(vXors, 4*i) );
    Vec_BitFree( vMapXorIns );
    return vXorRoots;
}
// collects XOR trees belonging to each of XOR roots
Vec_Int_t * Acec_RankTrees( Gia_Man_t * p, Vec_Int_t * vXors, Vec_Int_t * vXorRoots )
{
    Vec_Int_t * vDoubles = Vec_IntAlloc( 100 );
    int i, k, Entry;
    // map roots into their ranks
    Vec_Int_t * vRanks = Vec_IntStartFull( Gia_ManObjNum(p) ); 
    Vec_IntForEachEntry( vXorRoots, Entry, i )
        Vec_IntWriteEntry( vRanks, Entry, i );
    // map nodes into their ranks
    for ( i = Vec_IntSize(vXors)/4 - 1; i >= 0; i-- )
    {
        int Root = Vec_IntEntry( vXors, 4*i );
        int Rank = Vec_IntEntry( vRanks, Root );
        // skip XORs that are not part of any tree
        if ( Rank == -1 )
            continue;
        // iterate through XOR inputs
        for ( k = 1; k < 4; k++ )
        {
            int Node = Vec_IntEntry( vXors, 4*i+k );
            if ( Node == 0 ) // HA
                continue;
            Entry = Vec_IntEntry( vRanks, Node );
            if ( Entry == Rank ) // the same tree
                continue;
            if ( Entry == -1 )
                Vec_IntWriteEntry( vRanks, Node, Rank );
            else
                Vec_IntPush( vDoubles, Node );

            if ( Entry != -1 && Gia_ObjIsAnd(Gia_ManObj(p, Node)))
            printf( "Xor node %d belongs to Tree %d and Tree %d.\n", Node, Entry, Rank );
        }
    }
    // remove duplicated entries
    Vec_IntForEachEntry( vDoubles, Entry, i )
        Vec_IntWriteEntry( vRanks, Entry, -1 );
    Vec_IntFree( vDoubles );
    return vRanks;
}
// collects leaves of each XOR tree
Vec_Wec_t * Acec_FindXorLeaves( Gia_Man_t * p, Vec_Int_t * vXors, Vec_Int_t * vAdds, Vec_Int_t * vXorRoots, Vec_Int_t * vRanks, Vec_Wec_t ** pvAddBoxes )
{
    Vec_Bit_t * vMapXors = Acec_MapXorOuts2( p, vXors, vRanks );
    Vec_Int_t * vMapMajs = Acec_MapMajOuts2( p, vAdds, vRanks );
    Vec_Wec_t * vXorLeaves = Vec_WecStart( Vec_IntSize(vXorRoots) );
    Vec_Wec_t * vAddBoxes  = Vec_WecStart( Vec_IntSize(vXorRoots) );
    int i, k;
    for ( i = 0; 4*i < Vec_IntSize(vXors); i++ )
    {
        int Xor = Vec_IntEntry(vXors, 4*i);
        int Rank = Vec_IntEntry(vRanks, Xor);
        if ( Rank == -1 )
            continue;
        for ( k = 1; k < 4; k++ )
        {
            int Fanin = Vec_IntEntry(vXors, 4*i+k);
            //int RankFanin = Vec_IntEntry(vRanks, Fanin);
            if ( Fanin == 0 )
                continue;
            if ( Vec_BitEntry(vMapXors, Fanin) )
            {
                assert( Rank == Vec_IntEntry(vRanks, Fanin) );
                continue;
            }
//            if ( Vec_BitEntry(vMapXors, Fanin) && Rank == RankFanin )
//                continue;
            if ( Vec_IntEntry(vMapMajs, Fanin) == -1 ) // no adder driving this input
                Vec_WecPush( vXorLeaves, Rank, Fanin );
            else if ( Vec_IntEntry(vRanks, Xor) > 0 ) // save adder box
                Vec_WecPush( vAddBoxes, Rank-1, Vec_IntEntry(vMapMajs, Fanin) );
        }
    }
    Vec_BitFree( vMapXors );
    Vec_IntFree( vMapMajs );
    if ( pvAddBoxes )
        *pvAddBoxes = vAddBoxes;
    return vXorLeaves;
}
void Acec_CheckBoothPPs( Gia_Man_t * p, Vec_Wec_t * vLitLeaves )
{
    Vec_Bit_t * vMarked = Acec_MultMarkPPs( p );
    Vec_Int_t * vLevel;
    int i, k, iLit;
    Vec_WecForEachLevel( vLitLeaves, vLevel, i )
    {
        int CountPI = 0, CountB = 0, CountNB = 0;
        Vec_IntForEachEntry( vLevel, iLit, k )
            if ( !Gia_ObjIsAnd(Gia_ManObj(p, Abc_Lit2Var(iLit))) )
                CountPI++;
            else if ( Vec_BitEntry( vMarked, Abc_Lit2Var(iLit) ) )
                CountB++;
            else
                CountNB++;

        printf( "Rank %2d : Lits = %5d    PI = %d  Booth = %5d  Non-Booth = %5d\n", i, Vec_IntSize(vLevel), CountPI, CountB, CountNB );
    }
    Vec_BitFree( vMarked );
}
Acec_Box_t * Acec_FindBox( Gia_Man_t * p, Vec_Int_t * vAdds, Vec_Wec_t * vAddBoxes, Vec_Wec_t * vXorLeaves, Vec_Int_t * vXorRoots )
{
    extern Vec_Int_t * Acec_TreeCarryMap( Gia_Man_t * p, Vec_Int_t * vAdds, Vec_Wec_t * vBoxes );
    extern void Acec_TreePhases_rec( Gia_Man_t * p, Vec_Int_t * vAdds, Vec_Int_t * vMap, int Node, int fPhase, Vec_Bit_t * vVisit );
    extern void Acec_TreeVerifyPhases( Gia_Man_t * p, Vec_Int_t * vAdds, Vec_Wec_t * vBoxes );
    extern void Acec_TreeVerifyPhases2( Gia_Man_t * p, Vec_Int_t * vAdds, Vec_Wec_t * vBoxes );

    int MaxRank = Vec_WecSize( vAddBoxes );
    Vec_Bit_t * vVisit  = Vec_BitStart( Vec_IntSize(vAdds)/6 );
    Vec_Bit_t * vIsLeaf = Vec_BitStart( Gia_ManObjNum(p) );
    Vec_Bit_t * vIsRoot = Vec_BitStart( Gia_ManObjNum(p) );
    Vec_Int_t * vLevel, * vLevel2, * vMap;
    int i, j, k, Box, Node;

    Acec_Box_t * pBox = ABC_CALLOC( Acec_Box_t, 1 );
    pBox->pGia        = p;
    pBox->vAdds       = vAddBoxes; // Vec_WecDup( vAddBoxes );
    pBox->vLeafLits   = Vec_WecStart( MaxRank + 0 );
    pBox->vRootLits   = Vec_WecStart( MaxRank + 0 );

    assert( Vec_WecSize(vAddBoxes) == Vec_WecSize(vXorLeaves) );
    assert( Vec_WecSize(vAddBoxes) == Vec_IntSize(vXorRoots) );

    // collect boxes; mark inputs/outputs
    Vec_WecForEachLevel( pBox->vAdds, vLevel, i )
    Vec_IntForEachEntry( vLevel, Box, k )
    {
        Vec_BitWriteEntry( vIsLeaf, Vec_IntEntry(vAdds, 6*Box+0), 1 );
        Vec_BitWriteEntry( vIsLeaf, Vec_IntEntry(vAdds, 6*Box+1), 1 );
        Vec_BitWriteEntry( vIsLeaf, Vec_IntEntry(vAdds, 6*Box+2), 1 );
        Vec_BitWriteEntry( vIsRoot, Vec_IntEntry(vAdds, 6*Box+3), 1 );
        Vec_BitWriteEntry( vIsRoot, Vec_IntEntry(vAdds, 6*Box+4), 1 );
    }
    // sort each level
    Vec_WecForEachLevel( pBox->vAdds, vLevel, i )
        Vec_IntSort( vLevel, 0 );

    // set phases starting from roots
    vMap = Acec_TreeCarryMap( p, vAdds, pBox->vAdds );
    Vec_WecForEachLevelReverse( pBox->vAdds, vLevel, i )
        Vec_IntForEachEntry( vLevel, Box, k )
            if ( !Vec_BitEntry( vIsLeaf, Vec_IntEntry(vAdds, 6*Box+4) ) )
            {
                //printf( "Pushing phase of output %d of box %d\n", Vec_IntEntry(vAdds, 6*Box+4), Box );
                Acec_TreePhases_rec( p, vAdds, vMap, Vec_IntEntry(vAdds, 6*Box+4), Vec_IntEntry(vAdds, 6*Box+2) != 0, vVisit );
            }
    Acec_TreeVerifyPhases( p, vAdds, pBox->vAdds );
    Acec_TreeVerifyPhases2( p, vAdds, pBox->vAdds );
    Vec_BitFree( vVisit );
    Vec_IntFree( vMap );

    // collect inputs/outputs
    Vec_BitWriteEntry( vIsRoot, 0, 1 );
    Vec_WecForEachLevel( pBox->vAdds, vLevel, i )
        Vec_IntForEachEntry( vLevel, Box, j )
        {
            for ( k = 0; k < 3; k++ )
                if ( !Vec_BitEntry( vIsRoot, Vec_IntEntry(vAdds, 6*Box+k) ) )
                    Vec_WecPush( pBox->vLeafLits, i, Abc_Var2Lit(Vec_IntEntry(vAdds, 6*Box+k), Acec_SignBit2(vAdds, Box, k)) );
            for ( k = 3; k < 5; k++ )
                if ( !Vec_BitEntry( vIsLeaf, Vec_IntEntry(vAdds, 6*Box+k) ) )
                    Vec_WecPush( pBox->vRootLits, k == 4 ? i + 1 : i, Abc_Var2Lit(Vec_IntEntry(vAdds, 6*Box+k), Acec_SignBit2(vAdds, Box, k)) );
            if ( Vec_IntEntry(vAdds, 6*Box+2) == 0 && Acec_SignBit2(vAdds, Box, 2) )
                Vec_WecPush( pBox->vLeafLits, i, 1 );
        }
    Vec_BitFree( vIsLeaf );
    Vec_BitFree( vIsRoot );

    // collect last bit
    vLevel  = Vec_WecEntry( pBox->vLeafLits, Vec_WecSize(pBox->vLeafLits)-1 );
    vLevel2 = Vec_WecEntry( vXorLeaves, Vec_WecSize(vXorLeaves)-1 );
    if ( Vec_IntSize(vLevel) == 0 && Vec_IntSize(vLevel2) > 0 )
    {
        Vec_IntForEachEntry( vLevel2, Node, k )
            Vec_IntPush( vLevel, Abc_Var2Lit(Node, 0) );
    }
    vLevel  = Vec_WecEntry( pBox->vRootLits, Vec_WecSize(pBox->vRootLits)-1 );
    Vec_IntFill( vLevel, 1, Abc_Var2Lit(Vec_IntEntryLast(vXorRoots), 0) );

    // sort each level
    Vec_WecForEachLevel( pBox->vLeafLits, vLevel, i )
        Vec_IntSort( vLevel, 0 );
    Vec_WecForEachLevel( pBox->vRootLits, vLevel, i )
        Vec_IntSort( vLevel, 1 );

    //Acec_CheckBoothPPs( p, pBox->vLeafLits );
    return pBox;
}

Acec_Box_t * Acec_ProduceBox( Gia_Man_t * p, int fVerbose )
{
    extern void Acec_TreeVerifyConnections( Gia_Man_t * p, Vec_Int_t * vAdds, Vec_Wec_t * vBoxes );

    abctime clk = Abc_Clock();
    Acec_Box_t * pBox = NULL;
    Vec_Int_t * vXors, * vAdds = Ree_ManComputeCuts( p, &vXors, 0 );
    Vec_Int_t * vTemp, * vXorRoots = Acec_FindXorRoots( p, vXors ); 
    Vec_Int_t * vRanks = Acec_RankTrees( p, vXors, vXorRoots ); 
    Vec_Wec_t * vXorLeaves, * vAddBoxes = NULL; 

    Gia_ManLevelNum(p);

    //Acec_CheckXors( p, vXors );

    //Ree_ManPrintAdders( vAdds, 1 );
    if ( fVerbose )
        printf( "Detected %d full-adders and %d half-adders.  Found %d XOR-cuts.  ", Ree_ManCountFadds(vAdds), Vec_IntSize(vAdds)/6-Ree_ManCountFadds(vAdds), Vec_IntSize(vXors)/4 );
    if ( fVerbose )
        Abc_PrintTime( 1, "Time", Abc_Clock() - clk );

    vXorRoots = Acec_OrderTreeRoots( p, vAdds, vTemp = vXorRoots, vRanks );
    Vec_IntFree( vTemp );
    Vec_IntFree( vRanks );

    vRanks = Acec_RankTrees( p, vXors, vXorRoots ); 
    vXorLeaves = Acec_FindXorLeaves( p, vXors, vAdds, vXorRoots, vRanks, &vAddBoxes ); 
    Vec_IntFree( vRanks );

    //printf( "XOR roots after reordering: \n" );
    //Vec_IntPrint( vXorRoots );
    //printf( "XOR leaves: \n" );
    //Vec_WecPrint( vXorLeaves, 0 );
    //printf( "Adder boxes: \n" );
    //Vec_WecPrint( vAddBoxes, 0 );

    Acec_TreeVerifyConnections( p, vAdds, vAddBoxes );

    pBox = Acec_FindBox( p, vAdds, vAddBoxes, vXorLeaves, vXorRoots );
    //Vec_WecFree( vAddBoxes );

    if ( fVerbose )
        Acec_TreePrintBox( pBox, vAdds );

    Vec_IntFree( vXorRoots );
    Vec_WecFree( vXorLeaves );

    Vec_IntFree( vXors );
    Vec_IntFree( vAdds );

    return pBox;
}


/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Gia_ManTestXor( Gia_Man_t * p )
{
    Vec_Wrd_t * vSimsPi = Vec_WrdStartTruthTables( Gia_ManCiNum(p) );
    Vec_Wrd_t * vSims   = Gia_ManSimPatSimOut( p, vSimsPi, 1 );
    int n, i, nWords = Vec_WrdSize(vSimsPi) / Gia_ManCiNum(p);
    Gia_Obj_t * pObj; Vec_Wrd_t * vSims2;
    Gia_ManForEachAnd( p, pObj, i )
    {
        Gia_Obj_t Obj = *pObj;
        for ( n = 0; n < 2; n++ )
        {
            if ( n ) 
            {
                pObj->iDiff1 = pObj->iDiff0;
                pObj->fCompl1 = pObj->fCompl0;
            }
            else
            {
                pObj->iDiff0 = pObj->iDiff1;
                pObj->fCompl0 = pObj->fCompl1;
            }
            vSims2 = Gia_ManSimPatSimOut( p, vSimsPi, 1 );
            printf( "%2d %2d : %5d\n", i, n, Abc_TtCountOnesVecXor(Vec_WrdArray(vSims), Vec_WrdArray(vSims2), Vec_WrdSize(vSims2)) );
            Vec_WrdFree( vSims2 );
            *pObj = Obj;
        }
    }
    Vec_WrdFree( vSimsPi );
    Vec_WrdFree( vSims );
    nWords = 0;
}
// dfs collect xor internal nodes 

void recursiveXor(Gia_Man_t * p, Vec_Int_t * vInput, int xor2_traversal, Vec_Int_t * xor2_record) {
    Vec_Int_t * pStack = Vec_IntAlloc(1000);
    Vec_Int_t * pVisited = Vec_IntAlloc(1000);

    if (p == NULL || vInput == NULL || xor2_record == NULL) {
        printf("Error: NULL pointer passed to recursiveXor\n");
        Vec_IntFree(pStack);
        Vec_IntFree(pVisited);
        return;
    }

    if (xor2_traversal < 0 || xor2_traversal >= Gia_ManObjNum(p)) {
        printf("Warning: Invalid initial object ID %d in recursiveXor (valid range: 0-%d)\n", 
               xor2_traversal, Gia_ManObjNum(p) - 1);
        Vec_IntFree(pStack);
        Vec_IntFree(pVisited);
        return;
    }

    Vec_IntPush(pStack, xor2_traversal);

    while (Vec_IntSize(pStack) > 0) {
        int current_id = Vec_IntPop(pStack);

        if (Vec_IntFind(pVisited, current_id) != -1) {
            continue;
        }
        Vec_IntPush(pVisited, current_id);

        if (current_id < 0 || current_id >= Gia_ManObjNum(p)) {
            printf("Warning: Invalid object ID %d encountered in recursiveXor\n", current_id);
            continue;
        }

        Gia_Obj_t * pObj = Gia_ManObj(p, current_id);
        if (pObj == NULL) {
            printf("Warning: NULL object returned for ID %d in recursiveXor\n", current_id);
            continue;
        }

        if (!Gia_ObjIsAnd(pObj)) {
            Vec_IntPushUnique(xor2_record, current_id);
            continue;
        }

        int fanin0 = Gia_ObjFaninId0p(p, pObj);
        int fanin1 = Gia_ObjFaninId1p(p, pObj);

        if (fanin0 < 0 || fanin0 >= Gia_ManObjNum(p)) {
            printf("Warning: Invalid fanin0 ID %d for object %d\n", fanin0, current_id);
            continue;
        }
        if (fanin1 < 0 || fanin1 >= Gia_ManObjNum(p)) {
            printf("Warning: Invalid fanin1 ID %d for object %d\n", fanin1, current_id);
            continue;
        }

        Vec_IntPushUnique(xor2_record, fanin0);
        Vec_IntPushUnique(xor2_record, fanin1);

        if (Vec_IntFind(vInput, fanin0) == -1 && Vec_IntFind(pVisited, fanin0) == -1) {
            Vec_IntPush(pStack, fanin0);
        }
        if (Vec_IntFind(vInput, fanin1) == -1 && Vec_IntFind(pVisited, fanin1) == -1) {
            Vec_IntPush(pStack, fanin1);
        }
    }

    Vec_IntFree(pStack);
    Vec_IntFree(pVisited);
}



static int is_subset_of(Vec_Int_t *group1, Vec_Int_t *group2) {
    Vec_IntSort(group2, 0);
    for (int i = 0; i < Vec_IntSize(group1); i++) {
        int value = Vec_IntEntry(group1, i);
        if (Vec_IntFind(group2, value) == -1)
            return 0;
    }
    return 1;
}

int is_group_subset_of_list(Vec_Int_t *group, Vec_Int_t *list) {
    Vec_Int_t *current_group = Vec_IntAlloc(10);
    int is_subset = 0;
   
    for (int i = 0; i < Vec_IntSize(list); i++) {
        int value = Vec_IntEntry(list, i);
        if (value == 0) {
            if (is_subset_of(group, current_group)) {
                is_subset = 1;
                break;
            }
            Vec_IntClear(current_group);
        } else {
            Vec_IntPush(current_group, value);
        }
    }

    Vec_IntFree(current_group);
    return is_subset;
}

static void extract_group_heads(Vec_Int_t *grouped, Vec_Int_t *outHeads)
{
    if (!grouped || !outHeads) return;
    int n = Vec_IntSize(grouped);
    int expectHead = 1;
    for (int i = 0; i < n; ++i) {
        int v = Vec_IntEntry(grouped, i);
        if (v == 0) {
            expectHead = 1;
        } else if (expectHead) {
            Vec_IntPush(outHeads, v);
            expectHead = 0;
        }
    }
}

void Vec_IntSortUnique(Vec_Int_t *v)
{
    if (!v) return;
    if (Vec_IntSize(v) <= 1) return;
    Vec_IntSort(v, 0);
    int w = 0, n = Vec_IntSize(v);
    int prev = 0x7fffffff;
    for (int i = 0; i < n; ++i) {
        int x = Vec_IntEntry(v, i);
        if (i == 0 || x != prev) {
            Vec_IntWriteEntry(v, w++, x);
            prev = x;
        }
    }
    Vec_IntShrink(v, w);
}

void UniqueIdsInPlace(Gia_Man_t *pMan, Vec_Int_t *v) {
    if (!v) return;
    int n = Gia_ManObjNum(pMan);
    unsigned char *mark = (unsigned char *)ABC_CALLOC(unsigned char, n);
    int i, x, w = 0;
    Vec_IntForEachEntry(v, x, i) {
        if (x <= 0 || x >= n) continue;  
        if (mark[x]) continue;           
        mark[x] = 1;
        Vec_IntWriteEntry(v, w++, x);    
    }
    Vec_IntShrink(v, w);
    ABC_FREE(mark);
}

// ======================== printing to json ========================

static void print_vec_grouped_json(FILE *f, const char *name, Vec_Int_t *v) {
    fprintf(f, "  \"%s\": [\n", name);
    int first_group = 1;
    Vec_Int_t *grp = Vec_IntAlloc(16);
    for (int i = 0; i < Vec_IntSize(v); ++i) {
        int val = Vec_IntEntry(v, i);
        if (val == 0) {
            if (Vec_IntSize(grp) > 0) {
                if (!first_group) fprintf(f, ",\n");
                fprintf(f, "    [");
                for (int k = 0; k < Vec_IntSize(grp); ++k) {
                    if (k) fprintf(f, ", ");
                    fprintf(f, "%d", Vec_IntEntry(grp, k));
                }
                fprintf(f, "]");
                first_group = 0;
                Vec_IntClear(grp);
            }
        } else {
            Vec_IntPush(grp, val);
        }
    }
    if (Vec_IntSize(grp) > 0) {
        if (!first_group) fprintf(f, ",\n");
        fprintf(f, "    [");
        for (int k = 0; k < Vec_IntSize(grp); ++k) {
            if (k) fprintf(f, ", ");
            fprintf(f, "%d", Vec_IntEntry(grp, k));
        }
        fprintf(f, "]");
        first_group = 0;
        Vec_IntClear(grp);
    }
    fprintf(f, first_group ? "  ]" : "\n  ]");
    Vec_IntFree(grp);
}

static int vec_group_count_by_zero( Vec_Int_t *v) {
    int groups = 0;
    int in_group = 0;
    int n = Vec_IntSize((Vec_Int_t *)v);
    for (int i = 0; i < n; ++i) {
        int val = Vec_IntEntry((Vec_Int_t *)v, i);
        if (val == 0) {
            if (in_group) {
                groups++;
                in_group = 0;
            }
        } else {
            in_group = 1;
        }
    }
    if (in_group) groups++;
    return groups;
}

// =================== extract ha/fa labeling ========================
// ======================================================================
// VecInt result bundle for HA/FA extraction
// ======================================================================
typedef unsigned char u8;

typedef struct {
    Vec_Int_t *xor2_list;
    Vec_Int_t *xor3_list;
    Vec_Int_t *maj2_list;
    Vec_Int_t *maj3_list;
    Vec_Int_t *xor2_maj2_combined;
    Vec_Int_t *xor3_maj3_combined;
    Vec_Int_t *xor_all_list;
    Vec_Int_t *xor_remaining;
} VecIntResults1;

// Forward declarations of helpers used but not shown in the original snippet
// You should keep your own implementations or link from other files.
extern void recursiveXor(Gia_Man_t *p, Vec_Int_t *vInput, int iRoot, Vec_Int_t *vRec);
extern int  is_group_subset_of_list(Vec_Int_t *group, Vec_Int_t *list);
extern int  vec_group_count_by_zero(Vec_Int_t *v);
extern void print_vec_grouped_json(FILE *f, const char *name, Vec_Int_t *v);
extern void extract_group_heads(Vec_Int_t *vCombined, Vec_Int_t *vHeads);
extern void UniqueIdsInPlace(Gia_Man_t *p, Vec_Int_t *v);

// ======================================================================
// HA/FA labeling extraction
// ======================================================================
static VecIntResults1 extract_ha_fa_labeling(Gia_Man_t *p, Vec_Int_t *vAdds, Vec_Int_t *vXors)
{
    VecIntResults1 res;
    int nObjs = Gia_ManObjNum(p);

    res.xor2_list           = Vec_IntAlloc(nObjs);
    res.xor3_list           = Vec_IntAlloc(nObjs);
    res.maj2_list           = Vec_IntAlloc(nObjs);
    res.maj3_list           = Vec_IntAlloc(nObjs);
    res.xor_all_list        = Vec_IntAlloc(nObjs);
    res.xor_remaining       = Vec_IntAlloc(nObjs);
    res.xor2_maj2_combined  = Vec_IntAlloc(nObjs);
    res.xor3_maj3_combined  = Vec_IntAlloc(nObjs);

    // 1) Walk through vAdds and build HA/FA groups and lists
    for (int i = 0; 6 * i < Vec_IntSize(vAdds); i++) {
        int in1 = Vec_IntEntry(vAdds, 6 * i + 0);
        int in2 = Vec_IntEntry(vAdds, 6 * i + 1);
        int in3 = Vec_IntEntry(vAdds, 6 * i + 2);
        int iXor = Vec_IntEntry(vAdds, 6 * i + 3);
        int iMaj = Vec_IntEntry(vAdds, 6 * i + 4);

        if (in3 == 0) {
            // Half-adder: XOR2 + AND2
            Vec_Int_t *xor2_record = Vec_IntAlloc(10);
            Vec_Int_t *maj2_record = Vec_IntAlloc(10);
            Vec_Int_t *vInput      = Vec_IntAlloc(2);

            Vec_IntPush(xor2_record, iXor);
            Vec_IntPushTwo(vInput, in1, in2);
            recursiveXor(p, vInput, iXor, xor2_record);
            Vec_IntAppend(res.xor2_list, xor2_record);
            Vec_IntPush(res.xor2_list, 0);

            Vec_IntPush(maj2_record, iMaj);
            recursiveXor(p, vInput, iMaj, maj2_record);
            Vec_IntAppend(res.maj2_list, maj2_record);
            Vec_IntPush(res.maj2_list, 0);

            // Unique union and canonical combined group:
            // [in1, in2, middle..., iXor, iMaj, 0]
            Vec_Int_t *unionSet = Vec_IntAlloc(Vec_IntSize(xor2_record) + Vec_IntSize(maj2_record));
            Vec_Int_t *combined = Vec_IntAlloc(32);

            Vec_IntAppend(unionSet, xor2_record);
            Vec_IntAppend(unionSet, maj2_record);
            Vec_IntUniqify(unionSet);

            Vec_IntPushTwo(combined, in1, in2);
            for (int k = 0; k < Vec_IntSize(unionSet); ++k) {
                int v = Vec_IntEntry(unionSet, k);
                if (v == in1 || v == in2 || v == iXor || v == iMaj)
                    continue;
                Vec_IntPush(combined, v);
            }
            Vec_IntPush(combined, iXor);
            Vec_IntPush(combined, iMaj);
            Vec_IntPush(combined, 0);
            Vec_IntAppend(res.xor2_maj2_combined, combined);

            Vec_IntFree(unionSet);
            Vec_IntFree(xor2_record);
            Vec_IntFree(maj2_record);
            Vec_IntFree(vInput);
            Vec_IntFree(combined);
        } else {
            // Full-adder: XOR3 + MAJ3
            Vec_Int_t *xor3_record = Vec_IntAlloc(10);
            Vec_Int_t *maj3_record = Vec_IntAlloc(10);
            Vec_Int_t *vInput      = Vec_IntAlloc(3);

            Vec_IntPush(xor3_record, iXor);
            Vec_IntPushThree(vInput, in1, in2, in3);
            recursiveXor(p, vInput, iXor, xor3_record);
            Vec_IntAppend(res.xor3_list, xor3_record);
            Vec_IntPush(res.xor3_list, 0);

            Vec_IntPush(maj3_record, iMaj);
            recursiveXor(p, vInput, iMaj, maj3_record);
            Vec_IntAppend(res.maj3_list, maj3_record);
            Vec_IntPush(res.maj3_list, 0);

            // Unique union and canonical combined group:
            // [in1, in2, in3, middle..., iXor, iMaj, 0]
            Vec_Int_t *unionSet = Vec_IntAlloc(Vec_IntSize(xor3_record) + Vec_IntSize(maj3_record));
            Vec_Int_t *combined = Vec_IntAlloc(32);

            Vec_IntAppend(unionSet, xor3_record);
            Vec_IntAppend(unionSet, maj3_record);
            Vec_IntUniqify(unionSet);

            Vec_IntPushThree(combined, in1, in2, in3);
            for (int k = 0; k < Vec_IntSize(unionSet); ++k) {
                int v = Vec_IntEntry(unionSet, k);
                if (v == in1 || v == in2 || v == in3 || v == iXor || v == iMaj)
                    continue;
                Vec_IntPush(combined, v);
            }
            Vec_IntPush(combined, iXor);
            Vec_IntPush(combined, iMaj);
            Vec_IntPush(combined, 0);
            Vec_IntAppend(res.xor3_maj3_combined, combined);

            Vec_IntFree(unionSet);
            Vec_IntFree(xor3_record);
            Vec_IntFree(maj3_record);
            Vec_IntFree(vInput);
            Vec_IntFree(combined);
        }
    }

    // 2) Build all XOR groups from vXors (mixed XOR2/XOR3, 0-separated)

    // Actually construct xor_all_list from 0-separated vXors:
    // - If group has [iXor,in1,in2] -> XOR2, reuse HA XOR2 logic.
    // - If group has [iXor,in1,in2,in3,...] -> XOR3, reuse FA XOR3 logic.
    Vec_Int_t *cur = Vec_IntAlloc(16);
    for (int i = 0; i < Vec_IntSize(vXors); ++i) {
        int v = Vec_IntEntry(vXors, i);
        if (v == 0) {
            int sz = Vec_IntSize(cur);
            if (sz >= 3) {
                int iXor = Vec_IntEntry(cur, 0);
                int in1  = Vec_IntEntry(cur, 1);
                int in2  = Vec_IntEntry(cur, 2);
                int in3  = (sz >= 4) ? Vec_IntEntry(cur, 3) : 0;

                // printf("  XOR-group: iXor=%d in1=%d in2=%d in3=%d (sz=%d)\n",
                //        iXor, in1, in2, in3, sz);

                if (in3 == 0) {
                    // XOR2
                    Vec_Int_t *xor_record = Vec_IntAlloc(10);
                    Vec_Int_t *vInput     = Vec_IntAlloc(2);

                    Vec_IntPush(xor_record, iXor);
                    Vec_IntPushTwo(vInput, in1, in2);
                    recursiveXor(p, vInput, iXor, xor_record);

                    Vec_IntAppend(res.xor_all_list, xor_record);
                    Vec_IntPush(res.xor_all_list, 0);

                    Vec_IntFree(xor_record);
                    Vec_IntFree(vInput);
                } else {
                    // XOR3
                    Vec_Int_t *xor_record = Vec_IntAlloc(10);
                    Vec_Int_t *vInput     = Vec_IntAlloc(3);

                    Vec_IntPush(xor_record, iXor);
                    Vec_IntPushThree(vInput, in1, in2, in3);
                    recursiveXor(p, vInput, iXor, xor_record);

                    Vec_IntAppend(res.xor_all_list, xor_record);
                    Vec_IntPush(res.xor_all_list, 0);

                    Vec_IntFree(xor_record);
                    Vec_IntFree(vInput);
                }
            } else if (sz > 0) {
                printf("[vXors] WARNING: too short group (sz=%d): ", sz);
                for (int t = 0; t < sz; ++t)
                    printf("%d ", Vec_IntEntry(cur, t));
                printf("\n");
            }
            Vec_IntClear(cur);
        } else {
            Vec_IntPush(cur, v);
        }
    }
    // Handle last group without trailing 0, if any
    if (Vec_IntSize(cur) >= 3) {
        int sz   = Vec_IntSize(cur);
        int iXor = Vec_IntEntry(cur, 0);
        int in1  = Vec_IntEntry(cur, 1);
        int in2  = Vec_IntEntry(cur, 2);
        int in3  = (sz >= 4) ? Vec_IntEntry(cur, 3) : 0;

        // printf("  XOR-group(last, no 0): iXor=%d in1=%d in2=%d in3=%d (sz=%d)\n",
        //        iXor, in1, in2, in3, sz);

        if (in3 == 0) {
            Vec_Int_t *xor_record = Vec_IntAlloc(10);
            Vec_Int_t *vInput     = Vec_IntAlloc(2);
            Vec_IntPush(xor_record, iXor);
            Vec_IntPushTwo(vInput, in1, in2);
            recursiveXor(p, vInput, iXor, xor_record);
            Vec_IntAppend(res.xor_all_list, xor_record);
            Vec_IntPush(res.xor_all_list, 0);
            Vec_IntFree(xor_record);
            Vec_IntFree(vInput);
        } else {
            Vec_Int_t *xor_record = Vec_IntAlloc(10);
            Vec_Int_t *vInput     = Vec_IntAlloc(3);
            Vec_IntPush(xor_record, iXor);
            Vec_IntPushThree(vInput, in1, in2, in3);
            recursiveXor(p, vInput, iXor, xor_record);
            Vec_IntAppend(res.xor_all_list, xor_record);
            Vec_IntPush(res.xor_all_list, 0);
            Vec_IntFree(xor_record);
            Vec_IntFree(vInput);
        }
    }
    Vec_IntFree(cur);

    // 3) xor_remaining = xor_all_list - groups that are subsets of HA/FA xor groups
    Vec_Int_t *current_group = Vec_IntAlloc(10);
    for (int i = 0; i < Vec_IntSize(res.xor_all_list); i++) {
        int v = Vec_IntEntry(res.xor_all_list, i);
        if (v == 0) {
            if (!is_group_subset_of_list(current_group, res.xor2_list) &&
                !is_group_subset_of_list(current_group, res.xor3_list))
            {
                Vec_IntAppend(res.xor_remaining, current_group);
                Vec_IntPush(res.xor_remaining, 0);
            }
            Vec_IntClear(current_group);
        } else {
            Vec_IntPush(current_group, v);
        }
    }
    Vec_IntFree(current_group);

    return res;
}

// ======================================================================
// Phase/NPN records for HA/FA
// ======================================================================
typedef struct {
    int isHA;               // 1=HA, 0=FA
    int in1, in2, in3;      // inputs; for HA, in3==0
    int xRoot, mRoot;       // root node IDs
    u8  tx;                 // XOR truth (for 2-vars: low 4 bits)
    u8  tm;                 // MAJ/AND truth
    u8  inMask;             // input phase mask (HA: 2 bits; FA: 3 bits)
    u8  outXor;             // output phase XOR (0/1)
    u8  outMaj;             // output phase MAJ (0/1)
    u8  okXor;              // XOR side matched? (0/1)
    u8  okMaj;              // MAJ side matched? (0/1)
} NpnPhaseRec;

typedef struct {
    NpnPhaseRec* data;
    int size, cap;
} NpnPhaseVec;

static void NpvInit(NpnPhaseVec* V) { V->data=NULL; V->size=0; V->cap=0; }
static void NpvFree(NpnPhaseVec* V) { if (V->data) free(V->data); V->data=NULL; V->size=0; V->cap=0; }
static void NpvEnsure(NpnPhaseVec* V, int need) {
    if (V->cap >= need) return;
    int cap2 = V->cap ? (V->cap<<1) : 64;
    if (cap2 < need) cap2 = need;
    V->data = (NpnPhaseRec*)realloc(V->data, sizeof(NpnPhaseRec)*cap2);
    V->cap = cap2;
}
static void NpvPush(NpnPhaseVec* V, const NpnPhaseRec* r) {
    NpvEnsure(V, V->size+1);
    V->data[V->size++] = *r;
}

// Phase dump container to separate HA/FA for JSON
typedef struct {
    NpnPhaseVec ha; // only records with isHA=1
    NpnPhaseVec fa; // only records with isHA=0
} PhaseDump;

static void PhaseDumpInit(PhaseDump* P) { NpvInit(&P->ha); NpvInit(&P->fa); }
static void PhaseDumpFree(PhaseDump* P) { NpvFree(&P->ha); NpvFree(&P->fa); }

// ======================================================================
// Truth-table / NPN helpers
// ======================================================================
static inline word TtBaseVar(int i)      { return s_Truths6[i]; }
static inline word TtBaseXor2()          { return TtBaseVar(0) ^ TtBaseVar(1); }
static inline word TtBaseAnd2()          { return TtBaseVar(0) & TtBaseVar(1); }
static inline word TtBaseXor3()          { return TtBaseVar(0) ^ TtBaseVar(1) ^ TtBaseVar(2); }
static inline word TtBaseMaj3()          { word A=TtBaseVar(0),B=TtBaseVar(1),C=TtBaseVar(2); return (A&B)|(A&C)|(B&C); }
static inline word TtFlipVar(word t,int v){ return Abc_Tt6Flip(t,v); }
static inline word TtSwap(word t,int i,int j){ return Abc_Tt6SwapVars(t,i,j); }

static word TtApplyPerm(word t, int nVars, const int *perm) {
    int cur[6] = {0,1,2,3,4,5};
    word w = t;
    for (int i = 0; i < nVars; ++i) {
        int pos = i;
        for (int j = i; j < nVars; ++j)
            if (cur[j] == perm[i]) { pos = j; break; }
        if (pos != i) {
            w = TtSwap(w, i, pos);
            int tmp = cur[i]; cur[i]=cur[pos]; cur[pos]=tmp;
        }
    }
    return w;
}

static inline word TtApplyInPhase(word t, int nVars, unsigned mask) {
    for (int i = 0; i < nVars; ++i)
        if (mask & (1u << i))
            t = TtFlipVar(t, i);
    return t;
}

static inline u8 TruncTruthToNV(word w, int nVars) {
    u8 t = (u8)(w & 0xFF);
    return (nVars <= 2) ? (t & 0x0F) : t;
}

// apply perm+phase on 8-bit truth for nVars<=3
static inline u8 ApplyPermPhase8(u8 tt, int nVars, const int* perm, unsigned inMask) {
    u8 out = 0;
    int max = 1 << nVars;
    for (int idx = 0; idx < max; ++idx) {
        int x[3] = {0,0,0}, y[3] = {0,0,0};
        for (int v=0; v<nVars; ++v) x[v] = (idx >> v) & 1;
        for (int v=0; v<nVars; ++v) {
            int bit = x[v];
            if (inMask & (1u << v)) bit ^= 1;
            y[perm[v]] = bit;
        }
        int src = 0;
        for (int v=0; v<nVars; ++v) src |= (y[v] << v);
        out |= (((tt >> src) & 1) << idx);
    }
    return out;
}

static int GenPerms1(int perms[1][1]) { perms[0][0]=0; return 1; }
static int GenPerms2(int perms[2][2]) {
    perms[0][0]=0; perms[0][1]=1;
    perms[1][0]=1; perms[1][1]=0;
    return 2;
}
static int GenPerms3(int perms[6][3]) {
    int P[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    memcpy(perms, P, sizeof(P));
    return 6;
}
static inline void InvertPerm(const int* perm, int n, int* inv) {
    for(int i=0;i<n;++i) inv[perm[i]] = i;
}

static void PrintTruthN(word w, const char* name, int nVars) {
    (void)w; (void)name; (void)nVars;
    // u8 b = TruncTruthToNV(w, nVars);
    // int bits = 1 << nVars;
    // printf("[base] %s tt=0x%02X bits=", name, (int)b);
    // for (int i = 0; i < bits; ++i)
    //     printf("%d", (b >> i) & 1);
    // printf("\n");
}
static int s_printedBases = 0;
static void PrintBaseTruthsOnce() {
    if (s_printedBases) return;
    s_printedBases = 1;
    PrintTruthN(TtBaseXor2(), "XOR2", 2);
    PrintTruthN(TtBaseAnd2(), "AND2", 2);
    PrintTruthN(TtBaseXor3(), "XOR3", 3);
    PrintTruthN(TtBaseMaj3(), "MAJ3", 3);
}

static inline void DecodePermByKey(int nVars, int permKey, int* permOut) {
    if (nVars == 2) {
        int p2[2][2]; GenPerms2(p2);
        permOut[0]=p2[permKey][0]; permOut[1]=p2[permKey][1];
        return;
    }
    int p3[6][3]; GenPerms3(p3);
    permOut[0]=p3[permKey][0];
    permOut[1]=p3[permKey][1];
    permOut[2]=p3[permKey][2];
}

// ======================================================================
// NPN matching solution sets for XOR/MAJ
// ======================================================================
typedef struct { unsigned key; u8 outBits; } SolKV; // key encodes (perm, inMask); outBits: bit0 => outPhase 0, bit1 => outPhase 1
typedef struct { SolKV* data; int size, cap; } SolVec;

static void SolVecInit(SolVec* S)  { S->data=NULL; S->size=0; S->cap=0; }
static void SolVecFree(SolVec* S)  { if(S->data) free(S->data); S->data=NULL; S->size=0; S->cap=0; }
static void SolVecEnsure(SolVec* S, int need){
    if(S->cap >= need) return;
    int cap2 = S->cap ? (S->cap<<1) : 32;
    if(cap2 < need) cap2=need;
    S->data = (SolKV*)realloc(S->data, sizeof(SolKV)*cap2);
    S->cap  = cap2;
}
static int SolVecFind(SolVec* S, unsigned key){
    for(int i=0;i<S->size;++i)
        if(S->data[i].key == key) return i;
    return -1;
}
static void SolVecAdd(SolVec* S, unsigned key, u8 outBit){
    int idx = SolVecFind(S, key);
    if(idx < 0){
        SolVecEnsure(S, S->size+1);
        S->data[S->size].key = key;
        S->data[S->size].outBits = 0;
        idx = S->size++;
    }
    S->data[idx].outBits |= (1u << outBit);
}

// Enumerate all (perm, input-phase, output-phase) solutions that map tt_in
// to the reference XOR or MAJ function (up to NPN).
static void EnumMatchNPN_Collect(u8 tt_in, int nVars, int isMaj, SolVec* outSols, int* outNum)
{
    *outNum = 0;
    SolVecInit(outSols);

    u8 majRef = (nVars==3) ? (u8)(TtBaseMaj3() & 0xFF)
                : (nVars==2) ? (u8)(TtBaseAnd2() & 0x0F)
                : 0;
    u8 xorRef = (nVars==3) ? (u8)(TtBaseXor3() & 0xFF)
                : (nVars==2) ? (u8)(TtBaseXor2() & 0x0F)
                : 0;
    u8 ref = isMaj ? majRef : xorRef;

    int permsCount = (nVars==2) ? 2 : (nVars==3 ? 6 : 1);
    unsigned inMaskMax = (nVars==2) ? 4 : (nVars==3 ? 8 : 1);

    int perms2[2][2]; GenPerms2(perms2);
    int perms3[6][3]; GenPerms3(perms3);

    for (int pk = 0; pk < permsCount; ++pk) {
        const int* perm = (nVars==2) ? perms2[pk] : perms3[pk];
        for (unsigned mask = 0; mask < inMaskMax; ++mask) {
            u8 t = ApplyPermPhase8(tt_in, nVars, perm, mask);
            unsigned key = (nVars==2) ? ((pk<<2) | mask) : ((pk<<3) | mask);
            if (t == ref) {
                SolVecAdd(outSols, key, 0); (*outNum)++;
            } else {
                u8 maskNV = (nVars==2) ? 0x0F : 0xFF;
                if (((u8)((~t) & maskNV)) == ref) {
                    SolVecAdd(outSols, key, 1); (*outNum)++;
                }
            }
        }
    }
}

// ======================================================================
// GIA helpers
// ======================================================================
static const char* Gia_ObjTypeStr(Gia_Man_t* p, Gia_Obj_t* pObj) {
    if (!pObj)                 return "null";
    if (Gia_ObjIsConst0(pObj)) return "const0";
    if (Gia_ObjIsCi(pObj))     return "ci";
    if (Gia_ObjIsCo(pObj))     return "co";
    if (Gia_ObjIsAnd(pObj))    return "and";
    return "other";
}

static int Gia_DfsCheckConeAndOnly_rec(Gia_Man_t* p, int iObj, Vec_Int_t* vLeaves, int rootId)
{
    if (Gia_ObjIsTravIdCurrentId(p, iObj))
        return 1;
    Gia_ObjSetTravIdCurrentId(p, iObj);

    Gia_Obj_t* pObj = Gia_ManObj(p, iObj);
    if (!pObj) {
        if (g_acecVerbose >= 2)
            printf("[giaTruth] ERROR: null object id=%d (root=%d)\n", iObj, rootId);
        return 0;
    }
    if (Gia_ObjIsConst0(pObj))
        return 1;
    if (Vec_IntFind(vLeaves, iObj) >= 0)
        return 1;

    if (!Gia_ObjIsAnd(pObj)) {
        if (g_acecVerbose >= 2) {
            const char* typ = Gia_ObjTypeStr(p, pObj);
            int lit0=-1,lit1=-1,fan0=-1,fan1=-1;
            if (Gia_ObjIsCo(pObj)) {
                int fan = Gia_ObjFaninId0p(p, pObj);
                printf("[giaTruth] Non-AND in cone: root=%d obj=%d type=%s fanin=%d\n",
                       rootId, iObj, typ, fan);
            } else {
                if (Gia_ObjFaninId0p(p, pObj) >= 0) {
                    lit0 = Gia_ObjFaninLit0p(p, pObj);
                    fan0 = Abc_Lit2Var(lit0);
                }
                if (Gia_ObjFaninId1p(p, pObj) >= 0) {
                    lit1 = Gia_ObjFaninLit1p(p, pObj);
                    fan1 = Abc_Lit2Var(lit1);
                }
                printf("[giaTruth] Non-AND in cone: root=%d obj=%d type=%s fanins: ",
                       rootId, iObj, typ);
                if (fan0 >= 0) printf("f0=%d%s ", fan0, Abc_LitIsCompl(lit0)?"'":"");
                if (fan1 >= 0) printf("f1=%d%s", fan1, Abc_LitIsCompl(lit1)?"'":"");
                printf("\n");
            }
            printf("[giaTruth] Abort truth due to non-AND node above.\n");
        }
        return 0;
    }
    int fan0 = Abc_Lit2Var(Gia_ObjFaninLit0p(p, pObj));
    int fan1 = Abc_Lit2Var(Gia_ObjFaninLit1p(p, pObj));
    if (!Gia_DfsCheckConeAndOnly_rec(p, fan0, vLeaves, rootId))
        return 0;
    if (!Gia_DfsCheckConeAndOnly_rec(p, fan1, vLeaves, rootId))
        return 0;
    return 1;
}

static int Gia_DfsCheckConeAndOnly(Gia_Man_t* p, int rootId, Vec_Int_t* vLeaves, const char* role)
{
    Gia_Obj_t* pRoot = Gia_ManObj(p, rootId);
    if (!pRoot) {
        if (g_acecVerbose >= 2)
            printf("[giaTruth] %s root %d invalid (null)\n", role, rootId);
        return 0;
    }
    if (Gia_ObjIsConst0(pRoot) || Vec_IntFind(vLeaves, rootId) >= 0)
        return 1;

    Gia_ManIncrementTravId(p);
    int ok = Gia_DfsCheckConeAndOnly_rec(p, rootId, vLeaves, rootId);
    if (!ok) {
        if (g_acecVerbose >= 2)
            printf("[giaTruth] Abort %s truth compute for root=%d due to non-AND node above.\n",
                   role, rootId);
        return 0;
    }
    return 1;
}

// ======================================================================
// Parsed adder list from combined groups
// ======================================================================
typedef struct {
    int isHA;    // 1=HA, 0=FA
    int in1,in2,in3;
    int xRoot,mRoot;
} ParsedAdderItem;

typedef struct {
    ParsedAdderItem* items;
    int size, cap;
} ParsedAdderList;

static void PAL_Init(ParsedAdderList* L) {
    L->items=NULL; L->size=0; L->cap=0;
}
static void PAL_Push(ParsedAdderList* L, const ParsedAdderItem* it) {
    if (L->size == L->cap) {
        int cap2 = L->cap ? (L->cap<<1) : 128;
        L->items = (ParsedAdderItem*)realloc(L->items,sizeof(ParsedAdderItem)*cap2);
        L->cap = cap2;
    }
    L->items[L->size++] = *it;
}
static void PAL_Free(ParsedAdderList* L) {
    if (L->items) free(L->items);
    L->items=NULL; L->size=0; L->cap=0;
}

// Parse HA and FA combined groups into a structured list
//   xor2_maj2_combined: [in1, in2, mid..., iXor, iMaj, 0]
//   xor3_maj3_combined: [in1, in2, in3, mid..., iXor, iMaj, 0]
static void BuildAdderListFromExtracted(
    const Vec_Int_t* vHA,
    const Vec_Int_t* vFA,
    ParsedAdderList* outList
) {
    PAL_Init(outList);
    // HA groups
    if (vHA) {
        int n = Vec_IntSize((Vec_Int_t*)vHA), i = 0;
        while (i < n) {
            int in1 = Vec_IntEntry((Vec_Int_t*)vHA, i++); if (i>=n) break;
            int in2 = Vec_IntEntry((Vec_Int_t*)vHA, i++); if (i>=n) break;
            int prev=-1, prev2=-1, j=i;
            for (; j<n; ++j) {
                int v = Vec_IntEntry((Vec_Int_t*)vHA, j);
                if (v == 0) {
                    if (in1>0 && in2>0 && prev2>0 && prev>0) {
                        ParsedAdderItem it;
                        it.isHA = 1;
                        it.in1  = in1;
                        it.in2  = in2;
                        it.in3  = 0;
                        it.xRoot= prev2;
                        it.mRoot= prev;
                        PAL_Push(outList, &it);
                    }
                    j++;
                    break;
                }
                prev2 = prev; prev = v;
            }
            i = j;
        }
    }
    // FA groups
    if (vFA) {
        int n = Vec_IntSize((Vec_Int_t*)vFA), i = 0;
        while (i < n) {
            int in1 = Vec_IntEntry((Vec_Int_t*)vFA, i++); if (i>=n) break;
            int in2 = Vec_IntEntry((Vec_Int_t*)vFA, i++); if (i>=n) break;
            int in3 = Vec_IntEntry((Vec_Int_t*)vFA, i++); if (i>=n) break;
            int prev=-1, prev2=-1, j=i;
            for (; j<n; ++j) {
                int v = Vec_IntEntry((Vec_Int_t*)vFA, j);
                if (v == 0) {
                    if (in1>0 && in2>0 && in3>0 && prev2>0 && prev>0) {
                        ParsedAdderItem it;
                        it.isHA = 0;
                        it.in1  = in1;
                        it.in2  = in2;
                        it.in3  = in3;
                        it.xRoot= prev2;
                        it.mRoot= prev;
                        PAL_Push(outList, &it);
                    }
                    j++;
                    break;
                }
                prev2 = prev; prev = v;
            }
            i = j;
        }
    }
}

// ======================================================================
// NPN checker for HA/FA from extracted list
// ======================================================================
static int* Gia_CheckNpnAdder_FromExtracted(
    Gia_Man_t* p,
    const ParsedAdderList* L,
    NpnPhaseVec* outPhases
) {
    // PrintBaseTruthsOnce();
    if (!p || !L) {
        if (g_acecVerbose >= 2)
            printf("Invalid args.\n");
        return NULL;
    }
    if (p->vTtMemory == NULL) {
        if (g_acecVerbose >= 2)
            printf("TT memory not initialized.\n");
        return NULL;
    }
    NpvInit(outPhases);

    int* passFlags = (int*)calloc(L->size,sizeof(int));
    if (!passFlags) {
        if (g_acecVerbose >= 2)
            printf("OOM: passFlags\n");
        return NULL;
    }

    for (int idx = 0; idx < L->size; ++idx) {
        const ParsedAdderItem* A = &L->items[idx];
        int isHA  = A->isHA ? 1 : 0;
        int nVars = isHA ? 2 : 3;

        Vec_Int_t* vLeaves = Vec_IntAlloc(3);
        Vec_IntPush(vLeaves, A->in1);
        Vec_IntPush(vLeaves, A->in2);
        if (!isHA)
            Vec_IntPush(vLeaves, A->in3);

        Gia_Obj_t* pX = Gia_ManObj(p, A->xRoot);
        Gia_Obj_t* pM = Gia_ManObj(p, A->mRoot);
        if (!pX || !pM) {
            if (g_acecVerbose >= 2)
                printf("[%s] invalid roots x=%d m=%d\n", isHA?"HA":"FA", A->xRoot, A->mRoot);
            Vec_IntFree(vLeaves);
            NpnPhaseRec rec = {0};
            rec.isHA = isHA;
            rec.in1  = A->in1; rec.in2=A->in2; rec.in3=A->in3;
            rec.xRoot=A->xRoot; rec.mRoot=A->mRoot;
            NpvPush(outPhases, &rec);
            continue;
        }

        if (!Gia_DfsCheckConeAndOnly(p, A->xRoot, vLeaves, "xor")) {
            Vec_IntFree(vLeaves);
            NpnPhaseRec rec = {0};
            rec.isHA = isHA;
            rec.in1  = A->in1; rec.in2=A->in2; rec.in3=A->in3;
            rec.xRoot=A->xRoot; rec.mRoot=A->mRoot;
            NpvPush(outPhases, &rec);
            continue;
        }
        if (!Gia_DfsCheckConeAndOnly(p, A->mRoot, vLeaves, "maj")) {
            Vec_IntFree(vLeaves);
            NpnPhaseRec rec = {0};
            rec.isHA = isHA;
            rec.in1  = A->in1; rec.in2=A->in2; rec.in3=A->in3;
            rec.xRoot=A->xRoot; rec.mRoot=A->mRoot;
            NpvPush(outPhases, &rec);
            continue;
        }

        word* tx = Gia_ObjComputeTruthTableCut(p, pX, vLeaves);
        word* tm = Gia_ObjComputeTruthTableCut(p, pM, vLeaves);
        if (!tx || !tm) {
            if (g_acecVerbose >= 2)
                printf("[%s] TT compute failed\n", isHA?"HA":"FA");
            Vec_IntFree(vLeaves);
            NpnPhaseRec rec = {0};
            rec.isHA = isHA;
            rec.in1  = A->in1; rec.in2=A->in2; rec.in3=A->in3;
            rec.xRoot=A->xRoot; rec.mRoot=A->mRoot;
            NpvPush(outPhases, &rec);
            continue;
        }

        u8 tx8 = TruncTruthToNV(tx[0], nVars);
        u8 tm8 = TruncTruthToNV(tm[0], nVars);

        // if (isHA) {
        //     printf("[HA] leaves=(%d,%d) roots: xor=%d maj=%d\n",
        //            A->in1, A->in2, A->xRoot, A->mRoot);
        // } else {
        //     printf("[FA] leaves=(%d,%d,%d) roots: xor=%d maj=%d\n",
        //            A->in1, A->in2, A->in3, A->xRoot, A->mRoot);
        // }
        // printf("  raw truth: [maj=0x%02X xor=0x%02X]\n", (int)tm8, (int)tx8);

        SolVec majSols, xorSols; int nMaj=0, nXor=0;
        EnumMatchNPN_Collect(tm8, nVars, 1, &majSols, &nMaj);
        EnumMatchNPN_Collect(tx8, nVars, 0, &xorSols, &nXor);

        int interCount = 0;
        int haveBest   = 0;
        unsigned bestKey = 0;
        u8 bestOutMaj = 0, bestOutXor = 0;

        for (int i=0; i<majSols.size; ++i) {
            unsigned key = majSols.data[i].key;
            int j = SolVecFind(&xorSols, key);
            if (j < 0) continue;
            interCount++;

            u8 mBits = majSols.data[i].outBits;
            u8 xBits = xorSols.data[j].outBits;
            for (int om=0; om<=1; ++om) if (mBits & (1u<<om))
            for (int ox=0; ox<=1; ++ox) if (xBits & (1u<<ox)) {
                if (!haveBest) {
                    bestKey   = key;
                    bestOutMaj= (u8)om;
                    bestOutXor= (u8)ox;
                    haveBest  = 1;
                } else {
                    int pkA,imA,pkB,imB;
                    if (nVars==2) {
                        pkA = (bestKey>>2); imA=(bestKey&3);
                        pkB = (int)(key>>2); imB=(key&3);
                    } else {
                        pkA = (bestKey>>3); imA=(bestKey&7);
                        pkB = (int)(key>>3); imB=(key&7);
                    }
                    int better = (pkB<pkA) ||
                                 (pkB==pkA && (imB<imA ||
                                  (imB==imA && ((bestOutMaj && !om) ||
                                   (bestOutMaj==om && bestOutXor && !ox)))));
                    if (better) {
                        bestKey   = key;
                        bestOutMaj= (u8)om;
                        bestOutXor= (u8)ox;
                    }
                }
            }
        }

        // printf("  matches: maj=%d xor=%d inter=%d\n", nMaj, nXor, interCount);

        int okM=0, okX=0;
        unsigned inMask = 0;
        if (!haveBest) {
            // printf("  NO COMMON (perm,inPhase) -> reject\n");
        } else {
            int permKey, inMaskTmp;
            if (nVars==2){
                permKey  = (int)(bestKey>>2);
                inMaskTmp= (int)(bestKey&3);
            } else {
                permKey  = (int)(bestKey>>3);
                inMaskTmp= (int)(bestKey&7);
            }
            int perm[3]={0,1,2},invPerm[3]={0,1,2};
            DecodePermByKey(nVars, permKey, perm);
            InvertPerm(perm, nVars, invPerm);

            word basisM = (nVars==3)? TtBaseMaj3() : TtBaseAnd2();
            word basisX = (nVars==3)? TtBaseXor3() : TtBaseXor2();

            word mRaw = TtApplyInPhase(TtApplyPerm(basisM, nVars, invPerm), nVars, (unsigned)inMaskTmp);
            word xRaw = TtApplyInPhase(TtApplyPerm(basisX, nVars, invPerm), nVars, (unsigned)inMaskTmp);

            u8 maskNV = (nVars==2)? 0x0F : 0xFF;
            u8 majExp = (u8)((bestOutMaj ? (~mRaw) : mRaw) & maskNV);
            u8 xorExp = (u8)((bestOutXor ? (~xRaw) : xRaw) & maskNV);

            okM    = (majExp == tm8);
            okX    = (xorExp == tx8);
            inMask = (unsigned)inMaskTmp;

            // printf("  perm=[");
            // for (int v=0; v<nVars; ++v) {
            //     if (v) printf(",");
            //     printf("%d", perm[v]);
            // }
            // printf("] inPhase=");
            // if (nVars==2)
            //     printf("[%u,%u]\n", (inMask&1)?1u:0u, (inMask&2)?1u:0u);
            // else
            //     printf("[%u,%u,%u]\n",
            //            (inMask&1)?1u:0u, (inMask&2)?1u:0u, (inMask&4)?1u:0u);
            // printf("  outPhase: [xor=%u maj=%u]\n", (unsigned)bestOutXor, (unsigned)bestOutMaj);
            // printf("  check: MAJ%s XOR%s\n", okM?" OK":" BAD", okX?" OK":" BAD");
        }

        if (okM && okX)
            passFlags[idx] = 1;

        // Always record a phase record for JSON (ok flags indicate validity)
        NpnPhaseRec rec = {0};
        rec.isHA   = isHA;
        rec.in1    = A->in1;
        rec.in2    = A->in2;
        rec.in3    = A->in3;
        rec.xRoot  = A->xRoot;
        rec.mRoot  = A->mRoot;
        rec.tx     = tx8;
        rec.tm     = tm8;
        rec.inMask = (u8)inMask;
        rec.outXor = haveBest ? bestOutXor : 0;
        rec.outMaj = haveBest ? bestOutMaj : 0;
        rec.okXor  = okX ? 1 : 0;
        rec.okMaj  = okM ? 1 : 0;
        NpvPush(outPhases, &rec);

        SolVecFree(&majSols);
        SolVecFree(&xorSols);
        Vec_IntFree(vLeaves);
    }
    return passFlags;
}

// ======================================================================
// Multi-threaded version of NPN checking for HA/FA
// Note: Truth table computation is NOT thread-safe (uses shared p->vTtMemory)
// So we compute truth tables serially first, then parallelize NPN matching
// ======================================================================
#ifdef ABC_USE_PTHREADS

// Pre-computed truth table data for one adder
typedef struct {
    int idx;
    int isHA;
    int nVars;
    int in1, in2, in3;
    int xRoot, mRoot;
    u8 tx8;  // XOR truth table
    u8 tm8;  // MAJ truth table
    int valid;  // 1 if truth tables were computed successfully
} TruthTableData;

// Thread-local result for one adder
typedef struct {
    int idx;
    int passFlag;
    NpnPhaseRec rec;
} NpnThreadResult;

// Thread argument structure for NPN matching
typedef struct {
    const TruthTableData* ttData;
    int startIdx;
    int endIdx;
    NpnThreadResult* results;
} NpnMatchThreadArg;

// Worker function for NPN matching (thread-safe, no shared state)
static void* Gia_NpnMatch_Worker(void* arg)
{
    NpnMatchThreadArg* targ = (NpnMatchThreadArg*)arg;

    for (int idx = targ->startIdx; idx < targ->endIdx; ++idx) {
        const TruthTableData* tt = &targ->ttData[idx];
        NpnThreadResult* res = &targ->results[idx];
        res->idx = idx;
        res->passFlag = 0;
        memset(&res->rec, 0, sizeof(NpnPhaseRec));

        res->rec.isHA  = tt->isHA;
        res->rec.in1   = tt->in1;
        res->rec.in2   = tt->in2;
        res->rec.in3   = tt->in3;
        res->rec.xRoot = tt->xRoot;
        res->rec.mRoot = tt->mRoot;

        if (!tt->valid) {
            continue;
        }

        int nVars = tt->nVars;
        u8 tx8 = tt->tx8;
        u8 tm8 = tt->tm8;

        res->rec.tx = tx8;
        res->rec.tm = tm8;

        SolVec majSols, xorSols; int nMaj=0, nXor=0;
        EnumMatchNPN_Collect(tm8, nVars, 1, &majSols, &nMaj);
        EnumMatchNPN_Collect(tx8, nVars, 0, &xorSols, &nXor);

        int haveBest   = 0;
        unsigned bestKey = 0;
        u8 bestOutMaj = 0, bestOutXor = 0;

        for (int i=0; i<majSols.size; ++i) {
            unsigned key = majSols.data[i].key;
            int j = SolVecFind(&xorSols, key);
            if (j < 0) continue;

            u8 mBits = majSols.data[i].outBits;
            u8 xBits = xorSols.data[j].outBits;
            for (int om=0; om<=1; ++om) if (mBits & (1u<<om))
            for (int ox=0; ox<=1; ++ox) if (xBits & (1u<<ox)) {
                if (!haveBest) {
                    bestKey   = key;
                    bestOutMaj= (u8)om;
                    bestOutXor= (u8)ox;
                    haveBest  = 1;
                } else {
                    int pkA,imA,pkB,imB;
                    if (nVars==2) {
                        pkA = (bestKey>>2); imA=(bestKey&3);
                        pkB = (int)(key>>2); imB=(key&3);
                    } else {
                        pkA = (bestKey>>3); imA=(bestKey&7);
                        pkB = (int)(key>>3); imB=(key&7);
                    }
                    int better = (pkB<pkA) ||
                                 (pkB==pkA && (imB<imA ||
                                  (imB==imA && ((bestOutMaj && !om) ||
                                   (bestOutMaj==om && bestOutXor && !ox)))));
                    if (better) {
                        bestKey   = key;
                        bestOutMaj= (u8)om;
                        bestOutXor= (u8)ox;
                    }
                }
            }
        }

        int okM=0, okX=0;
        unsigned inMask = 0;
        if (haveBest) {
            int permKey, inMaskTmp;
            if (nVars==2){
                permKey  = (int)(bestKey>>2);
                inMaskTmp= (int)(bestKey&3);
            } else {
                permKey  = (int)(bestKey>>3);
                inMaskTmp= (int)(bestKey&7);
            }
            int perm[3]={0,1,2},invPerm[3]={0,1,2};
            DecodePermByKey(nVars, permKey, perm);
            InvertPerm(perm, nVars, invPerm);

            word basisM = (nVars==3)? TtBaseMaj3() : TtBaseAnd2();
            word basisX = (nVars==3)? TtBaseXor3() : TtBaseXor2();

            word mRaw = TtApplyInPhase(TtApplyPerm(basisM, nVars, invPerm), nVars, (unsigned)inMaskTmp);
            word xRaw = TtApplyInPhase(TtApplyPerm(basisX, nVars, invPerm), nVars, (unsigned)inMaskTmp);

            u8 maskNV = (nVars==2)? 0x0F : 0xFF;
            u8 majExp = (u8)((bestOutMaj ? (~mRaw) : mRaw) & maskNV);
            u8 xorExp = (u8)((bestOutXor ? (~xRaw) : xRaw) & maskNV);

            okM    = (majExp == tm8);
            okX    = (xorExp == tx8);
            inMask = (unsigned)inMaskTmp;
        }

        if (okM && okX)
            res->passFlag = 1;

        res->rec.inMask = (u8)inMask;
        res->rec.outXor = haveBest ? bestOutXor : 0;
        res->rec.outMaj = haveBest ? bestOutMaj : 0;
        res->rec.okXor  = okX ? 1 : 0;
        res->rec.okMaj  = okM ? 1 : 0;

        SolVecFree(&majSols);
        SolVecFree(&xorSols);
    }
    return NULL;
}

// Multi-threaded version of Gia_CheckNpnAdder_FromExtracted
// Phase 1: Compute truth tables serially (not thread-safe)
// Phase 2: NPN matching in parallel (thread-safe)
static int* Gia_CheckNpnAdder_FromExtracted_MT(
    Gia_Man_t* p,
    const ParsedAdderList* L,
    NpnPhaseVec* outPhases,
    int nThreads
) {
    if (!p || !L) {
        if (g_acecVerbose >= 2)
            printf("Invalid args.\n");
        return NULL;
    }
    if (p->vTtMemory == NULL) {
        if (g_acecVerbose >= 2)
            printf("TT memory not initialized.\n");
        return NULL;
    }
    NpvInit(outPhases);

    int* passFlags = (int*)calloc(L->size, sizeof(int));
    if (!passFlags) {
        if (g_acecVerbose >= 2)
            printf("OOM: passFlags\n");
        return NULL;
    }

    if (L->size == 0) {
        return passFlags;
    }

    // Phase 1: Compute truth tables serially
    TruthTableData* ttData = (TruthTableData*)calloc(L->size, sizeof(TruthTableData));
    if (!ttData) {
        free(passFlags);
        return NULL;
    }

    for (int idx = 0; idx < L->size; ++idx) {
        const ParsedAdderItem* A = &L->items[idx];
        TruthTableData* tt = &ttData[idx];

        tt->idx = idx;
        tt->isHA = A->isHA ? 1 : 0;
        tt->nVars = tt->isHA ? 2 : 3;
        tt->in1 = A->in1;
        tt->in2 = A->in2;
        tt->in3 = A->in3;
        tt->xRoot = A->xRoot;
        tt->mRoot = A->mRoot;
        tt->valid = 0;

        Vec_Int_t* vLeaves = Vec_IntAlloc(3);
        Vec_IntPush(vLeaves, A->in1);
        Vec_IntPush(vLeaves, A->in2);
        if (!tt->isHA)
            Vec_IntPush(vLeaves, A->in3);

        Gia_Obj_t* pX = Gia_ManObj(p, A->xRoot);
        Gia_Obj_t* pM = Gia_ManObj(p, A->mRoot);
        if (!pX || !pM) {
            Vec_IntFree(vLeaves);
            continue;
        }

        if (!Gia_DfsCheckConeAndOnly(p, A->xRoot, vLeaves, "xor") ||
            !Gia_DfsCheckConeAndOnly(p, A->mRoot, vLeaves, "maj")) {
            Vec_IntFree(vLeaves);
            continue;
        }

        word* tx = Gia_ObjComputeTruthTableCut(p, pX, vLeaves);
        word* tm = Gia_ObjComputeTruthTableCut(p, pM, vLeaves);
        if (tx && tm) {
            tt->tx8 = TruncTruthToNV(tx[0], tt->nVars);
            tt->tm8 = TruncTruthToNV(tm[0], tt->nVars);
            tt->valid = 1;
        }

        Vec_IntFree(vLeaves);
    }

    // Phase 2: NPN matching in parallel
    // Limit threads to available work
    if (nThreads > L->size) nThreads = L->size;
    if (nThreads < 1) nThreads = 1;

    NpnThreadResult* results = (NpnThreadResult*)calloc(L->size, sizeof(NpnThreadResult));
    if (!results) {
        free(ttData);
        free(passFlags);
        return NULL;
    }

    pthread_t* threads = (pthread_t*)malloc(nThreads * sizeof(pthread_t));
    NpnMatchThreadArg* args = (NpnMatchThreadArg*)malloc(nThreads * sizeof(NpnMatchThreadArg));

    int chunkSize = (L->size + nThreads - 1) / nThreads;
    for (int t = 0; t < nThreads; ++t) {
        args[t].ttData = ttData;
        args[t].startIdx = t * chunkSize;
        args[t].endIdx = (t + 1) * chunkSize;
        if (args[t].endIdx > L->size) args[t].endIdx = L->size;
        args[t].results = results;

        pthread_create(&threads[t], NULL, Gia_NpnMatch_Worker, &args[t]);
    }

    // Wait for all threads
    for (int t = 0; t < nThreads; ++t) {
        pthread_join(threads[t], NULL);
    }

    // Collect results
    for (int i = 0; i < L->size; ++i) {
        passFlags[i] = results[i].passFlag;
        NpvPush(outPhases, &results[i].rec);
    }

    free(threads);
    free(args);
    free(results);
    free(ttData);

    return passFlags;
}

#endif // ABC_USE_PTHREADS

// ======================================================================
// Filter combined HA/FA groups by passFlags
// ======================================================================
typedef struct {
    int isHA;
    int xRoot;
    int mRoot;
    int pass;
} HmPass;

// Filter combined = [in..., middle..., xRoot, mRoot, 0] by (isHA,xRoot,mRoot)->pass
static Vec_Int_t* FilterCombinedByPass(
    const Vec_Int_t* combined,
    const HmPass* map,
    int mapSize,
    int isHA
){
    Vec_Int_t* out = Vec_IntAlloc(Vec_IntSize((Vec_Int_t*)combined));
    int n = Vec_IntSize((Vec_Int_t*)combined);
    int i = 0;

    while (i < n) {
        int start = i;
        int in1 = Vec_IntEntry((Vec_Int_t*)combined, i++); if (i>=n) break;
        int in2 = Vec_IntEntry((Vec_Int_t*)combined, i++); if (i>=n) break;
        if (!isHA) {
            // FA has third input but we just skip reading it here for filtering
            int in3 = Vec_IntEntry((Vec_Int_t*)combined, i++);
            (void)in3;
        }

        int prev=-1, prev2=-1, j=i;
        for (; j<n; ++j) {
            int v = Vec_IntEntry((Vec_Int_t*)combined, j);
            if (v == 0) {
                int iXor = prev2;
                int iMaj = prev;
                int keep = 0;
                for (int t=0; t<mapSize; ++t) {
                    if (map[t].isHA == isHA &&
                        map[t].xRoot == iXor &&
                        map[t].mRoot == iMaj)
                    {
                        keep = map[t].pass;
                        break;
                    }
                }
                if (keep) {
                    for (int k=start; k<=j; ++k)
                        Vec_IntPush(out, Vec_IntEntry((Vec_Int_t*)combined, k));
                }
                j++;
                break;
            }
            prev2 = prev;
            prev  = v;
        }
        i = j;
    }
    return out;
}

// ======================================================================
// Group extraction from combined vectors (0-separated)
// ======================================================================
static Vec_Wec_t* groups_from_combined(const Vec_Int_t* combined) {
    Vec_Wec_t* W = Vec_WecStart(16);
    if (!combined)
        return W;

    int n = Vec_IntSize((Vec_Int_t*)combined);
    Vec_Int_t* cur = Vec_IntAlloc(32);
    for (int i=0; i<n; ++i) {
        int v = Vec_IntEntry((Vec_Int_t*)combined, i);
        if (v == 0) {
            Vec_Int_t* slot = Vec_WecPushLevel(W);
            Vec_IntAppend(slot, cur);
            Vec_IntClear(cur);
        } else {
            Vec_IntPush(cur, v);
        }
    }
    Vec_IntFree(cur);
    return W;
}

static const NpnPhaseRec* find_phase_rec(
    const NpnPhaseVec* V,
    int isHA,
    int xRoot,
    int mRoot
) {
    if (!V || !V->data) return NULL;
    for (int i=0; i<V->size; ++i) {
        const NpnPhaseRec* r = &V->data[i];
        if (r->isHA == isHA && r->xRoot == xRoot && r->mRoot == mRoot)
            return r;
    }
    return NULL;
}

// ======================================================================
// Phase records for XOR-remaining (support 2-input and 3-input XOR)
// ======================================================================
typedef struct {
    int  nVars;           // 2 or 3
    int  xorRoot;
    int  in1, in2, in3;   // leaves used for truth/NPN (for 2-input, in3=0)
    int  gLast1, gLast2;  // last two nodes of the group (for debugging/JSON)
    u8   tt;              // 4-bit truth if nVars=2, 8-bit truth if nVars=3
    u8   inMask;          // 2 or 3 bits of input phase mask (low bits)
    u8   outXor;          // output phase: 0/1
    u8   ok;              // NPN match success
} XorPhaseRec;

typedef struct {
    XorPhaseRec* data;
    int size, cap;
} XorPhaseVec;

static void XpvInit(XorPhaseVec* V){ V->data=NULL; V->size=0; V->cap=0; }
static void XpvFree(XorPhaseVec* V){ if(V->data) free(V->data); V->data=NULL; V->size=0; V->cap=0; }
static void XpvEnsure(XorPhaseVec* V, int need){
    if (V->cap >= need) return;
    int cap2 = V->cap ? (V->cap<<1) : 64;
    if (cap2 < need) cap2 = need;
    V->data = (XorPhaseRec*)realloc(V->data, sizeof(XorPhaseRec)*cap2);
    V->cap  = cap2;
}
static void XpvPush(XorPhaseVec* V, const XorPhaseRec* r){
    XpvEnsure(V, V->size+1);
    V->data[V->size++] = *r;
}

// Parsed XOR-remaining item (can be XOR2 or XOR3)
typedef struct {
    int nVars;            // 2 or 3
    int xorRoot;
    int in1, in2, in3;    // leaves for truth/NPN (for 2-input, in3=0)
    int gLast1, gLast2;   // last two nodes of the group
} ParsedXorItem;

typedef struct {
    ParsedXorItem* items;
    int size, cap;
} ParsedXorList;

static void PXL_Init(ParsedXorList* L) {
    L->items=NULL; L->size=0; L->cap=0;
}
static void PXL_Free(ParsedXorList* L) {
    if(L->items) free(L->items);
    L->items=NULL; L->size=0; L->cap=0;
}
static void PXL_Push(ParsedXorList* L, const ParsedXorItem* it) {
    if(L->size == L->cap) {
        int cap2 = L->cap ? (L->cap<<1) : 128;
        L->items = (ParsedXorItem*)realloc(L->items, sizeof(ParsedXorItem)*cap2);
        L->cap   = cap2;
    }
    L->items[L->size++] = *it;
}

// Build XOR list from xor_remaining groups.
//
// Rule:
//   - Let sz = number of nodes in the group (no trailing 0 here).
//   - If sz > 5 => treat this group as XOR3 (3-variable).
//   - Else       => treat as XOR2 (2-variable).
//
// For both cases:
//   - xorRoot = first node of the group
//   - gLast1/gLast2 = last two nodes of the group
//   - For XOR2: in1,in2 = gLast1,gLast2; in3 = 0
//   - For XOR3: in1,in2,in3 = the first 3 nodes after xorRoot
//              (i.e. entries 1,2,3 of the group, if available)
//
static void BuildXorListFromRemaining(
    Gia_Man_t* p,
    const Vec_Int_t* xor_remaining,
    ParsedXorList* out
){
    (void)p; // currently unused, kept for future extension
    PXL_Init(out);
    Vec_Wec_t* W = groups_from_combined(xor_remaining);

    for (int gi=0; gi<Vec_WecSize(W); ++gi) {
        Vec_Int_t* G = Vec_WecEntry(W, gi);
        int sz = Vec_IntSize(G);
        if (sz < 3) // need at least [xorRoot, x, y]
            continue;

        ParsedXorItem it;
        memset(&it, 0, sizeof(it));

        it.xorRoot = Vec_IntEntry(G, 0);
        it.gLast1  = Vec_IntEntry(G, sz - 2);
        it.gLast2  = Vec_IntEntry(G, sz - 1);

        if (sz > 5) {
            // Treat as XOR3
            it.nVars = 3;
            it.in1 = (sz > 1) ? Vec_IntEntry(G, 1) : 0;
            it.in2 = (sz > 2) ? Vec_IntEntry(G, 2) : 0;
            it.in3 = (sz > 3) ? Vec_IntEntry(G, 3) : 0;
        } else {
            // Treat as XOR2
            it.nVars = 2;
            it.in1 = it.gLast1;
            it.in2 = it.gLast2;
            it.in3 = 0;
        }

        if (it.in1 <= 0 || it.in2 <= 0)
            continue;

        PXL_Push(out, &it);
    }
    Vec_WecFree(W);
}

// Compute XOR NPN phase for xor_remaining list.
// If nVars=2, treat as XOR2 (4-bit truth).
// If nVars=3, treat as XOR3 (8-bit truth, same style as FA xor3).
static void CheckXorRemainingPhases(
    Gia_Man_t* p,
    const ParsedXorList* XL,
    XorPhaseVec* out
){
    XpvInit(out);
    if (!p || !XL) return;

    // PrintBaseTruthsOnce();
    if (p->vTtMemory == NULL) {
        printf("TT memory not initialized for xor_remaining.\n");
        return;
    }

    for (int i=0; i<XL->size; ++i) {
        const ParsedXorItem* X = &XL->items[i];

        XorPhaseRec rec;
        memset(&rec, 0, sizeof(rec));

        rec.nVars  = X->nVars;
        rec.xorRoot= X->xorRoot;
        rec.in1    = X->in1;
        rec.in2    = X->in2;
        rec.in3    = X->in3;
        rec.gLast1 = X->gLast1;
        rec.gLast2 = X->gLast2;

        int nVars = X->nVars;
        if (nVars != 2 && nVars != 3) {
            rec.ok = 0;
            XpvPush(out, &rec);
            continue;
        }

        if (X->xorRoot <= 0 || X->in1 <= 0 || X->in2 <= 0 ||
            (nVars == 3 && X->in3 <= 0))
        {
            rec.ok = 0;
            XpvPush(out, &rec);
            continue;
        }

        Vec_Int_t* vLeaves = Vec_IntAlloc(3);
        Vec_IntPush(vLeaves, X->in1);
        Vec_IntPush(vLeaves, X->in2);
        if (nVars == 3)
            Vec_IntPush(vLeaves, X->in3);

        if (!Gia_DfsCheckConeAndOnly(p, X->xorRoot, vLeaves, "xor_remaining")) {
            Vec_IntFree(vLeaves);
            rec.ok = 0;
            XpvPush(out, &rec);
            continue;
        }

        Gia_Obj_t* pX = Gia_ManObj(p, X->xorRoot);
        if (!pX) {
            Vec_IntFree(vLeaves);
            rec.ok = 0;
            XpvPush(out, &rec);
            continue;
        }

        word* tx = Gia_ObjComputeTruthTableCut(p, pX, vLeaves);
        if (!tx) {
            Vec_IntFree(vLeaves);
            rec.ok = 0;
            XpvPush(out, &rec);
            continue;
        }

        u8 tx8 = TruncTruthToNV(tx[0], nVars);
        rec.tt = tx8;

        SolVec xorSols; int nXor=0;
        EnumMatchNPN_Collect(tx8, nVars, 0, &xorSols, &nXor);

        int haveBest = 0;
        unsigned bestKey = 0;
        u8 bestOut = 0;

        for (int s=0; s<xorSols.size; ++s) {
            unsigned key   = xorSols.data[s].key;
            u8       bits  = xorSols.data[s].outBits;
            for (int ox=0; ox<=1; ++ox) if (bits & (1u<<ox)) {
                if (!haveBest) {
                    bestKey = key;
                    bestOut = (u8)ox;
                    haveBest= 1;
                } else {
                    int pkA, imA, pkB, imB;
                    if (nVars == 2) {
                        pkA = (bestKey>>2); imA=(bestKey&3);
                        pkB = (int)(key>>2); imB=(key&3);
                    } else {
                        pkA = (bestKey>>3); imA=(bestKey&7);
                        pkB = (int)(key>>3); imB=(key&7);
                    }
                    int better = (pkB<pkA) ||
                                 (pkB==pkA && (imB<imA ||
                                  (imB==imA && bestOut && !ox)));
                    if (better) {
                        bestKey = key;
                        bestOut = (u8)ox;
                    }
                }
            }
        }

        if (!haveBest) {
            rec.ok     = 0;
            rec.inMask = 0;
            rec.outXor = 0;
        } else {
            unsigned inMaskTmp;
            int permKey;

            if (nVars == 2) {
                permKey   = (int)(bestKey >> 2);
                inMaskTmp = (unsigned)(bestKey & 3);
            } else {
                permKey   = (int)(bestKey >> 3);
                inMaskTmp = (unsigned)(bestKey & 7);
            }

            rec.inMask = (u8)inMaskTmp;
            rec.outXor = bestOut;

            int perm[3]={0,1,2}, invPerm[3]={0,1,2};
            DecodePermByKey(nVars, permKey, perm);
            InvertPerm(perm, nVars, invPerm);

            word basisX = (nVars == 3) ? TtBaseXor3() : TtBaseXor2();
            word xRaw   = TtApplyInPhase(
                              TtApplyPerm(basisX, nVars, invPerm),
                              nVars,
                              (unsigned)inMaskTmp
                         );

            u8 maskNV = (nVars == 2) ? 0x0F : 0xFF;
            u8 xorExp = (u8)((bestOut ? (~xRaw) : xRaw) & maskNV);
            rec.ok = (xorExp == tx8) ? 1 : 0;
        }

        SolVecFree(&xorSols);
        Vec_IntFree(vLeaves);
        XpvPush(out, &rec);
    }
}

// ======================================================================
// Populate structured HA/FA/XOR output (replaces JSON writing)
// ======================================================================

// Helper to allocate and copy an int array
static int* copy_int_array(const int* src, int count) {
    if (count <= 0 || !src) return NULL;
    int* dst = (int*)malloc(sizeof(int) * count);
    if (dst) memcpy(dst, src, sizeof(int) * count);
    return dst;
}

// Helper to allocate int array from Vec_Int_t
static int* copy_vec_int_to_array(Vec_Int_t* v, int* out_count) {
    if (!v) { *out_count = 0; return NULL; }
    int n = Vec_IntSize(v);
    *out_count = n;
    if (n <= 0) return NULL;
    int* arr = (int*)malloc(sizeof(int) * n);
    if (arr) {
        for (int i = 0; i < n; i++) arr[i] = Vec_IntEntry(v, i);
    }
    return arr;
}

// Populate GiaHalfAdder array from extraction results
static void populate_half_adders(
    GiaDetectingResults* results,
    const VecIntResults1* er,
    const ParsedAdderList* LOrder,
    const PhaseDump* phases
) {
    // Count HAs from LOrder
    int ha_count = 0;
    for (int i = 0; i < LOrder->size; i++) {
        if (LOrder->items[i].isHA) ha_count++;
    }

    if (ha_count <= 0) {
        results->half_adders = NULL;
        results->half_adders_count = 0;
        return;
    }

    results->half_adders = (GiaHalfAdder*)calloc(ha_count, sizeof(GiaHalfAdder));
    results->half_adders_count = ha_count;

    if (!results->half_adders) {
        results->half_adders_count = 0;
        return;
    }

    // Build HA groups from combined
    Vec_Wec_t* ha_groupsW = groups_from_combined(er->xor2_maj2_combined);

    int id = 0;
    for (int i = 0; i < LOrder->size && id < ha_count; i++) {
        const ParsedAdderItem* A = &LOrder->items[i];
        if (!A->isHA) continue;

        GiaHalfAdder* ha = &results->half_adders[id];
        ha->id = id;

        // Find matching nodes group
        Vec_Int_t* nodes = NULL;
        for (int gw = 0; gw < Vec_WecSize(ha_groupsW); gw++) {
            Vec_Int_t* G = Vec_WecEntry(ha_groupsW, gw);
            int sz = Vec_IntSize(G);
            if (sz >= 4) {
                int xr = Vec_IntEntry(G, sz - 2);
                int mr = Vec_IntEntry(G, sz - 1);
                if (xr == A->xRoot && mr == A->mRoot) {
                    nodes = G;
                    break;
                }
            }
        }

        // Copy nodes
        if (nodes) {
            ha->nodes = copy_vec_int_to_array(nodes, &ha->nodes_count);
        }

        // Find phase record
        const NpnPhaseRec* r = find_phase_rec(&phases->ha, 1, A->xRoot, A->mRoot);
        if (r) {
            // tt: 2 bytes [tx, tm]
            ha->tt_count = 2;
            ha->tt = (unsigned char*)malloc(2);
            if (ha->tt) {
                ha->tt[0] = r->tx;
                ha->tt[1] = r->tm;
            }

            // input_phase: 2 elements
            ha->input_phase_count = 2;
            ha->input_phase = (int*)malloc(2 * sizeof(int));
            if (ha->input_phase) {
                ha->input_phase[0] = (r->inMask & 1) ? 1 : 0;
                ha->input_phase[1] = (r->inMask & 2) ? 1 : 0;
            }

            // output_phase: 2 elements [outXor, outMaj]
            ha->output_phase_count = 2;
            ha->output_phase = (int*)malloc(2 * sizeof(int));
            if (ha->output_phase) {
                ha->output_phase[0] = r->outXor;
                ha->output_phase[1] = r->outMaj;
            }
        } else {
            // Default values when no phase record
            ha->tt_count = 2;
            ha->tt = (unsigned char*)calloc(2, 1);
            ha->input_phase_count = 2;
            ha->input_phase = (int*)calloc(2, sizeof(int));
            ha->output_phase_count = 2;
            ha->output_phase = (int*)calloc(2, sizeof(int));
        }

        id++;
    }

    Vec_WecFree(ha_groupsW);
}

// Populate GiaFullAdder array from extraction results
static void populate_full_adders(
    GiaDetectingResults* results,
    const VecIntResults1* er,
    const ParsedAdderList* LOrder,
    const PhaseDump* phases
) {
    // Count FAs from LOrder
    int fa_count = 0;
    for (int i = 0; i < LOrder->size; i++) {
        if (!LOrder->items[i].isHA) fa_count++;
    }

    if (fa_count <= 0) {
        results->full_adders = NULL;
        results->full_adders_count = 0;
        return;
    }

    results->full_adders = (GiaFullAdder*)calloc(fa_count, sizeof(GiaFullAdder));
    results->full_adders_count = fa_count;

    if (!results->full_adders) {
        results->full_adders_count = 0;
        return;
    }

    // Build FA groups from combined
    Vec_Wec_t* fa_groupsW = groups_from_combined(er->xor3_maj3_combined);

    int id = 0;
    for (int i = 0; i < LOrder->size && id < fa_count; i++) {
        const ParsedAdderItem* A = &LOrder->items[i];
        if (A->isHA) continue;

        GiaFullAdder* fa = &results->full_adders[id];
        fa->id = id;

        // Find matching nodes group
        Vec_Int_t* nodes = NULL;
        for (int gw = 0; gw < Vec_WecSize(fa_groupsW); gw++) {
            Vec_Int_t* G = Vec_WecEntry(fa_groupsW, gw);
            int sz = Vec_IntSize(G);
            if (sz >= 5) {
                int xr = Vec_IntEntry(G, sz - 2);
                int mr = Vec_IntEntry(G, sz - 1);
                if (xr == A->xRoot && mr == A->mRoot) {
                    nodes = G;
                    break;
                }
            }
        }

        // Copy nodes
        if (nodes) {
            fa->nodes = copy_vec_int_to_array(nodes, &fa->nodes_count);
        }

        // Find phase record
        const NpnPhaseRec* r = find_phase_rec(&phases->fa, 0, A->xRoot, A->mRoot);
        if (r) {
            // tt: 2 bytes [tx, tm]
            fa->tt_count = 2;
            fa->tt = (unsigned char*)malloc(2);
            if (fa->tt) {
                fa->tt[0] = r->tx;
                fa->tt[1] = r->tm;
            }

            // input_phase: 3 elements
            fa->input_phase_count = 3;
            fa->input_phase = (int*)malloc(3 * sizeof(int));
            if (fa->input_phase) {
                fa->input_phase[0] = (r->inMask & 1) ? 1 : 0;
                fa->input_phase[1] = (r->inMask & 2) ? 1 : 0;
                fa->input_phase[2] = (r->inMask & 4) ? 1 : 0;
            }

            // output_phase: 2 elements [outXor, outMaj]
            fa->output_phase_count = 2;
            fa->output_phase = (int*)malloc(2 * sizeof(int));
            if (fa->output_phase) {
                fa->output_phase[0] = r->outXor;
                fa->output_phase[1] = r->outMaj;
            }
        } else {
            // Default values when no phase record
            fa->tt_count = 2;
            fa->tt = (unsigned char*)calloc(2, 1);
            fa->input_phase_count = 3;
            fa->input_phase = (int*)calloc(3, sizeof(int));
            fa->output_phase_count = 2;
            fa->output_phase = (int*)calloc(2, sizeof(int));
        }

        id++;
    }

    Vec_WecFree(fa_groupsW);
}

// Populate GiaXorRemain array from extraction results
static void populate_xor_remains(
    GiaDetectingResults* results,
    const VecIntResults1* er,
    const XorPhaseVec* xrPhases
) {
    // Build XOR remaining groups
    Vec_Wec_t* xrW = groups_from_combined(er->xor_remaining);

    // Count valid groups
    int xor_count = 0;
    for (int gi = 0; gi < Vec_WecSize(xrW); gi++) {
        Vec_Int_t* G = Vec_WecEntry(xrW, gi);
        if (Vec_IntSize(G) > 0) xor_count++;
    }

    if (xor_count <= 0) {
        results->xor_remains = NULL;
        results->xor_remains_count = 0;
        Vec_WecFree(xrW);
        return;
    }

    results->xor_remains = (GiaXorRemain*)calloc(xor_count, sizeof(GiaXorRemain));
    results->xor_remains_count = xor_count;

    if (!results->xor_remains) {
        results->xor_remains_count = 0;
        Vec_WecFree(xrW);
        return;
    }

    int id = 0;
    for (int gi = 0; gi < Vec_WecSize(xrW) && id < xor_count; gi++) {
        Vec_Int_t* G = Vec_WecEntry(xrW, gi);
        int sz = Vec_IntSize(G);
        if (sz <= 0) continue;

        GiaXorRemain* xr = &results->xor_remains[id];
        xr->id = id;

        // Copy nodes
        xr->nodes = copy_vec_int_to_array(G, &xr->nodes_count);

        // Find phase record by id
        if (xrPhases && xrPhases->data && id < xrPhases->size) {
            const XorPhaseRec* r = &xrPhases->data[id];

            // tt: 1 byte
            xr->tt_count = 1;
            xr->tt = (unsigned char*)malloc(1);
            if (xr->tt) {
                xr->tt[0] = r->tt;
            }

            // input_phase: 2 or 3 elements depending on nVars
            if (r->nVars == 2) {
                xr->input_phase_count = 2;
                xr->input_phase = (int*)malloc(2 * sizeof(int));
                if (xr->input_phase) {
                    xr->input_phase[0] = (r->inMask & 1) ? 1 : 0;
                    xr->input_phase[1] = (r->inMask & 2) ? 1 : 0;
                }
            } else {
                xr->input_phase_count = 3;
                xr->input_phase = (int*)malloc(3 * sizeof(int));
                if (xr->input_phase) {
                    xr->input_phase[0] = (r->inMask & 1) ? 1 : 0;
                    xr->input_phase[1] = (r->inMask & 2) ? 1 : 0;
                    xr->input_phase[2] = (r->inMask & 4) ? 1 : 0;
                }
            }

            // output_phase: 1 element
            xr->output_phase_count = 1;
            xr->output_phase = (int*)malloc(sizeof(int));
            if (xr->output_phase) {
                xr->output_phase[0] = r->outXor;
            }
        } else {
            // Default values when no phase record
            xr->tt_count = 1;
            xr->tt = (unsigned char*)calloc(1, 1);
            xr->input_phase_count = 2;
            xr->input_phase = (int*)calloc(2, sizeof(int));
            xr->output_phase_count = 1;
            xr->output_phase = (int*)calloc(1, sizeof(int));
        }

        id++;
    }

    Vec_WecFree(xrW);
}

// ======================================================================
// write_detect_results_json: JSON writer (HA/FA + xor_remaining)
// ======================================================================
static int write_detect_results_json(
    const char *path,
    VecIntResults1 *er,
    Gia_Man_t *p,
    Vec_Int_t *vAdds,
    const PhaseDump* phases,
    const ParsedAdderList* LOrder,
    const XorPhaseVec* xrPhases
){
    FILE *f = fopen(path, "w");
    if (!f) {
        perror("fopen failed for detect_faha_output.json");
        return 0;
    }

    // HA/FA counts from filtered combined groups
    int ha_count = vec_group_count_by_zero(er->xor2_maj2_combined);
    int fa_count = vec_group_count_by_zero(er->xor3_maj3_combined);

    // Gather all HA/FA roots from vAdds
    Vec_Int_t *roots = Vec_IntAlloc(2 * (Vec_IntSize(vAdds)/6 + 1));
    int nObjs = Gia_ManObjNum(p);
    for (int i = 0; 6 * i < Vec_IntSize(vAdds); i++) {
        int iXor = Vec_IntEntry(vAdds, 6 * i + 3);
        int iMaj = Vec_IntEntry(vAdds, 6 * i + 4);
        if (iXor > 0 && iXor < nObjs) Vec_IntPush(roots, iXor);
        if (iMaj > 0 && iMaj < nObjs) Vec_IntPush(roots, iMaj);
    }
    Vec_IntSortUnique(roots);

    fprintf(f, "{\n");
    fprintf(f, "  \"ha_count\": %d,\n", ha_count);
    fprintf(f, "  \"fa_count\": %d,\n", fa_count);

    fprintf(f, "  \"roots_all_sorted\": [");
    for (int i = 0; i < Vec_IntSize(roots); i++) {
        if (i) fprintf(f, ", ");
        fprintf(f, "%d", Vec_IntEntry(roots, i));
    }
    fprintf(f, "],\n");

    // Raw grouped vectors (compatibility)
    print_vec_grouped_json(f, "xor2_list", er->xor2_list); fprintf(f, ",\n");
    print_vec_grouped_json(f, "xor3_list", er->xor3_list); fprintf(f, ",\n");
    print_vec_grouped_json(f, "maj2_list", er->maj2_list); fprintf(f, ",\n");
    print_vec_grouped_json(f, "maj3_list", er->maj3_list); fprintf(f, ",\n");
    print_vec_grouped_json(f, "half adders (raw)", er->xor2_maj2_combined); fprintf(f, ",\n");
    print_vec_grouped_json(f, "full adders (raw)", er->xor3_maj3_combined); fprintf(f, ",\n");
    print_vec_grouped_json(f, "xor_all_list", er->xor_all_list); fprintf(f, ",\n");
    print_vec_grouped_json(f, "xor_remaining (excluding_xor used in ha & fa)", er->xor_remaining); fprintf(f, ",\n");

    // Build HA/FA group arrays (nodes per adder) aligned to LOrder
    Vec_Wec_t* fa_groupsW = groups_from_combined(er->xor3_maj3_combined);
    Vec_Wec_t* ha_groupsW = groups_from_combined(er->xor2_maj2_combined);

    // Full adders list: [{id,xor_root,maj_root,nodes:[...]}]
    fprintf(f, "  \"full adders\": [\n");
    if (LOrder) {
        int first = 1;
        int id = 0;
        for (int i=0; i<LOrder->size; ++i) {
            const ParsedAdderItem* A = &LOrder->items[i];
            if (A->isHA) continue;

            Vec_Int_t* nodes = NULL;
            for (int gw=0; gw<Vec_WecSize(fa_groupsW); ++gw) {
                Vec_Int_t* G = Vec_WecEntry(fa_groupsW, gw);
                int sz = Vec_IntSize(G);
                if (sz >= 5) {
                    int xr = Vec_IntEntry(G, sz-2);
                    int mr = Vec_IntEntry(G, sz-1);
                    if (xr == A->xRoot && mr == A->mRoot) {
                        nodes = G;
                        break;
                    }
                }
            }

            if (!first) fprintf(f, ",\n");
            fprintf(f, "    {\"id\": %d, \"xor_root\": %d, \"maj_root\": %d, \"nodes\": [",
                    id, A->xRoot, A->mRoot);
            if (nodes) {
                for (int k=0; k<Vec_IntSize(nodes); ++k) {
                    if (k) fprintf(f, ", ");
                    fprintf(f, "%d", Vec_IntEntry(nodes,k));
                }
            }
            fprintf(f, "]}");
            first = 0;
            id++;
        }
        if (!first) fprintf(f, "\n");
    }
    fprintf(f, "  ],\n");

    // Half adders list: [{id,xor_root,maj_root,nodes:[...]}]
    fprintf(f, "  \"half adders\": [\n");
    if (LOrder) {
        int first = 1;
        int id = 0;
        for (int i=0; i<LOrder->size; ++i) {
            const ParsedAdderItem* A = &LOrder->items[i];
            if (!A->isHA) continue;

            Vec_Int_t* nodes = NULL;
            for (int gw=0; gw<Vec_WecSize(ha_groupsW); ++gw) {
                Vec_Int_t* G = Vec_WecEntry(ha_groupsW, gw);
                int sz = Vec_IntSize(G);
                if (sz >= 4) {
                    int xr = Vec_IntEntry(G, sz-2);
                    int mr = Vec_IntEntry(G, sz-1);
                    if (xr == A->xRoot && mr == A->mRoot) {
                        nodes = G;
                        break;
                    }
                }
            }

            if (!first) fprintf(f, ",\n");
            fprintf(f, "    {\"id\": %d, \"xor_root\": %d, \"maj_root\": %d, \"nodes\": [",
                    id, A->xRoot, A->mRoot);
            if (nodes) {
                for (int k=0; k<Vec_IntSize(nodes); ++k) {
                    if (k) fprintf(f, ", ");
                    fprintf(f, "%d", Vec_IntEntry(nodes,k));
                }
            }
            fprintf(f, "]}");
            first = 0;
            id++;
        }
        if (!first) fprintf(f, "\n");
    }
    fprintf(f, "  ],\n");

    // HA phases: aligned to final HA order
    fprintf(f, "  \"ha_tt_phase\": [\n");
    if (phases && LOrder) {
        int printed = 0;
        int id = 0;
        for (int i=0; i<LOrder->size; ++i) {
            const ParsedAdderItem* A = &LOrder->items[i];
            if (!A->isHA) continue;
            const NpnPhaseRec* r = find_phase_rec(&phases->ha, 1, A->xRoot, A->mRoot);
            if (printed) fprintf(f, ",\n");
            if (r) {
                fprintf(f,
                     "    {\"id\": %d, \"tt\": [\"0x%02X\",\"0x%02X\"], ",
                     id, (int)r->tx, (int)r->tm);
                unsigned a = (r->inMask & 1) ? 1u : 0u;
                unsigned b = (r->inMask & 2) ? 1u : 0u;
                fprintf(f, "\"input_phase\": [%u,%u], ", a, b);
                fprintf(f, "\"input_nodes\": [%d,%d], ",
                        A->in1, A->in2);
                fprintf(f,
                     "\"output_phase\": [%u,%u], \"check\": [%u,%u], "
                     "\"xor_root\": %d, \"maj_root\": %d}",
                     (unsigned)r->outXor, (unsigned)r->outMaj,
                     (unsigned)r->okXor,  (unsigned)r->okMaj,
                     A->xRoot, A->mRoot
                );
            } else {
                fprintf(f,
                    "    {\"id\": %d, \"tt\": [\"0x00\",\"0x00\"], "
                    "\"input_phase\": [0,0], "
                    "\"input_nodes\": [%d,%d], "
                    "\"output_phase\": [0,0], \"check\": [0,0], "
                    "\"xor_root\": %d, \"maj_root\": %d}",
                    id, A->in1, A->in2, A->xRoot, A->mRoot
                );
            }
            printed = 1;
            id++;
        }
        if (printed) fprintf(f, "\n");
    }
    fprintf(f, "  ],\n");

    // FA phases
    fprintf(f, "  \"fa_tt_phase\": [\n");
    if (phases && LOrder) {
        int printed = 0;
        int id = 0;
        for (int i=0; i<LOrder->size; ++i) {
            const ParsedAdderItem* A = &LOrder->items[i];
            if (A->isHA) continue;
            const NpnPhaseRec* r = find_phase_rec(&phases->fa, 0, A->xRoot, A->mRoot);
            if (printed) fprintf(f, ",\n");
            if (r) {
                fprintf(f,
                     "    {\"id\": %d, \"tt\": [\"0x%02X\",\"0x%02X\"], ",
                     id, (int)r->tx, (int)r->tm);
                unsigned a = (r->inMask & 1) ? 1u : 0u;
                unsigned b = (r->inMask & 2) ? 1u : 0u;
                unsigned c = (r->inMask & 4) ? 1u : 0u;
                fprintf(f, "\"input_phase\": [%u,%u,%u], ", a,b,c);
                fprintf(f, "\"input_nodes\": [%d,%d,%d], ",
                        A->in1, A->in2, A->in3);
                fprintf(f,
                     "\"output_phase\": [%u,%u], \"check\": [%u,%u], "
                     "\"xor_root\": %d, \"maj_root\": %d}",
                     (unsigned)r->outXor, (unsigned)r->outMaj,
                     (unsigned)r->okXor,  (unsigned)r->okMaj,
                     A->xRoot, A->mRoot
                );
            } else {
                fprintf(f,
                    "    {\"id\": %d, \"tt\": [\"0x00\",\"0x00\"], "
                    "\"input_phase\": [0,0,0], "
                    "\"input_nodes\": [%d,%d,%d], "
                    "\"output_phase\": [0,0], \"check\": [0,0], "
                    "\"xor_root\": %d, \"maj_root\": %d}",
                    id, A->in1, A->in2, A->in3, A->xRoot, A->mRoot
                );
            }
            printed = 1;
            id++;
        }
        if (printed) fprintf(f, "\n");
    }
    fprintf(f, "  ],\n");

    // xor_remaining_xors: [{id,xor_root,nodes:[...]}]
    fprintf(f, "  \"xor_remaining_xors\": [\n");
    {
        Vec_Wec_t* xrW = groups_from_combined(er->xor_remaining);
        int first = 1;
        int id = 0;
        for (int gi=0; gi<Vec_WecSize(xrW); ++gi) {
            Vec_Int_t* G = Vec_WecEntry(xrW, gi);
            int sz = Vec_IntSize(G);
            if (sz <= 0) continue;
            int xor_root = Vec_IntEntry(G, 0);
            if (!first) fprintf(f, ",\n");
            fprintf(f, "    {\"id\": %d, \"xor_root\": %d, \"nodes\": [", id, xor_root);
            for (int k=0; k<sz; ++k) {
                if (k) fprintf(f, ", ");
                fprintf(f, "%d", Vec_IntEntry(G,k));
            }
            fprintf(f, "]}");
            first = 0;
            id++;
        }
        if (!first) fprintf(f, "\n");
        Vec_WecFree(xrW);
    }
    fprintf(f, "  ],\n");

    // xor_remaining_phase: XOR2 or XOR3 depending on nVars
    fprintf(f, "  \"xor_remaining_phase\": [\n");
    if (xrPhases && xrPhases->data) {
        for (int id=0; id<xrPhases->size; ++id) {
            const XorPhaseRec* r = &xrPhases->data[id];
            if (id) fprintf(f, ",\n");
            fprintf(f, "    {\"id\": %d, ", id);
            fprintf(f, "\"n_inputs\": %d, ", r->nVars);
            fprintf(f, "\"tt\": \"0x%02X\", ", (int)r->tt);

            if (r->nVars == 2) {
                unsigned a = (r->inMask & 1) ? 1u : 0u;
                unsigned b = (r->inMask & 2) ? 1u : 0u;
                fprintf(f, "\"input_phase\": [%u,%u], ", a,b);
                fprintf(f, "\"input_nodes\": [%d,%d], ",
                        r->in1, r->in2);
            } else {
                unsigned a = (r->inMask & 1) ? 1u : 0u;
                unsigned b = (r->inMask & 2) ? 1u : 0u;
                unsigned c = (r->inMask & 4) ? 1u : 0u;
                fprintf(f, "\"input_phase\": [%u,%u,%u], ", a,b,c);
                fprintf(f, "\"input_nodes\": [%d,%d,%d], ",
                        r->in1, r->in2, r->in3);
            }

            fprintf(f, "\"group_last_two_nodes\": [%d,%d], ",
                    r->gLast1, r->gLast2);

            fprintf(f, "\"output_phase\": [%u], ", (unsigned)r->outXor);
            fprintf(f, "\"check\": [%u], ", (unsigned)r->ok);
            fprintf(f, "\"xor_root\": %d}", r->xorRoot);
        }
        if (xrPhases->size)
            fprintf(f, "\n");
    }
    fprintf(f, "  ]\n");

    fprintf(f, "}\n");

    Vec_IntFree(roots);
    Vec_WecFree(fa_groupsW);
    Vec_WecFree(ha_groupsW);

    if (fclose(f) != 0) {
        perror("fclose failed for detect_faha_output.json");
        return 0;
    }
    return 1;
}

// ======================================================================
// Main HA/FA detection driver
// ======================================================================

extern Vec_Int_t * Ree_ManComputeCuts( Gia_Man_t * p, Vec_Int_t ** pvXors, int fVerbose );

// Main function: detect HA/FA and XOR-remaining, compute phases, and emit JSON
GiaDetectingResults* Gia_detecting_ha_fa(Gia_Man_t *p)
{
    if (!p) {
        printf("Input Gia_Man_t is NULL\n");
        return NULL;
    }

    GiaDetectingResults *results = (GiaDetectingResults *)malloc(sizeof(GiaDetectingResults));
    if (!results) {
        printf("malloc GiaDetectingResults failed\n");
        return NULL;
    }
    memset(results, 0, sizeof(GiaDetectingResults));

    Vec_Int_t *vXors = NULL;
    Vec_Int_t *vAdds = Ree_ManComputeCuts(p, &vXors, 0);
    Vec_Int_t *vAdds_copy = Vec_IntDup(vAdds);
    Vec_Int_t *vXors_copy = Vec_IntDup(vXors);

    // 1) Extract HA/FA candidates and XOR lists
    VecIntResults1 er = extract_ha_fa_labeling(p, vAdds_copy, vXors_copy);

    // 2) Build parsed adder list from unfiltered combined groups
    ParsedAdderList L;
    BuildAdderListFromExtracted(er.xor2_maj2_combined, er.xor3_maj3_combined, &L);

    // 3) Truth-table engine + NPN check for HA/FA
    Gia_ObjComputeTruthTableStart(p, 16);
    NpnPhaseVec phasesAll; NpvInit(&phasesAll);
#ifdef ABC_USE_PTHREADS
    // Use multi-threaded version with up to 16 threads (NPN matching is parallelized)
    int* passFlags = Gia_CheckNpnAdder_FromExtracted_MT(p, &L, &phasesAll, 16);
#else
    int* passFlags = Gia_CheckNpnAdder_FromExtracted(p, &L, &phasesAll);
#endif

    // 3.5) XOR remaining: build XOR list (2-input or 3-input) + phase records
    ParsedXorList XL;
    BuildXorListFromRemaining(p, er.xor_remaining, &XL);
    XorPhaseVec xrPhases; XpvInit(&xrPhases);
    CheckXorRemainingPhases(p, &XL, &xrPhases);

    Gia_ObjComputeTruthTableStop(p);

    // 4) Build pass map for filtering HA/FA groups
    HmPass* map = NULL;
    int mapSize = 0;
    if (passFlags) {
        mapSize = L.size;
        map = (HmPass*)malloc(sizeof(HmPass)*mapSize);
        for (int i=0; i<L.size; ++i) {
            map[i].isHA  = L.items[i].isHA;
            map[i].xRoot = L.items[i].xRoot;
            map[i].mRoot = L.items[i].mRoot;
            map[i].pass  = passFlags[i] ? 1 : 0;
        }
    }

    // 5) Filter combined groups by passFlags
    Vec_Int_t* xor2_maj2_filtered =
        passFlags ? FilterCombinedByPass(er.xor2_maj2_combined, map, mapSize, 1)
                  : Vec_IntDup(er.xor2_maj2_combined);
    Vec_Int_t* xor3_maj3_filtered =
        passFlags ? FilterCombinedByPass(er.xor3_maj3_combined, map, mapSize, 0)
                  : Vec_IntDup(er.xor3_maj3_combined);

    if (passFlags) {
        Vec_IntFree(er.xor2_maj2_combined);
        Vec_IntFree(er.xor3_maj3_combined);
        er.xor2_maj2_combined = xor2_maj2_filtered;
        er.xor3_maj3_combined = xor3_maj3_filtered;
    } else {
        Vec_IntFree(xor2_maj2_filtered);
        Vec_IntFree(xor3_maj3_filtered);
    }

    // 6) Prepare result vectors
    results->xor2_list          = Vec_IntDup(er.xor2_list);
    results->xor3_list          = Vec_IntDup(er.xor3_list);
    results->maj2_list          = Vec_IntDup(er.maj2_list);
    results->maj3_list          = Vec_IntDup(er.maj3_list);
    results->xor_all_list       = Vec_IntDup(er.xor_all_list);
    results->xor_remaining      = Vec_IntDup(er.xor_remaining);
    results->xor2_maj2_combined = Vec_IntDup(er.xor2_maj2_combined);
    results->xor3_maj3_combined = Vec_IntDup(er.xor3_maj3_combined);

    // 7) Roots from filtered HA/FA groups
    results->ha_roots  = Vec_IntAlloc(128);
    results->fa_roots  = Vec_IntAlloc(128);
    results->all_roots = Vec_IntAlloc(256);

    extract_group_heads(results->xor2_maj2_combined, results->ha_roots);
    extract_group_heads(results->xor3_maj3_combined, results->fa_roots);

    {
        int x,i;
        Vec_IntForEachEntry(results->ha_roots, x, i)
            Vec_IntPush(results->all_roots, x);
        Vec_IntForEachEntry(results->fa_roots, x, i)
            Vec_IntPush(results->all_roots, x);
    }
    {
        int w=0, n=Vec_IntSize(results->all_roots);
        for (int i=0; i<n; ++i) {
            int x = Vec_IntEntry(results->all_roots,i);
            if (x > 0) Vec_IntWriteEntry(results->all_roots, w++, x);
        }
        Vec_IntShrink(results->all_roots, w);
        Vec_IntSortUnique(results->all_roots);
        UniqueIdsInPlace(p, results->all_roots);
        UniqueIdsInPlace(p, results->ha_roots);
        UniqueIdsInPlace(p, results->fa_roots);
    }

    // 7.5) Rebuild ParsedAdderList from FILTERED groups (for aligned output)
    ParsedAdderList L_filtered;
    BuildAdderListFromExtracted(er.xor2_maj2_combined, er.xor3_maj3_combined, &L_filtered);

    // 8) Split HA/FA phase records into PhaseDump
    PhaseDump phaseDump; PhaseDumpInit(&phaseDump);
    for (int i=0; i<phasesAll.size; ++i) {
        NpnPhaseRec r = phasesAll.data[i];
        if (r.isHA) NpvPush(&phaseDump.ha, &r);
        else        NpvPush(&phaseDump.fa, &r);
    }

    // 8.5) Populate structured HA/FA/XOR arrays for direct C++/Talisman consumption
    {
        VecIntResults1 er_for_populate = {
            .xor2_list          = results->xor2_list,
            .xor3_list          = results->xor3_list,
            .maj2_list          = results->maj2_list,
            .maj3_list          = results->maj3_list,
            .xor_all_list       = results->xor_all_list,
            .xor_remaining      = results->xor_remaining,
            .xor2_maj2_combined = results->xor2_maj2_combined,
            .xor3_maj3_combined = results->xor3_maj3_combined
        };
        populate_half_adders(results, &er_for_populate, &L_filtered, &phaseDump);
        populate_full_adders(results, &er_for_populate, &L_filtered, &phaseDump);
        populate_xor_remains(results, &er_for_populate, &xrPhases);
    }

    if (g_acecVerbose >= 2)
        printf("Detection complete: half_adders=%d, full_adders=%d, xor_remains=%d\n",
               results->half_adders_count,
               results->full_adders_count,
               results->xor_remains_count);

    // 9) Cleanup
    if (passFlags) free(passFlags);
    if (map)       free(map);
    PAL_Free(&L);
    PAL_Free(&L_filtered);
    NpvFree(&phasesAll);
    PhaseDumpFree(&phaseDump);
    PXL_Free(&XL);
    XpvFree(&xrPhases);

    Vec_IntFree(vAdds);
    Vec_IntFree(vXors);
    Vec_IntFree(vAdds_copy);
    Vec_IntFree(vXors_copy);
    Vec_IntFree(er.xor2_list);
    Vec_IntFree(er.xor3_list);
    Vec_IntFree(er.maj2_list);
    Vec_IntFree(er.maj3_list);
    Vec_IntFree(er.xor_all_list);
    Vec_IntFree(er.xor_remaining);
    Vec_IntFree(er.xor2_maj2_combined);
    Vec_IntFree(er.xor3_maj3_combined);

    return results;
}

// Free results
void Gia_free_detecting_results(GiaDetectingResults *results) {
    if (!results) return;
    if (results->xor2_list)           Vec_IntFree(results->xor2_list);
    if (results->xor3_list)           Vec_IntFree(results->xor3_list);
    if (results->maj2_list)           Vec_IntFree(results->maj2_list);
    if (results->maj3_list)           Vec_IntFree(results->maj3_list);
    if (results->xor_all_list)        Vec_IntFree(results->xor_all_list);
    if (results->xor_remaining)       Vec_IntFree(results->xor_remaining);
    if (results->xor2_maj2_combined)  Vec_IntFree(results->xor2_maj2_combined);
    if (results->xor3_maj3_combined)  Vec_IntFree(results->xor3_maj3_combined);
    if (results->ha_roots)            Vec_IntFree(results->ha_roots);
    if (results->fa_roots)            Vec_IntFree(results->fa_roots);
    if (results->all_roots)           Vec_IntFree(results->all_roots);

    // Free half adders array
    if (results->half_adders) {
        for (int i = 0; i < results->half_adders_count; i++) {
            GiaHalfAdder* ha = &results->half_adders[i];
            if (ha->nodes)        free(ha->nodes);
            if (ha->tt)           free(ha->tt);
            if (ha->input_phase)  free(ha->input_phase);
            if (ha->output_phase) free(ha->output_phase);
        }
        free(results->half_adders);
    }

    // Free full adders array
    if (results->full_adders) {
        for (int i = 0; i < results->full_adders_count; i++) {
            GiaFullAdder* fa = &results->full_adders[i];
            if (fa->nodes)        free(fa->nodes);
            if (fa->tt)           free(fa->tt);
            if (fa->input_phase)  free(fa->input_phase);
            if (fa->output_phase) free(fa->output_phase);
        }
        free(results->full_adders);
    }

    // Free xor remains array
    if (results->xor_remains) {
        for (int i = 0; i < results->xor_remains_count; i++) {
            GiaXorRemain* xr = &results->xor_remains[i];
            if (xr->nodes)        free(xr->nodes);
            if (xr->tt)           free(xr->tt);
            if (xr->input_phase)  free(xr->input_phase);
            if (xr->output_phase) free(xr->output_phase);
        }
        free(results->xor_remains);
    }

    free(results);
}

// Extended version of Gia_detecting_ha_fa that also returns vAdds and vXors for reuse
// This avoids redundant Ree_ManComputeCuts calls in Gia_PoLocalKCover_reverse
GiaDetectingResults* Gia_detecting_ha_fa_ext(Gia_Man_t *p, Vec_Int_t ** ppvAdds, Vec_Int_t ** ppvXors)
{
    if (!p) {
        printf("Input Gia_Man_t is NULL\n");
        if (ppvAdds) *ppvAdds = NULL;
        if (ppvXors) *ppvXors = NULL;
        return NULL;
    }

    GiaDetectingResults *results = (GiaDetectingResults *)malloc(sizeof(GiaDetectingResults));
    if (!results) {
        printf("malloc GiaDetectingResults failed\n");
        if (ppvAdds) *ppvAdds = NULL;
        if (ppvXors) *ppvXors = NULL;
        return NULL;
    }
    memset(results, 0, sizeof(GiaDetectingResults));

    Vec_Int_t *vXors = NULL;
    Vec_Int_t *vAdds = Ree_ManComputeCuts(p, &vXors, 0);
    Vec_Int_t *vAdds_copy = Vec_IntDup(vAdds);
    Vec_Int_t *vXors_copy = Vec_IntDup(vXors);

    // 1) Extract HA/FA candidates and XOR lists
    VecIntResults1 er = extract_ha_fa_labeling(p, vAdds_copy, vXors_copy);

    // 2) Build parsed adder list from unfiltered combined groups
    ParsedAdderList L;
    BuildAdderListFromExtracted(er.xor2_maj2_combined, er.xor3_maj3_combined, &L);

    // 3) Truth-table engine + NPN check for HA/FA
    Gia_ObjComputeTruthTableStart(p, 16);
    NpnPhaseVec phasesAll; NpvInit(&phasesAll);
#ifdef ABC_USE_PTHREADS
    // Use multi-threaded version with up to 16 threads (NPN matching is parallelized)
    int* passFlags = Gia_CheckNpnAdder_FromExtracted_MT(p, &L, &phasesAll, 16);
#else
    int* passFlags = Gia_CheckNpnAdder_FromExtracted(p, &L, &phasesAll);
#endif

    // 3.5) XOR remaining: build XOR list (2-input or 3-input) + phase records
    ParsedXorList XL;
    BuildXorListFromRemaining(p, er.xor_remaining, &XL);
    XorPhaseVec xrPhases; XpvInit(&xrPhases);
    CheckXorRemainingPhases(p, &XL, &xrPhases);

    Gia_ObjComputeTruthTableStop(p);

    // 4) Build pass map for filtering HA/FA groups
    HmPass* map = NULL;
    int mapSize = 0;
    if (passFlags) {
        mapSize = L.size;
        map = (HmPass*)malloc(sizeof(HmPass)*mapSize);
        for (int i=0; i<L.size; ++i) {
            map[i].isHA  = L.items[i].isHA;
            map[i].xRoot = L.items[i].xRoot;
            map[i].mRoot = L.items[i].mRoot;
            map[i].pass  = passFlags[i] ? 1 : 0;
        }
    }

    // 5) Filter combined groups by passFlags
    Vec_Int_t* xor2_maj2_filtered =
        passFlags ? FilterCombinedByPass(er.xor2_maj2_combined, map, mapSize, 1)
                  : Vec_IntDup(er.xor2_maj2_combined);
    Vec_Int_t* xor3_maj3_filtered =
        passFlags ? FilterCombinedByPass(er.xor3_maj3_combined, map, mapSize, 0)
                  : Vec_IntDup(er.xor3_maj3_combined);

    if (passFlags) {
        Vec_IntFree(er.xor2_maj2_combined);
        Vec_IntFree(er.xor3_maj3_combined);
        er.xor2_maj2_combined = xor2_maj2_filtered;
        er.xor3_maj3_combined = xor3_maj3_filtered;
    } else {
        Vec_IntFree(xor2_maj2_filtered);
        Vec_IntFree(xor3_maj3_filtered);
    }

    // 6) Prepare result vectors
    results->xor2_list          = Vec_IntDup(er.xor2_list);
    results->xor3_list          = Vec_IntDup(er.xor3_list);
    results->maj2_list          = Vec_IntDup(er.maj2_list);
    results->maj3_list          = Vec_IntDup(er.maj3_list);
    results->xor_all_list       = Vec_IntDup(er.xor_all_list);
    results->xor_remaining      = Vec_IntDup(er.xor_remaining);
    results->xor2_maj2_combined = Vec_IntDup(er.xor2_maj2_combined);
    results->xor3_maj3_combined = Vec_IntDup(er.xor3_maj3_combined);

    // 7) Roots from filtered HA/FA groups
    results->ha_roots  = Vec_IntAlloc(128);
    results->fa_roots  = Vec_IntAlloc(128);
    results->all_roots = Vec_IntAlloc(256);

    extract_group_heads(results->xor2_maj2_combined, results->ha_roots);
    extract_group_heads(results->xor3_maj3_combined, results->fa_roots);

    {
        int x,i;
        Vec_IntForEachEntry(results->ha_roots, x, i)
            Vec_IntPush(results->all_roots, x);
        Vec_IntForEachEntry(results->fa_roots, x, i)
            Vec_IntPush(results->all_roots, x);
    }
    {
        int w=0, n=Vec_IntSize(results->all_roots);
        for (int i=0; i<n; ++i) {
            int x = Vec_IntEntry(results->all_roots,i);
            if (x > 0) Vec_IntWriteEntry(results->all_roots, w++, x);
        }
        Vec_IntShrink(results->all_roots, w);
        Vec_IntSortUnique(results->all_roots);
        UniqueIdsInPlace(p, results->all_roots);
        UniqueIdsInPlace(p, results->ha_roots);
        UniqueIdsInPlace(p, results->fa_roots);
    }

    // 7.5) Rebuild ParsedAdderList from FILTERED groups (for aligned output)
    ParsedAdderList L_filtered;
    BuildAdderListFromExtracted(er.xor2_maj2_combined, er.xor3_maj3_combined, &L_filtered);

    // 8) Split HA/FA phase records into PhaseDump
    PhaseDump phaseDump; PhaseDumpInit(&phaseDump);
    for (int i=0; i<phasesAll.size; ++i) {
        NpnPhaseRec r = phasesAll.data[i];
        if (r.isHA) NpvPush(&phaseDump.ha, &r);
        else        NpvPush(&phaseDump.fa, &r);
    }

    // 8.5) Populate structured HA/FA/XOR arrays for direct C++/Talisman consumption
    {
        VecIntResults1 er_for_populate = {
            .xor2_list          = results->xor2_list,
            .xor3_list          = results->xor3_list,
            .maj2_list          = results->maj2_list,
            .maj3_list          = results->maj3_list,
            .xor_all_list       = results->xor_all_list,
            .xor_remaining      = results->xor_remaining,
            .xor2_maj2_combined = results->xor2_maj2_combined,
            .xor3_maj3_combined = results->xor3_maj3_combined
        };
        populate_half_adders(results, &er_for_populate, &L_filtered, &phaseDump);
        populate_full_adders(results, &er_for_populate, &L_filtered, &phaseDump);
        populate_xor_remains(results, &er_for_populate, &xrPhases);
    }

    if (g_acecVerbose >= 2)
        printf("Detection complete (ext): half_adders=%d, full_adders=%d, xor_remains=%d\n",
               results->half_adders_count,
               results->full_adders_count,
               results->xor_remains_count);

    // 9) Cleanup internal structures
    if (passFlags) free(passFlags);
    if (map)       free(map);
    PAL_Free(&L);
    PAL_Free(&L_filtered);
    NpvFree(&phasesAll);
    PhaseDumpFree(&phaseDump);
    PXL_Free(&XL);
    XpvFree(&xrPhases);

    // Return vAdds and vXors to caller for reuse (caller is responsible for freeing)
    if (ppvAdds) {
        *ppvAdds = vAdds;
    } else {
        Vec_IntFree(vAdds);
    }
    if (ppvXors) {
        *ppvXors = vXors;
    } else {
        Vec_IntFree(vXors);
    }

    Vec_IntFree(vAdds_copy);
    Vec_IntFree(vXors_copy);
    Vec_IntFree(er.xor2_list);
    Vec_IntFree(er.xor3_list);
    Vec_IntFree(er.maj2_list);
    Vec_IntFree(er.maj3_list);
    Vec_IntFree(er.xor_all_list);
    Vec_IntFree(er.xor_remaining);
    Vec_IntFree(er.xor2_maj2_combined);
    Vec_IntFree(er.xor3_maj3_combined);

    return results;
}


// ======================== Reverse-cut local K-cover (cleaned, keep U2 only + JSON) ========================

// --------------------- Vec helpers ---------------------
static Vec_Int_t * ExportBitsetIds_NoZero( Vec_Bit_t * vBits )
{
    Vec_Int_t * v = Vec_IntAlloc( 1024 );
    int n = Vec_BitSize( vBits );
    for ( int i = 0; i < n; ++i )
        if ( i != 0 && Vec_BitEntry( vBits, i ) )
            Vec_IntPush( v, i );
    Vec_IntSort( v, 0 );
    Vec_IntUniqify( v );
    return v;
}

static Vec_Int_t * ExportBitsetIds_All( Vec_Bit_t * vBits )
{
    Vec_Int_t * v = Vec_IntAlloc( 1024 );
    int n = Vec_BitSize( vBits );
    for ( int i = 0; i < n; ++i )
        if ( Vec_BitEntry( vBits, i ) )
            Vec_IntPush( v, i );
    Vec_IntSort( v, 0 );
    Vec_IntUniqify( v );
    return v;
}

// --------------------- Boundary helpers ---------------------
static inline int Gia_ObjIsBoundaryStop( Gia_Man_t * p, Vec_Bit_t * vBoundaryBits, int id )
{
    int nObjs = Gia_ManObjNum(p);
    if ( id <= 0 || id >= nObjs ) return 1;
    if ( vBoundaryBits && Vec_BitEntry( vBoundaryBits, id ) ) return 1;
    Gia_Obj_t * pObj = Gia_ManObj( p, id );
    if ( pObj == NULL ) return 1;
    if ( Gia_ObjIsConst0( pObj ) ) return 1;
    if ( Gia_ObjIsPi( p, pObj ) ) return 1;
    return 0;
}

// --------------------- HA/FA mapping ---------------------
typedef struct {
    Vec_Bit_t *isRoot;
    Vec_Int_t *root2pair;
    Vec_Wec_t *pairInputs;
    Vec_Wec_t *pairRoots;
    int pairCnt;
} RMaps2;

static RMaps2 * BuildMapsFromAdds2( Gia_Man_t * p, Vec_Int_t * vAdds )
{
    RMaps2 *M = (RMaps2 *)ABC_CALLOC( RMaps2, 1 );
    int nObjs = Gia_ManObjNum( p );
    int nPairs = Vec_IntSize( vAdds ) / 6;
    M->isRoot     = Vec_BitStart( nObjs );
    M->root2pair  = Vec_IntStartFull( nObjs );
    M->pairInputs = Vec_WecStart( nPairs );
    M->pairRoots  = Vec_WecStart( nPairs );
    M->pairCnt    = nPairs;

    for (int i = 0; 6*i < Vec_IntSize(vAdds); i++ )
    {
        int in0 = Vec_IntEntry( vAdds, 6*i+0 );
        int in1 = Vec_IntEntry( vAdds, 6*i+1 );
        int in2 = Vec_IntEntry( vAdds, 6*i+2 );
        int iX  = Vec_IntEntry( vAdds, 6*i+3 );
        int iM  = Vec_IntEntry( vAdds, 6*i+4 );

        Vec_Int_t *vIn = Vec_WecEntry( M->pairInputs, i );
        if (in0 > 0 && in0 < nObjs) Vec_IntPush( vIn, in0 );
        if (in1 > 0 && in1 < nObjs) Vec_IntPush( vIn, in1 );
        if (in2 > 0 && in2 < nObjs) Vec_IntPush( vIn, in2 );

        Vec_Int_t *vRt = Vec_WecEntry( M->pairRoots, i );
        if (iX > 0 && iX < nObjs) { Vec_IntPush(vRt, iX); Vec_BitWriteEntry(M->isRoot, iX, 1); Vec_IntWriteEntry(M->root2pair, iX, i); }
        if (iM > 0 && iM < nObjs) { Vec_IntPush(vRt, iM); Vec_BitWriteEntry(M->isRoot, iM, 1); Vec_IntWriteEntry(M->root2pair, iM, i); }
    }
    return M;
}

static void FreeMaps2( RMaps2 *M )
{
    if (!M) return;
    if (M->isRoot)     Vec_BitFree(M->isRoot);
    if (M->root2pair)  Vec_IntFree(M->root2pair);
    if (M->pairInputs) Vec_WecFree(M->pairInputs);
    if (M->pairRoots)  Vec_WecFree(M->pairRoots);
    ABC_FREE(M);
}

// --------------------- Pair cone collection ---------------------
static int CollectPairConeNodes( Gia_Man_t * p, Vec_Int_t * vPairRoots, Vec_Bit_t * vBoundaryInputs, Vec_Bit_t * vHitBits )
{
    int nObjs = Gia_ManObjNum(p);
    Vec_Int_t * st = Vec_IntAlloc(64);
    int i, r, added = 0;
    Vec_IntForEachEntry( vPairRoots, r, i )
        if ( r > 0 && r < nObjs ) Vec_IntPush( st, r );

    while ( Vec_IntSize( st ) )
    {
        int u = Vec_IntPop( st );
        if ( u <= 0 || u >= nObjs ) continue;
        if ( Vec_BitEntry( vHitBits, u ) ) continue;
        Vec_BitWriteEntry( vHitBits, u, 1 );
        added++;

        if ( Gia_ObjIsBoundaryStop( p, vBoundaryInputs, u ) )
            continue;

        Gia_Obj_t * o = Gia_ManObj( p, u );
        if ( !o ) continue;

        if ( Gia_ObjIsAnd( o ) )
        {
            int f0 = Gia_ObjFaninId0( o, u );
            int f1 = Gia_ObjFaninId1( o, u );
            if ( f0 > 0 && f0 < nObjs ) Vec_IntPush( st, f0 );
            if ( f1 > 0 && f1 < nObjs ) Vec_IntPush( st, f1 );
        }
        else if ( Gia_ObjIsCo( o ) )
        {
            int f0 = Gia_ObjFaninId0( o, u );
            if ( f0 > 0 && f0 < nObjs ) Vec_IntPush( st, f0 );
        }
    }
    Vec_IntFree( st );
    return added;
}

// --------------------- Prefer-new pair finder ---------------------
static void FindPoFirstHitPair2_PrefNew(
    Gia_Man_t * p,
    int poId,
    RMaps2 * M,
    Vec_Bit_t * isPairHit,
    Vec_Bit_t * bannedHitNodes,
    int * pRootId,
    int * pPairId
)
{
    *pRootId = -1; *pPairId = -1;
    int nObjs = Gia_ManObjNum(p);
    if (poId <= 0 || poId >= nObjs) return;

    Vec_Int_t * q   = Vec_IntAlloc(256);
    Vec_Bit_t * vis = Vec_BitStart(nObjs);

    Vec_IntPush(q, poId);
    Vec_BitWriteEntry(vis, poId, 1);

    int q_head = 0;
    while (q_head < Vec_IntSize(q)) {
        int level_begin = q_head;
        int level_end   = Vec_IntSize(q);

        int bestNewRoot = -1, bestNewPair = -1;
        int bestOldRoot = -1, bestOldPair = -1;

        for (int idx = level_begin; idx < level_end; ++idx) {
            int u = Vec_IntEntry(q, idx);
            if ( Vec_BitEntry(bannedHitNodes, u) )
                continue;
            if ( Vec_BitEntry(M->isRoot, u) ) {
                int pid = Vec_IntEntry(M->root2pair, u);
                if (pid >= 0 && pid < M->pairCnt) {
                    if (!Vec_BitEntry(isPairHit, pid)) {
                        if (bestNewPair < 0) { bestNewRoot = u; bestNewPair = pid; }
                    } else {
                        if (bestOldPair < 0) { bestOldRoot = u; bestOldPair = pid; }
                    }
                }
            }
        }
        if (bestNewPair >= 0) { *pRootId = bestNewRoot; *pPairId = bestNewPair; Vec_IntFree(q); Vec_BitFree(vis); return; }
        if (bestOldPair >= 0) { *pRootId = bestOldRoot; *pPairId = bestOldPair; Vec_IntFree(q); Vec_BitFree(vis); return; }

        for (int idx = level_begin; idx < level_end; ++idx) {
            int u = Vec_IntEntry(q, idx);
            Gia_Obj_t * o = Gia_ManObj(p, u);
            if (!o) continue;

            if (Gia_ObjIsCo(o)) {
                int a = Gia_ObjFaninId0(o, u);
                if (a > 0 && a < nObjs && !Vec_BitEntry(vis, a)) {
                    Vec_IntPush(q, a); Vec_BitWriteEntry(vis, a, 1);
                }
            } else if (Gia_ObjIsAnd(o)) {
                int a = Gia_ObjFaninId0(o, u);
                int b = Gia_ObjFaninId1(o, u);
                if (a > 0 && a < nObjs && !Vec_BitEntry(vis, a)) {
                    Vec_IntPush(q, a); Vec_BitWriteEntry(vis, a, 1);
                }
                if (b > 0 && b < nObjs && !Vec_BitEntry(vis, b)) {
                    Vec_IntPush(q, b); Vec_BitWriteEntry(vis, b, 1);
                }
            }
        }
        q_head = level_end;
    }
    Vec_IntFree(q);
    Vec_BitFree(vis);
}

// --------------------- Per-PO U2 only ---------------------
static Vec_Bit_t * BuildPoBorderBitset( int nObjs, Vec_Int_t * vInputs )
{
    Vec_Bit_t * v = Vec_BitStart( nObjs );
    int i, x; Vec_IntForEachEntry( vInputs, x, i )
        if ( x > 0 && x < nObjs )
            Vec_BitWriteEntry( v, x, 1 );
    return v;
}

static void ComputeU2ForPo_POInputs_UsingGiaLevel_GlobalBorder(
    Gia_Man_t * p,
    int poId,
    Vec_Bit_t * vPoBorders,
    Vec_Bit_t * vBorderGlobal,
    int minInputLevel,
    int hasInputLevel,
    Vec_Bit_t * vMiddleLocal,
    Vec_Bit_t * vNewBorders,
    Vec_Bit_t * vLevelBorders
)
{
    const int nObjs = Gia_ManObjNum(p);

    Vec_Int_t * st = Vec_IntAlloc(64);
    Vec_Bit_t * vis = Vec_BitStart(nObjs);

    Vec_IntPush(st, poId);

    while ( Vec_IntSize(st) )
    {
        int u = Vec_IntPop(st);
        if ( u <= 0 || u >= nObjs ) continue;
        if ( Vec_BitEntry(vis, u) ) continue;
        Vec_BitWriteEntry(vis, u, 1);

        Gia_Obj_t * o = Gia_ManObj(p, u);
        if ( !o ) continue;

        int isPoBorder   = (vPoBorders     && Vec_BitEntry(vPoBorders,   u));
        int isGlobBorder = (vBorderGlobal  && Vec_BitEntry(vBorderGlobal,u));
        int isPIorConst  = Gia_ObjIsPi(p, o) || Gia_ObjIsConst0(o);

        if ( isPoBorder || isGlobBorder || isPIorConst )
        {
            if ( !isGlobBorder )
                Vec_BitWriteEntry(vNewBorders, u, 1);
            continue;
        }

        Vec_BitWriteEntry(vMiddleLocal, u, 1);

        if ( hasInputLevel )
        {
            int lv = Gia_ObjLevel(p, o);
            if ( lv < minInputLevel )
            {
                Vec_BitWriteEntry(vNewBorders,   u, 1);
                Vec_BitWriteEntry(vLevelBorders, u, 1);
                continue;
            }
        }

        if ( Gia_ObjIsAnd(o) )
        {
            int f0 = Gia_ObjFaninId0(o, u);
            int f1 = Gia_ObjFaninId1(o, u);
            if ( f0 > 0 && f0 < nObjs ) Vec_IntPush(st, f0);
            if ( f1 > 0 && f1 < nObjs ) Vec_IntPush(st, f1);
        }
        else if ( Gia_ObjIsCo( o ) )
        {
            int f0 = Gia_ObjFaninId0( o, u );
            if ( f0 > 0 && f0 < nObjs ) Vec_IntPush(st, f0);
        }
    }

    Vec_BitFree(vis);
    Vec_IntFree(st);
}

static int FixClosure_AddMissingFaninToBorder( Gia_Man_t * p, Vec_Bit_t * vBorderAll, Vec_Bit_t * vMiddleAll )
{
    int nObjs = Gia_ManObjNum( p );
    int totalAdded = 0;
    int changed;

    do {
        changed = 0;

        Gia_Obj_t * o; int id;
        Gia_ManForEachObj( p, o, id )
        {
            int f0, f1;

            if ( id <= 0 || id >= nObjs )
                continue;

            if ( !Vec_BitEntry( vMiddleAll, id ) )
                continue;

            if ( Gia_ObjIsPi( p, o ) || Gia_ObjIsConst0( o ) )
                continue;

            if ( Gia_ObjIsAnd( o ) )
            {
                f0 = Gia_ObjFaninId0( o, id );
                f1 = Gia_ObjFaninId1( o, id );

                if ( f0 > 0 && f0 < nObjs )
                {
                    int inMiddle = Vec_BitEntry( vMiddleAll, f0 );
                    int inBorder = Vec_BitEntry( vBorderAll, f0 );
                    if ( !inMiddle && !inBorder )
                    {
                        Vec_BitWriteEntry( vBorderAll, f0, 1 );
                        changed++; totalAdded++;
                    }
                }
                if ( f1 > 0 && f1 < nObjs )
                {
                    int inMiddle = Vec_BitEntry( vMiddleAll, f1 );
                    int inBorder = Vec_BitEntry( vBorderAll, f1 );
                    if ( !inMiddle && !inBorder )
                    {
                        Vec_BitWriteEntry( vBorderAll, f1, 1 );
                        changed++; totalAdded++;
                    }
                }
            }
            else if ( Gia_ObjIsCo( o ) )
            {
                f0 = Gia_ObjFaninId0( o, id );
                if ( f0 > 0 && f0 < nObjs )
                {
                    int inMiddle = Vec_BitEntry( vMiddleAll, f0 );
                    int inBorder = Vec_BitEntry( vBorderAll, f0 );
                    if ( !inMiddle && !inBorder )
                    {
                        Vec_BitWriteEntry( vBorderAll, f0, 1 );
                        changed++; totalAdded++;
                    }
                }
            }
        }

        if ( changed )
        {
            int j;
            for ( j = 1; j < nObjs; ++j )
            {
                if ( Vec_BitEntry( vBorderAll, j ) && Vec_BitEntry( vMiddleAll, j ) )
                    Vec_BitWriteEntry( vMiddleAll, j, 0 );
            }
        }
    } while ( changed > 0 );

    return totalAdded;
}

static void RunClosureByTraversalU2_POInputs( Gia_Man_t * p, Vec_Wec_t * poInputs, Vec_Bit_t ** ppBorderAll, Vec_Bit_t ** ppMiddleAll )
{
    const int nObjs = Gia_ManObjNum(p);

    Gia_ManLevelNum(p);

    Vec_Bit_t * vMiddleAll = Vec_BitStart( nObjs );
    Vec_Bit_t * vBorderAll = Vec_BitStart( nObjs );

    Gia_Obj_t * pPo; int i;
    Gia_ManForEachPo( p, pPo, i )
    {
        int poId = Gia_ObjId( p, pPo );

        Vec_Int_t * vIn = Vec_WecEntry( poInputs, i );
        Vec_Bit_t * vPoB = BuildPoBorderBitset( nObjs, vIn );

        int hasMin = 0;
        int minLv  = 0;
        {
            int has = 0;
            int mn  = ABC_INFINITY;
            int k, x;
            Vec_IntForEachEntry( vIn, x, k )
            {
                if ( x <= 0 || x >= nObjs ) continue;
                Gia_Obj_t * o = Gia_ManObj( p, x );
                int lv = Gia_ObjLevel( p, o );
                if ( lv < mn ) { mn = lv; has = 1; }
            }
            hasMin = has;
            minLv  = has ? mn : 0;
        }

        Vec_Bit_t * vMidLocal    = Vec_BitStart( nObjs );
        Vec_Bit_t * vBrdNew      = Vec_BitStart( nObjs );
        Vec_Bit_t * vLevelBrdNew = Vec_BitStart( nObjs );

        ComputeU2ForPo_POInputs_UsingGiaLevel_GlobalBorder(
            p, poId, vPoB, vBorderAll, minLv, hasMin, vMidLocal, vBrdNew, vLevelBrdNew
        );

        for (int j = 0; j < nObjs; ++j)
            if ( Vec_BitEntry(vBrdNew, j) )
                Vec_BitWriteEntry(vBorderAll, j, 1);

        for (int j = 0; j < nObjs; ++j)
            if ( Vec_BitEntry(vMidLocal, j) && !Vec_BitEntry(vBorderAll, j) )
                Vec_BitWriteEntry(vMiddleAll, j, 1);

        for (int j = 0; j < nObjs; ++j)
            if ( Vec_BitEntry(vBorderAll, j) && Vec_BitEntry(vMiddleAll, j) )
                Vec_BitWriteEntry(vMiddleAll, j, 0);

        Vec_BitFree( vPoB );
        Vec_BitFree( vMidLocal );
        Vec_BitFree( vBrdNew );
        Vec_BitFree( vLevelBrdNew );
    }

    *ppBorderAll = vBorderAll;
    *ppMiddleAll = vMiddleAll;
}

// --------------------- Entry ---------------------
extern Vec_Int_t * Ree_ManComputeCuts( Gia_Man_t * p, Vec_Int_t ** pvXors, int fVerbose );

Vec_Int_t * Gia_PoLocalKCover_reverse( Gia_Man_t * p, int k, Vec_Int_t ** ppBoundaryNodes )
{
    *ppBoundaryNodes = NULL;

    if ( p == NULL ) {
        *ppBoundaryNodes = Vec_IntAlloc(0);
        return Vec_IntAlloc(0);
    }
    if ( k < 0 ) k = 0;
    int nObjs = Gia_ManObjNum(p);

    // A) Build HA/FA info
    Vec_Int_t *vXors = NULL;
    Vec_Int_t *vAdds = Ree_ManComputeCuts( p, &vXors, 0 );
    if (!vAdds) {
        if (vXors) Vec_IntFree(vXors);
        *ppBoundaryNodes = Vec_IntAlloc(0);
        return Vec_IntAlloc(0);
    }
    RMaps2 *M = BuildMapsFromAdds2( p, vAdds );

    Vec_Bit_t * isPairHit      = Vec_BitStart( M->pairCnt );
    Vec_Bit_t * bannedHitNodes = Vec_BitStart( nObjs );

    // B) For each PO, pick first-hit pair
    Vec_Bit_t * vBoundaryInputs = Vec_BitStart( nObjs );
    Vec_Bit_t * vHitAllBits     = Vec_BitStart( nObjs );
    Vec_Wec_t * poInputs        = Vec_WecStart( Gia_ManPoNum(p) );

    Gia_Obj_t *pPo; int i;
    Gia_ManForEachPo( p, pPo, i ) {
        int poId = Gia_ObjId( p, pPo );
        int rootId = -1, pairId = -1;

        FindPoFirstHitPair2_PrefNew( p, poId, M, isPairHit, bannedHitNodes, &rootId, &pairId );

        if (pairId >= 0) {
            Vec_BitWriteEntry(isPairHit, pairId, 1);

            Vec_Int_t *vIn = Vec_WecEntry( M->pairInputs, pairId );
            Vec_Int_t *vRt = Vec_WecEntry( M->pairRoots,  pairId );

            int t, x;
            Vec_IntForEachEntry(vIn, x, t) {
                if (x > 0 && x < nObjs) Vec_BitWriteEntry(vBoundaryInputs, x, 1);
            }

            // Save per-PO inputs
            Vec_Int_t * vSlot = Vec_WecEntry( poInputs, i );
            int tt, xx;
            Vec_IntForEachEntry( vIn, xx, tt )
                if ( xx > 0 && xx < nObjs )
                    Vec_IntPush( vSlot, xx );
            Vec_IntSort( vSlot, 0 );
            Vec_IntUniqify( vSlot );

            // Pair nodes for this hit
            Vec_Bit_t * vTmpHit = Vec_BitStart( nObjs );
            (void)CollectPairConeNodes( p, vRt, vBoundaryInputs, vTmpHit );

            // Merge into global hit set and ban-list
            (void)CollectPairConeNodes( p, vRt, vBoundaryInputs, vHitAllBits );
            for ( int j = 0; j < Vec_BitSize(vTmpHit); ++j )
                if ( Vec_BitEntry(vTmpHit, j) )
                    Vec_BitWriteEntry( bannedHitNodes, j, 1 );

            Vec_BitFree( vTmpHit );
        }
    }

    Vec_Int_t *vBoundaryList = ExportBitsetIds_All( vBoundaryInputs );

    // U2 computation
    Vec_Bit_t * vBorderAll = NULL;
    Vec_Bit_t * vMiddleAll = NULL;
    RunClosureByTraversalU2_POInputs( p, poInputs, &vBorderAll, &vMiddleAll );
    (void)FixClosure_AddMissingFaninToBorder( p, vBorderAll, vMiddleAll );

    // Build U2 all nodes (Border ∪ Middle)
    Vec_Bit_t * vU2All = Vec_BitStart( nObjs );
    for (int j = 0; j < nObjs; ++j) {
        if ( Vec_BitEntry(vBorderAll, j) || Vec_BitEntry(vMiddleAll, j) )
            Vec_BitWriteEntry(vU2All, j, 1);
    }

    // Export results
    Vec_Int_t * vFsaNodes = ExportBitsetIds_NoZero( vU2All );
    Vec_Int_t * vBrdNodes = ExportBitsetIds_NoZero( vBorderAll );

    // If boundary is empty: return all nodes as FSA
    if ( Vec_IntSize(vBoundaryList) == 0 )
    {
        Vec_IntFree(vFsaNodes);
        vFsaNodes = Vec_IntAlloc( nObjs );
        Gia_Obj_t * pObj; int id;
        Gia_ManForEachObj( p, pObj, id )
            Vec_IntPush( vFsaNodes, id );
        Vec_IntSort( vFsaNodes, 0 );
        Vec_IntUniqify( vFsaNodes );
    }

    *ppBoundaryNodes = vBrdNodes;

    // Free
    Vec_IntFree(vAdds);
    Vec_IntFree(vXors);
    FreeMaps2(M);
    Vec_BitFree(isPairHit);
    Vec_BitFree(bannedHitNodes);
    Vec_BitFree(vBoundaryInputs);
    Vec_BitFree(vHitAllBits);
    Vec_IntFree(vBoundaryList);
    Vec_WecFree(poInputs);
    Vec_BitFree(vBorderAll);
    Vec_BitFree(vMiddleAll);
    Vec_BitFree(vU2All);

    return vFsaNodes;
}

// Optimized version that reuses pre-computed vAdds from Gia_detecting_ha_fa_ext
// This avoids redundant Ree_ManComputeCuts calls
Vec_Int_t * Gia_PoLocalKCover_reverse_with_adds( Gia_Man_t * p, int k, Vec_Int_t ** ppBoundaryNodes, Vec_Int_t * vAddsIn )
{
    *ppBoundaryNodes = NULL;

    if ( p == NULL ) {
        *ppBoundaryNodes = Vec_IntAlloc(0);
        return Vec_IntAlloc(0);
    }
    if ( k < 0 ) k = 0;
    int nObjs = Gia_ManObjNum(p);

    // Use provided vAdds or compute new one
    Vec_Int_t *vAdds;
    int freeAdds = 0;
    if (vAddsIn) {
        vAdds = Vec_IntDup(vAddsIn);  // Make a copy since we'll modify/free it
        freeAdds = 1;
    } else {
        Vec_Int_t *vXors = NULL;
        vAdds = Ree_ManComputeCuts( p, &vXors, 0 );
        if (!vAdds) {
            if (vXors) Vec_IntFree(vXors);
            *ppBoundaryNodes = Vec_IntAlloc(0);
            return Vec_IntAlloc(0);
        }
        Vec_IntFree(vXors);
        freeAdds = 1;
    }

    RMaps2 *M = BuildMapsFromAdds2( p, vAdds );

    Vec_Bit_t * isPairHit      = Vec_BitStart( M->pairCnt );
    Vec_Bit_t * bannedHitNodes = Vec_BitStart( nObjs );

    // B) For each PO, pick first-hit pair
    Vec_Bit_t * vBoundaryInputs = Vec_BitStart( nObjs );
    Vec_Bit_t * vHitAllBits     = Vec_BitStart( nObjs );
    Vec_Wec_t * poInputs        = Vec_WecStart( Gia_ManPoNum(p) );

    Gia_Obj_t *pPo; int i;
    Gia_ManForEachPo( p, pPo, i ) {
        int poId = Gia_ObjId( p, pPo );
        int rootId = -1, pairId = -1;

        FindPoFirstHitPair2_PrefNew( p, poId, M, isPairHit, bannedHitNodes, &rootId, &pairId );

        if (pairId >= 0) {
            Vec_BitWriteEntry(isPairHit, pairId, 1);

            Vec_Int_t *vIn = Vec_WecEntry( M->pairInputs, pairId );
            Vec_Int_t *vRt = Vec_WecEntry( M->pairRoots,  pairId );

            int t, x;
            Vec_IntForEachEntry(vIn, x, t) {
                if (x > 0 && x < nObjs) Vec_BitWriteEntry(vBoundaryInputs, x, 1);
            }

            // Save per-PO inputs
            Vec_Int_t * vSlot = Vec_WecEntry( poInputs, i );
            int tt, xx;
            Vec_IntForEachEntry( vIn, xx, tt )
                if ( xx > 0 && xx < nObjs )
                    Vec_IntPush( vSlot, xx );
            Vec_IntSort( vSlot, 0 );
            Vec_IntUniqify( vSlot );

            // Pair nodes for this hit
            Vec_Bit_t * vTmpHit = Vec_BitStart( nObjs );
            (void)CollectPairConeNodes( p, vRt, vBoundaryInputs, vTmpHit );

            // Merge into global hit set and ban-list
            (void)CollectPairConeNodes( p, vRt, vBoundaryInputs, vHitAllBits );
            for ( int j = 0; j < Vec_BitSize(vTmpHit); ++j )
                if ( Vec_BitEntry(vTmpHit, j) )
                    Vec_BitWriteEntry( bannedHitNodes, j, 1 );

            Vec_BitFree( vTmpHit );
        }
    }

    Vec_Int_t *vBoundaryList = ExportBitsetIds_All( vBoundaryInputs );

    // U2 computation
    Vec_Bit_t * vBorderAll = NULL;
    Vec_Bit_t * vMiddleAll = NULL;
    RunClosureByTraversalU2_POInputs( p, poInputs, &vBorderAll, &vMiddleAll );
    (void)FixClosure_AddMissingFaninToBorder( p, vBorderAll, vMiddleAll );

    // Build U2 all nodes (Border ∪ Middle)
    Vec_Bit_t * vU2All = Vec_BitStart( nObjs );
    for (int j = 0; j < nObjs; ++j) {
        if ( Vec_BitEntry(vBorderAll, j) || Vec_BitEntry(vMiddleAll, j) )
            Vec_BitWriteEntry(vU2All, j, 1);
    }

    // Export results
    Vec_Int_t * vFsaNodes = ExportBitsetIds_NoZero( vU2All );
    Vec_Int_t * vBrdNodes = ExportBitsetIds_NoZero( vBorderAll );

    // If boundary is empty: return all nodes as FSA
    if ( Vec_IntSize(vBoundaryList) == 0 )
    {
        Vec_IntFree(vFsaNodes);
        vFsaNodes = Vec_IntAlloc( nObjs );
        Gia_Obj_t * pObj; int id;
        Gia_ManForEachObj( p, pObj, id )
            Vec_IntPush( vFsaNodes, id );
        Vec_IntSort( vFsaNodes, 0 );
        Vec_IntUniqify( vFsaNodes );
    }

    *ppBoundaryNodes = vBrdNodes;

    // Free
    if (freeAdds) Vec_IntFree(vAdds);
    FreeMaps2(M);
    Vec_BitFree(isPairHit);
    Vec_BitFree(bannedHitNodes);
    Vec_BitFree(vBoundaryInputs);
    Vec_BitFree(vHitAllBits);
    Vec_IntFree(vBoundaryList);
    Vec_WecFree(poInputs);
    Vec_BitFree(vBorderAll);
    Vec_BitFree(vMiddleAll);
    Vec_BitFree(vU2All);

    return vFsaNodes;
}

// ======================== Reverse-cut local K-cover with reporting (enhanced, crash-safe, prefer-new-pairs) ========================

// ======================== Reverse-cut local K-cover with reporting (enhanced, crash-safe, prefer-new-pairs) ========================

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


ABC_NAMESPACE_IMPL_END

