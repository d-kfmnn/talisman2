/**CFile****************************************************************

  FileName    [acecLinear.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [CEC for arithmetic circuits.]

  Synopsis    [Detection of affine (linear) cuts with k <= 6.]

  Author      [Codex]
  Affiliation [OpenAI]
  Date        [Ver. 1.2. Updated - Jan 05, 2026.]

  Notes (important changes vs naive version):
    - Adds trivial singleton cut for each AND node (robust connectivity).
    - Maintains canonical sorted unique leaves for every cut.
    - Adds per-node cut limit and containment-based pruning.
    - Affine check accepts both Truth and its complemented version (Truth^mask),
      which makes detection robust to polarity conventions.

***********************************************************************/

#include "acecInt.h"
#include "misc/vec/vecWec.h"
#include "misc/vec/vecWrd.h"
#include "misc/util/utilTruth.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

ABC_NAMESPACE_IMPL_START

////////////////////////////////////////////////////////////////////////
///                        PARAMETERS/DEFS                           ///
////////////////////////////////////////////////////////////////////////

#define REE_LIN_MAX_CUT      6
#define REE_LIN_SEG_SIZE     (REE_LIN_MAX_CUT + 1) // [count, leaf0..leaf5]
#define REE_LIN_CUT_LIMIT    64   // max cuts per node (tune as needed)

#ifndef REE_LIN_USE_DERIV_CHECK
#define REE_LIN_USE_DERIV_CHECK 0
#endif

// Standard 6-var variable truth tables (LSB-first enumeration), for optional recon.
static word s_ReeTruths6[REE_LIN_MAX_CUT] = {
    ABC_CONST(0xAAAAAAAAAAAAAAAA), // x0
    ABC_CONST(0xCCCCCCCCCCCCCCCC), // x1
    ABC_CONST(0xF0F0F0F0F0F0F0F0), // x2
    ABC_CONST(0xFF00FF00FF00FF00), // x3
    ABC_CONST(0xFFFF0000FFFF0000), // x4
    ABC_CONST(0xFFFFFFFF00000000)  // x5
};

////////////////////////////////////////////////////////////////////////
///                     INTERNAL HELPERS                             ///
////////////////////////////////////////////////////////////////////////

static inline word Ree_LinTtMask( int nVars )
{
    if ( nVars >= 6 )
        return ~(word)0;
    return (((word)1) << (1 << nVars)) - 1;
}

#if REE_LIN_USE_DERIV_CHECK
static inline word Ree_LinTtFlipVar( word tt, int nVars, int iVar )
{
    (void)nVars;
    return Abc_Tt6Flip( tt, iVar );
}
static inline int Ree_LinTruthIsAffineDeriv( word tt, int nVars )
{
    int i, j;
    word mask = Ree_LinTtMask( nVars );
    tt &= mask;
    for ( i = 0; i < nVars; i++ )
    for ( j = i + 1; j < nVars; j++ )
    {
        word t_i  = Ree_LinTtFlipVar( tt, nVars, i ) & mask;
        word t_j  = Ree_LinTtFlipVar( tt, nVars, j ) & mask;
        word t_ij = Ree_LinTtFlipVar( t_i, nVars, j ) & mask;
        word d2   = (tt ^ t_i ^ t_j ^ t_ij) & mask;
        if ( d2 != 0 )
            return 0;
    }
    return 1;
}
#endif

static int Ree_LinIntCmp( const void * pa, const void * pb )
{
    int a = *(const int *)pa;
    int b = *(const int *)pb;
    return (a > b) - (a < b);
}

// Collect internal nodes on paths from root to leaves (excludes leaves themselves).
static void Ree_LinCollectNodesBetween( Gia_Man_t * p, Vec_Int_t * vLeaves, int iRoot, Vec_Int_t * vNodes, Vec_Int_t * vStack, Vec_Int_t * vVisited )
{
    Vec_IntClear( vNodes );
    Vec_IntClear( vStack );
    Vec_IntClear( vVisited );
    Vec_IntPush( vStack, iRoot );
    while ( Vec_IntSize(vStack) )
    {
        int iObj = Vec_IntPop( vStack );
        Gia_Obj_t * pObj;
        if ( Vec_IntFind( vVisited, iObj ) != -1 )
            continue;
        Vec_IntPush( vVisited, iObj );
        if ( Vec_IntFind( vLeaves, iObj ) != -1 )
            continue; // stop at leaves
        if ( iObj != iRoot )
            Vec_IntPushUnique( vNodes, iObj );
        pObj = Gia_ManObj( p, iObj );
        if ( pObj == NULL || !Gia_ObjIsAnd( pObj ) )
            continue;
        Vec_IntPush( vStack, Gia_ObjFaninId0p( p, pObj ) );
        Vec_IntPush( vStack, Gia_ObjFaninId1p( p, pObj ) );
    }
}

// Canonicalize leaves: sort ascending and unique.
// Returns new size.
static inline int Ree_LinCanonLeaves( int * pLeaves, int nLeaves )
{
    int i, k;
    if ( nLeaves <= 1 )
        return nLeaves;
    // sort (ascending)
    qsort( pLeaves, (size_t)nLeaves, sizeof(int), Ree_LinIntCmp );
    // unique
    for ( i = 1, k = 1; i < nLeaves; i++ )
        if ( pLeaves[i] != pLeaves[i-1] )
            pLeaves[k++] = pLeaves[i];
    return k;
}

static inline int Ree_LinLeavesAreSortedUnique( int const * pLeaves, int nLeaves )
{
    int i;
    for ( i = 1; i < nLeaves; i++ )
        if ( pLeaves[i-1] >= pLeaves[i] )
            return 0;
    return 1;
}

// Merge two cuts (both must be sorted unique), produce sorted unique union in pOut.
// pCut layout: [nLeaves, leaf0..], with padding ignored.
// Returns nLeaves in union, or 0 if exceeds limit.
static inline int Ree_LinCutMergeOne( int const * pCut0, int const * pCut1, int * pOutLeaves )
{
    int i0 = 0, i1 = 0;
    int n0 = pCut0[0], n1 = pCut1[0];
    int r  = 0;

    // both input leaf arrays must be sorted unique
    // (keep these asserts enabled during bring-up; can be removed for speed)
    assert( n0 >= 0 && n0 <= REE_LIN_MAX_CUT );
    assert( n1 >= 0 && n1 <= REE_LIN_MAX_CUT );
    assert( Ree_LinLeavesAreSortedUnique( pCut0 + 1, n0 ) );
    assert( Ree_LinLeavesAreSortedUnique( pCut1 + 1, n1 ) );

    while ( i0 < n0 || i1 < n1 )
    {
        int v;
        if ( i1 == n1 || (i0 < n0 && pCut0[1+i0] < pCut1[1+i1]) )
            v = pCut0[1+i0++];
        else if ( i0 == n0 || pCut1[1+i1] < pCut0[1+i0] )
            v = pCut1[1+i1++];
        else
            v = pCut0[1+i0], i0++, i1++;

        if ( r && pOutLeaves[r-1] == v )
            continue;
        if ( r == REE_LIN_MAX_CUT )
            return 0;
        pOutLeaves[r++] = v;
    }
    return r;
}

// Push cut into vCuts using fixed segment layout.
// Caller guarantees pLeaves sorted unique of size nLeaves.
static inline void Ree_LinCutPush( Vec_Int_t * vCuts, int nLeaves, int const * pLeaves )
{
    int k;
    Vec_IntPush( vCuts, nLeaves );
    for ( k = 0; k < REE_LIN_MAX_CUT; k++ )
        Vec_IntPush( vCuts, k < nLeaves ? pLeaves[k] : -1 );
}

#define Ree_LinForEachCut( vCuts, pCut, i ) \
    for ( i = 0, pCut = Vec_IntArray(vCuts); i < Vec_IntSize(vCuts); i += REE_LIN_SEG_SIZE, pCut = Vec_IntArray(vCuts) + i )

static inline int Ree_LinCutEqual( int const * pCut, int nLeaves, int const * pLeaves )
{
    int k;
    if ( pCut[0] != nLeaves )
        return 0;
    for ( k = 0; k < nLeaves; k++ )
        if ( pCut[1+k] != pLeaves[k] )
            return 0;
    return 1;
}

// Returns 1 if A ⊆ B (both sorted unique), where sizes nA<=nB.
static inline int Ree_LinIsSubset( int const * A, int nA, int const * B, int nB )
{
    int i = 0, j = 0;
    assert( nA <= nB );
    while ( i < nA && j < nB )
    {
        if ( A[i] == B[j] ) { i++; j++; continue; }
        if ( A[i] >  B[j] ) { j++; continue; }
        return 0; // A[i] < B[j] => missing element
    }
    return i == nA;
}

// Containment pruning + duplication check.
// - If existing cut C ⊆ New, we may want to keep smaller C and drop New.
// - If New ⊆ existing cut C, we drop C (New is better or equal for enumeration).
// In practice, for truth-based checking, smaller cuts are preferable, but we still
// keep multiple sizes; containment pruning reduces redundancy.
static inline int Ree_LinCutInsertPrune( Vec_Int_t * vCuts, int nLeaves, int const * pLeaves )
{
    int * pCut;
    int i, k;
    int idx = 0; // write index in segments

    // First pass: detect duplicate, and also decide if New is dominated by any existing cut.
    Ree_LinForEachCut( vCuts, pCut, i )
    {
        int n = pCut[0];
        if ( n == nLeaves && Ree_LinCutEqual( pCut, nLeaves, pLeaves ) )
            return 0; // duplicate => do not insert

        // If there exists C ⊆ New, prefer keeping C and drop New.
        // Only check when n <= nLeaves.
        if ( n <= nLeaves )
        {
            int tmp[REE_LIN_MAX_CUT];
            for ( k = 0; k < n; k++ ) tmp[k] = pCut[1+k];
            if ( Ree_LinIsSubset( tmp, n, pLeaves, nLeaves ) )
                return 0; // New dominated
        }
    }

    // Second pass: remove existing cuts that are supersets of New (New ⊆ C).
    // We compact in place by copying kept segments forward.
    idx = 0;
    Ree_LinForEachCut( vCuts, pCut, i )
    {
        int n = pCut[0];
        int drop = 0;
        if ( nLeaves <= n )
        {
            int tmp[REE_LIN_MAX_CUT];
            for ( k = 0; k < n; k++ ) tmp[k] = pCut[1+k];
            if ( Ree_LinIsSubset( pLeaves, nLeaves, tmp, n ) )
                drop = 1;
        }
        if ( drop )
            continue;

        if ( idx != i )
            memmove( Vec_IntArray(vCuts) + idx, pCut, REE_LIN_SEG_SIZE * sizeof(int) );
        idx += REE_LIN_SEG_SIZE;
    }
    if ( idx != Vec_IntSize(vCuts) )
        Vec_IntShrink( vCuts, idx );

    // Enforce limit: if too many cuts, do not insert (simple policy).
    // (More sophisticated: keep smallest cuts, or random drop, etc.)
    if ( Vec_IntSize(vCuts) / REE_LIN_SEG_SIZE >= REE_LIN_CUT_LIMIT )
        return 0;

    Ree_LinCutPush( vCuts, nLeaves, pLeaves );
    return 1;
}

// Build sorted table of all affine truth tables for nVars.
static Vec_Wrd_t * Ree_LinBuildAffineTable( int nVars )
{
    int maskLimit = 1 << nVars;
    int expected  = 1 << (nVars + 1);
    Vec_Wrd_t * vRes = Vec_WrdAlloc( expected );
    word m = Ree_LinTtMask( nVars );
    int cMask, cConst;

    for ( cConst = 0; cConst < 2; cConst++ )
    for ( cMask = 0; cMask < maskLimit; cMask++ )
    {
        word Truth = cConst ? m : 0;
        int v;
        for ( v = 0; v < nVars; v++ )
            if ( (cMask >> v) & 1 )
                Truth ^= (s_ReeTruths6[v] & m);
        Truth &= m;
        Vec_WrdPush( vRes, Truth );
    }
    Vec_WrdSort( vRes, 0 );
    assert( Vec_WrdSize(vRes) == expected );
    return vRes;
}

static Vec_Wrd_t ** Ree_LinTables()
{
    static Vec_Wrd_t * vTables[REE_LIN_MAX_CUT + 1];
    static int fInit = 0;
    int k;
    if ( fInit )
        return vTables;
    for ( k = 1; k <= REE_LIN_MAX_CUT; k++ )
        vTables[k] = Ree_LinBuildAffineTable( k );
    fInit = 1;
    return vTables;
}

static inline int Ree_LinTruthIsAffineTable( Vec_Wrd_t * vTable, word Truth )
{
    int Left = 0, Right = Vec_WrdSize(vTable) - 1;
    while ( Left <= Right )
    {
        int Mid = (Left + Right) >> 1;
        word Value = Vec_WrdEntry( vTable, Mid );
        if ( Truth < Value )
            Right = Mid - 1;
        else if ( Truth > Value )
            Left = Mid + 1;
        else
            return 1;
    }
    return 0;
}

static inline void Ree_LinRecordCut( Ree_LinearRes_t * pRes, int iRoot, int const * pLeaves, int nLeaves, word Truth, Vec_Int_t * vNodes )
{
    int k, id;
    Vec_IntPush( pRes->vRoots, iRoot );
    Vec_WrdPush( pRes->vTruths, Truth );
    {
        Vec_Int_t * vCut = Vec_WecPushLevel( pRes->vLeaves );
        for ( k = 0; k < nLeaves; k++ )
            Vec_IntPush( vCut, pLeaves[k] );
    }
    {
        Vec_Int_t * vMid = Vec_WecPushLevel( pRes->vNodes );
        Vec_IntForEachEntry( vNodes, id, k )
            Vec_IntPush( vMid, id );
    }
}

static void Ree_LinCheckCuts( Gia_Man_t * p, int iObj, Vec_Int_t * vCuts, Vec_Int_t * vLeafTemp, Vec_Wrd_t ** vTables, Ree_LinearRes_t * pRes, Vec_Int_t * vNodeTemp, Vec_Int_t * vStack, Vec_Int_t * vVisited )
{
    int * pCut, i, k;
    Gia_Obj_t * pObj = Gia_ManObj( p, iObj );

    // We can analyze CIs and ANDs; others usually have no cuts anyway.
    if ( !Gia_ObjIsAnd(pObj) && !Gia_ObjIsCi(pObj) )
        return;

    Ree_LinForEachCut( vCuts, pCut, i )
    {
        int nLeaves = pCut[0];
        if ( nLeaves <= 0 || nLeaves > REE_LIN_MAX_CUT )
            continue;
        if ( nLeaves <= 1 )
            continue; // skip singleton cuts

        Vec_IntClear( vLeafTemp );
        for ( k = 0; k < nLeaves; k++ )
            Vec_IntPush( vLeafTemp, pCut[1+k] );

        // compute cut truth
        {
            word Truth = *Gia_ObjComputeTruthTableCut( p, pObj, vLeafTemp );
            word m = Ree_LinTtMask( nLeaves );
            Truth &= m;

            // polarity-robust affine check: Truth or complemented Truth is affine
            if ( !Ree_LinTruthIsAffineTable( vTables[nLeaves], Truth ) )
            {
                word TruthC = (Truth ^ m) & m;
                if ( !Ree_LinTruthIsAffineTable( vTables[nLeaves], TruthC ) )
                    continue;
                Truth = TruthC;
            }

#if REE_LIN_USE_DERIV_CHECK
            if ( !Ree_LinTruthIsAffineDeriv( Truth, nLeaves ) )
                continue;
#endif
            Ree_LinCollectNodesBetween( p, vLeafTemp, iObj, vNodeTemp, vStack, vVisited );
            Ree_LinRecordCut( pRes, iObj, Vec_IntArray(vLeafTemp), nLeaves, Truth, vNodeTemp );
        }
    }
}

////////////////////////////////////////////////////////////////////////
///                        EXPORTED FUNCTIONS                         ///
////////////////////////////////////////////////////////////////////////

Ree_LinearRes_t * Ree_ManComputeLinearCuts( Gia_Man_t * p, int fVerbose )
{
    Vec_Wec_t * vNodeCuts = Vec_WecStart( Gia_ManObjNum(p) );
    Vec_Wrd_t ** vTables = Ree_LinTables();
    Vec_Int_t * vTempLeaves = Vec_IntAlloc( REE_LIN_MAX_CUT );
    Vec_Int_t * vTempNodes  = Vec_IntAlloc( 64 );
    Vec_Int_t * vTempStack  = Vec_IntAlloc( 64 );
    Vec_Int_t * vTempVisit  = Vec_IntAlloc( 64 );
    Ree_LinearRes_t * pRes = ABC_ALLOC( Ree_LinearRes_t, 1 );
    Gia_Obj_t * pObj;
    int i, * pCut0, * pCut1;
    int TmpLeaves[REE_LIN_MAX_CUT];
    int nRecorded = 0, nVisited = 0, nCutsAnd = 0;

    pRes->vRoots  = Vec_IntAlloc( 256 );
    pRes->vTruths = Vec_WrdAlloc( 256 );
    pRes->vLeaves = Vec_WecAlloc( 256 );
    pRes->vNodes  = Vec_WecAlloc( 256 );

    // Prepare truth infrastructure for k <= REE_LIN_MAX_CUT.
    Gia_ObjComputeTruthTableStart( p, REE_LIN_MAX_CUT );

    // Initialize CI cuts: singleton {ci} (canonical by definition)
    Gia_ManForEachCi( p, pObj, i )
    {
        int leaf = Gia_ObjId( p, pObj );
        Vec_Int_t * vCuts = Vec_WecEntry( vNodeCuts, leaf );
        TmpLeaves[0] = leaf;
        Ree_LinCutInsertPrune( vCuts, 1, TmpLeaves );
    }

    // (Optional) Constant object: we do not create cuts for const0 by default.
    // If you want, you can add empty cut {} but Gia_ObjComputeTruthTableCut expects leaves,
    // so keep it simple and skip.

    // Iterate AND nodes in topo order; build cuts by merging fanin cuts.
    Gia_ManForEachAnd( p, pObj, i )
    {
        Vec_Int_t * vCuts  = Vec_WecEntry( vNodeCuts, i );
        Vec_Int_t * vCuts0 = Vec_WecEntry( vNodeCuts, Gia_ObjFaninId0(pObj, i) );
        Vec_Int_t * vCuts1 = Vec_WecEntry( vNodeCuts, Gia_ObjFaninId1(pObj, i) );
        int i0, i1;

        // IMPORTANT: add trivial cut {i} like acecRe.c does (robustness)
        TmpLeaves[0] = i;
        Ree_LinCutInsertPrune( vCuts, 1, TmpLeaves );

        // Merge fanin cuts
        Ree_LinForEachCut( vCuts0, pCut0, i0 )
        Ree_LinForEachCut( vCuts1, pCut1, i1 )
        {
            int nLeaves = Ree_LinCutMergeOne( pCut0, pCut1, TmpLeaves );
            if ( nLeaves == 0 )
                continue;

            // Canonicalize defensively (should already be sorted unique from merge)
            nLeaves = Ree_LinCanonLeaves( TmpLeaves, nLeaves );
            if ( nLeaves == 0 || nLeaves > REE_LIN_MAX_CUT )
                continue;
            assert( Ree_LinLeavesAreSortedUnique( TmpLeaves, nLeaves ) );

            Ree_LinCutInsertPrune( vCuts, nLeaves, TmpLeaves );
        }

        nCutsAnd += Vec_IntSize(vCuts) / REE_LIN_SEG_SIZE;
    }

    // Analyze all nodes: compute cut truth and test affine
    Gia_ManForEachObj( p, pObj, i )
    {
        Vec_Int_t * vCuts = Vec_WecEntry( vNodeCuts, i );
        if ( Vec_IntSize( vCuts ) == 0 )
            continue;
        nVisited += Vec_IntSize( vCuts ) / REE_LIN_SEG_SIZE;
        Ree_LinCheckCuts( p, i, vCuts, vTempLeaves, vTables, pRes, vTempNodes, vTempStack, vTempVisit );
    }

    nRecorded = Vec_IntSize( pRes->vRoots );

    Gia_ObjComputeTruthTableStop( p );
    Vec_WecFree( vNodeCuts );
    Vec_IntFree( vTempLeaves );
    Vec_IntFree( vTempNodes );
    Vec_IntFree( vTempStack );
    Vec_IntFree( vTempVisit );

    if ( fVerbose )
    {
        printf( "Affine cuts (k<=%d): visited %d, recorded %d.  AvgCuts/AND = %.2f (limit %d)\n",
            REE_LIN_MAX_CUT, nVisited, nRecorded,
            Gia_ManAndNum(p) ? 1.0*nCutsAnd/Gia_ManAndNum(p) : 0.0,
            REE_LIN_CUT_LIMIT );
    }
    return pRes;
}

void Ree_LinearResFree( Ree_LinearRes_t * pRes )
{
    if ( pRes == NULL )
        return;
    Vec_IntFree( pRes->vRoots );
    Vec_WrdFree( pRes->vTruths );
    Vec_WecFree( pRes->vLeaves );
    Vec_WecFree( pRes->vNodes );
    ABC_FREE( pRes );
}

ABC_NAMESPACE_IMPL_END
