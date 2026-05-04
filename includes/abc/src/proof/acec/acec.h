/**CFile****************************************************************

  FileName    [acec.h]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [CEC for arithmetic circuits.]

  Synopsis    [External declarations.]

  Author      [Alan Mishchenko]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - June 20, 2005.]

  Revision    [$Id: acec.h,v 1.00 2005/06/20 00:00:00 alanmi Exp $]

***********************************************************************/

#ifndef ABC__proof__acec__acec_h
#define ABC__proof__acec__acec_h


////////////////////////////////////////////////////////////////////////
///                          INCLUDES                                ///
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
///                         PARAMETERS                               ///
////////////////////////////////////////////////////////////////////////

typedef struct Vec_Int_t_ Vec_Int_t;
typedef struct Vec_Wec_t_ Vec_Wec_t;
typedef struct Vec_Wrd_t_ Vec_Wrd_t;

ABC_NAMESPACE_HEADER_START

// C-compatible structures for HA/FA/XOR detection results
// These match the C++ structures in parser.hpp for direct memory passing

// Half Adder structure
typedef struct GiaHalfAdder_ {
    int      id;
    int*     nodes;              // node IDs array
    int      nodes_count;        // number of nodes
    unsigned char* tt;           // truth table bytes
    int      tt_count;           // number of tt bytes (typically 2)
    int*     input_phase;        // input phase array
    int      input_phase_count;  // number of input phases (typically 2)
    int*     output_phase;       // output phase array
    int      output_phase_count; // number of output phases (typically 2)
} GiaHalfAdder;

// Full Adder structure
typedef struct GiaFullAdder_ {
    int      id;
    int*     nodes;              // node IDs array
    int      nodes_count;        // number of nodes
    unsigned char* tt;           // truth table bytes
    int      tt_count;           // number of tt bytes (typically 2)
    int*     input_phase;        // input phase array
    int      input_phase_count;  // number of input phases (typically 3)
    int*     output_phase;       // output phase array
    int      output_phase_count; // number of output phases (typically 2)
} GiaFullAdder;

// XOR Remaining structure
typedef struct GiaXorRemain_ {
    int      id;
    int*     nodes;              // node IDs array
    int      nodes_count;        // number of nodes
    unsigned char* tt;           // truth table bytes
    int      tt_count;           // number of tt bytes (typically 1)
    int*     input_phase;        // input phase array
    int      input_phase_count;  // number of input phases (2 or 3)
    int*     output_phase;       // output phase array
    int      output_phase_count; // number of output phases (typically 1)
} GiaXorRemain;

typedef struct GiaDetectingResults_ {
  Vec_Int_t *xor2_list;
  Vec_Int_t *xor3_list;
  Vec_Int_t *maj2_list;
  Vec_Int_t *maj3_list;
  Vec_Int_t *xor_all_list;
  Vec_Int_t *xor_remaining;
  Vec_Int_t *xor2_maj2_combined;
  Vec_Int_t *xor3_maj3_combined;
  Vec_Int_t *ha_roots;
  Vec_Int_t *fa_roots;
  Vec_Int_t *all_roots;

  // New structured output for direct C++/Talisman consumption
  GiaHalfAdder* half_adders;      // array of half adders
  int           half_adders_count;
  GiaFullAdder* full_adders;      // array of full adders
  int           full_adders_count;
  GiaXorRemain* xor_remains;      // array of xor remaining
  int           xor_remains_count;
} GiaDetectingResults;

// Linear/affine cut collection (single-output affine functions, k <= 6)
typedef struct Ree_LinearRes_t_
{
  Vec_Int_t * vRoots;   // root node IDs, one per recorded cut
  Vec_Wec_t * vLeaves;  // leaves for each cut (parallel to vRoots)
  Vec_Wrd_t * vTruths;  // truth tables for each cut
  Vec_Wec_t * vNodes;   // internal nodes between leaves and root
} Ree_LinearRes_t;
////////////////////////////////////////////////////////////////////////
///                         BASIC TYPES                              ///
////////////////////////////////////////////////////////////////////////

// combinational equivalence checking parameters
typedef struct Acec_ParCec_t_ Acec_ParCec_t;
struct Acec_ParCec_t_
{
    int              nBTLimit;      // conflict limit at a node
    int              TimeLimit;     // the runtime limit in seconds
    int              fMiter;        // input circuit is a miter
    int              fDualOutput;   // dual-output miter
    int              fTwoOutput;    // two-output miter
    int              fBooth;        // expecting Booth multiplier
    int              fSilent;       // print no messages
    int              fVeryVerbose;  // verbose stats
    int              fVerbose;      // verbose stats
    int              iOutFail;      // the number of failed output
};

////////////////////////////////////////////////////////////////////////
///                      MACRO DEFINITIONS                           ///
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
///                             ITERATORS                            ///
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
///                    FUNCTION DECLARATIONS                         ///
////////////////////////////////////////////////////////////////////////

/*=== acecCl.c ========================================================*/
extern Gia_Man_t *   Acec_ManDecla( Gia_Man_t * pGia, int fBooth, int fVerbose );
/*=== acecCore.c ========================================================*/
extern void          Acec_ManCecSetDefaultParams( Acec_ParCec_t * p );
extern int           Acec_Solve( Gia_Man_t * pGia0, Gia_Man_t * pGia1, Acec_ParCec_t * pPars );
/*=== acecFadds.c ========================================================*/
extern Vec_Int_t *   Gia_ManDetectFullAdders( Gia_Man_t * p, int fVerbose, Vec_Int_t ** vCutsXor2 );
extern Vec_Int_t *   Gia_ManDetectHalfAdders( Gia_Man_t * p, int fVerbose );
/*=== acecOrder.c ========================================================*/
extern Vec_Int_t *   Gia_PolynReorder( Gia_Man_t * pGia, int fVerbose, int fVeryVerbose );
extern Vec_Int_t *   Gia_PolynFindOrder( Gia_Man_t * pGia, Vec_Int_t * vFadds, Vec_Int_t * vHadds, int fVerbose, int fVeryVerbose );
/*=== acecPolyn.c ========================================================*/
extern void          Gia_PolynBuild( Gia_Man_t * pGia, Vec_Int_t * vOrder, int fSigned, int fVerbose, int fVeryVerbose );
/*=== acecRe.c ========================================================*/
extern Vec_Int_t *   Ree_ManComputeCuts( Gia_Man_t * p, Vec_Int_t ** pvXors, int fVerbose );
extern int           Ree_ManCountFadds( Vec_Int_t * vAdds );
extern void          Ree_ManPrintAdders( Vec_Int_t * vAdds, int fVerbose );
/*=== acecLinear.c ======================================================*/
extern Ree_LinearRes_t * Ree_ManComputeLinearCuts( Gia_Man_t * p, int fVerbose );
extern void             Ree_LinearResFree( Ree_LinearRes_t * pRes );
/*=== acecTree.c ========================================================*/
extern Gia_Man_t *   Acec_Normalize( Gia_Man_t * pGia, int fBooth, int fVerbose );
extern GiaDetectingResults *Gia_detecting_ha_fa(Gia_Man_t * pGia) ;
// Extended version that also returns vAdds and vXors for reuse
extern GiaDetectingResults *Gia_detecting_ha_fa_ext(Gia_Man_t * pGia, Vec_Int_t ** ppvAdds, Vec_Int_t ** ppvXors);
extern void Gia_free_detecting_results(GiaDetectingResults *results);
extern Vec_Int_t * CollectRemainingNodesAfterTfiRemoval(Gia_Man_t *p);
extern void Vec_IntSortUnique(Vec_Int_t *v);
extern void DumpMarkedTfi(Gia_Man_t *p, const char *tag);
extern  void MarkTfi_FromRoots(Gia_Man_t *p, Vec_Int_t *vRoots);
extern void UniqueIdsInPlace(Gia_Man_t *pMan, Vec_Int_t *v);
extern  Vec_Int_t* Gia_PoLocalKCover_reverse( Gia_Man_t * p, int k, Vec_Int_t ** ppBoundaryNodes );
// Optimized version that reuses pre-computed vAdds from Gia_detecting_ha_fa_ext
extern  Vec_Int_t* Gia_PoLocalKCover_reverse_with_adds( Gia_Man_t * p, int k, Vec_Int_t ** ppBoundaryNodes, Vec_Int_t * vAdds );
ABC_NAMESPACE_HEADER_END


#endif

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////
