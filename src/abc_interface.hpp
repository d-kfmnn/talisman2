// abc_simple_interface.hpp
#ifndef ABC_SIMPLE_INTERFACE_HPP
#define ABC_SIMPLE_INTERFACE_HPP

// ABC is a C library, needs extern "C"
extern "C" {
    #include "misc/util/abc_global.h"
    #include "misc/vec/vec.h"
    #include "misc/vec/vecWec.h"
    #include "misc/vec/vecWrd.h"
    #include "aig/gia/gia.h"
    #include "base/main/main.h"
    #include "base/main/mainInt.h"
    #include "proof/acec/acec.h"
}

#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cstdint>

// Include talisman parser structures
#include "parser.hpp"

namespace abc_interface {

// ======================================================================
// Cached vAdds for reuse between detect_ha_fa and get_kcover_nodes
// ======================================================================
struct DetectionCache {
    Vec_Int_t* vAdds = nullptr;
    Vec_Int_t* vXors = nullptr;

    void clear() {
        if (vAdds) { Vec_IntFree(vAdds); vAdds = nullptr; }
        if (vXors) { Vec_IntFree(vXors); vXors = nullptr; }
    }

    ~DetectionCache() { clear(); }
};

// Global cache for reuse (thread-local if needed)
inline DetectionCache& get_detection_cache() {
    static DetectionCache cache;
    return cache;
}

// ======================================================================
// Fast Möbius Transform: Truth Table -> ANF (Algebraic Normal Form)
// ======================================================================
// Based on the algorithm from "Boolean Functions for Cryptography and Coding Theory"
// by Claude Carlet. For an n-variable boolean function:
//   ANF: f(x) = ⊕_{I⊆{1,...,n}} a_I * x_I
// where x_I = ∏_{i∈I} x_i
//
// The Fast Möbius Transform computes ANF coefficients from truth table:
//   a_I = ⊕_{x: supp(x)⊆I} f(x)
//
// For AFFINE functions, the ANF has only constant and linear terms:
//   f(x) = a_0 ⊕ a_1*x_1 ⊕ a_2*x_2 ⊕ ... ⊕ a_n*x_n

// In-place Fast Möbius Transform
// Input: truth table tt of size 2^nVars (as array or vector)
// Output: ANF coefficients in same array
// After transform: anf[u] = coefficient a_{supp(u)} in ANF
inline void fast_mobius_transform(std::vector<int>& tt, int nVars) {
    int size = 1 << nVars;
    for (int i = 0; i < nVars; i++) {
        int step = 1 << i;
        for (int j = 0; j < size; j++) {
            if ((j >> i) & 1) {
                tt[j] ^= tt[j ^ step];
            }
        }
    }
}

// Convert 64-bit truth table to ANF coefficients array
// For k-variable function, truth table has 2^k entries
// Returns: vector of size 2^k where result[i] is coefficient for monomial x^i
inline std::vector<int> truth_to_anf(uint64_t truth, int nVars) {
    int size = 1 << nVars;
    std::vector<int> anf(size);

    // Extract truth table bits
    for (int i = 0; i < size; i++) {
        anf[i] = (truth >> i) & 1;
    }

    // Apply Fast Möbius Transform
    fast_mobius_transform(anf, nVars);

    return anf;
}

// Extract linear polynomial coefficients from ANF of an affine function
// For affine function: f(x) = a_0 ⊕ a_1*x_1 ⊕ ... ⊕ a_n*x_n
// ANF indices: 0 -> constant, 1 -> x_1, 2 -> x_2, 4 -> x_3, 8 -> x_4, etc.
// Returns: (constant_term, [coeff_x1, coeff_x2, ...])
inline void extract_affine_coeffs(
    const std::vector<int>& anf,
    int nVars,
    int& constant_coeff,
    std::vector<int>& linear_coeffs
) {
    constant_coeff = (anf.size() > 0) ? anf[0] : 0;
    linear_coeffs.clear();
    linear_coeffs.reserve(nVars);

    for (int i = 0; i < nVars; i++) {
        int idx = 1 << i;  // Index for variable x_i (2^i)
        linear_coeffs.push_back((idx < (int)anf.size()) ? anf[idx] : 0);
    }
}

// Compute ANF polynomial for a LinearCut and populate its coefficient fields
// The polynomial represents: root = constant_coeff ⊕ Σ(linear_coeffs[i] * leaves[i])
inline void compute_linear_cut_polynomial(talisman::LinearCut& cut) {
    int nVars = static_cast<int>(cut.leaves.size());
    if (nVars == 0 || nVars > 6) {
        cut.constant_coeff = 0;
        cut.linear_coeffs.clear();
        return;
    }

    // Convert truth table to ANF
    std::vector<int> anf = truth_to_anf(cut.truth, nVars);

    // Extract affine coefficients
    extract_affine_coeffs(anf, nVars, cut.constant_coeff, cut.linear_coeffs);
}

// Verify that a truth table represents an affine function
// An affine function has only constant and degree-1 terms in its ANF
inline bool is_affine_anf(const std::vector<int>& anf, int nVars) {
    int size = 1 << nVars;
    for (int i = 0; i < size; i++) {
        // Check if this index corresponds to a non-linear term
        // Non-linear terms have popcount(i) > 1
        int popcount = __builtin_popcount(i);
        if (popcount > 1 && anf[i] != 0) {
            return false;  // Has non-linear term
        }
    }
    return true;
}

// Initialize ABC framework
inline void abc_init() {
    Abc_Start();
}

// Cleanup ABC framework
inline void abc_cleanup() {
    Abc_Stop();
}

// Load AIG file, return Gia_Man_t* pointer
inline Gia_Man_t* load_aig(const std::string& filename) {
    Gia_Man_t* pGia = Gia_AigerRead(const_cast<char*>(filename.c_str()), 0, 0, 0);
    if (!pGia) {
        throw std::runtime_error("Failed to load AIG file: " + filename);
    }
    return pGia;
}

// Free Gia_Man_t
inline void free_gia(Gia_Man_t* pGia) {
    if (pGia) {
        Gia_ManStop(pGia);
    }
}

// ======================================================================
// Conversion functions: GiaDetectingResults -> talisman structures
// ======================================================================

// Convert GiaHalfAdder array to std::vector<talisman::HalfAdder>
inline void convert_half_adders(
    const GiaHalfAdder* src,
    int count,
    std::vector<talisman::HalfAdder>& dst
) {
    dst.clear();
    dst.reserve(count);

    for (int i = 0; i < count; i++) {
        const GiaHalfAdder& ha = src[i];
        talisman::HalfAdder th;

        th.id = ha.id;

        // Copy nodes
        if (ha.nodes && ha.nodes_count > 0) {
            th.nodes.assign(ha.nodes, ha.nodes + ha.nodes_count);
        }

        // Copy tt
        if (ha.tt && ha.tt_count > 0) {
            th.tt.assign(ha.tt, ha.tt + ha.tt_count);
        }

        // Copy input_phase
        if (ha.input_phase && ha.input_phase_count > 0) {
            th.input_phase.assign(ha.input_phase, ha.input_phase + ha.input_phase_count);
        }

        // Copy output_phase
        if (ha.output_phase && ha.output_phase_count > 0) {
            th.output_phase.assign(ha.output_phase, ha.output_phase + ha.output_phase_count);
        }

        dst.push_back(std::move(th));
    }
}

// Convert GiaFullAdder array to std::vector<talisman::FullAdder>
inline void convert_full_adders(
    const GiaFullAdder* src,
    int count,
    std::vector<talisman::FullAdder>& dst
) {
    dst.clear();
    dst.reserve(count);

    for (int i = 0; i < count; i++) {
        const GiaFullAdder& fa = src[i];
        talisman::FullAdder tf;

        tf.id = fa.id;

        // Copy nodes
        if (fa.nodes && fa.nodes_count > 0) {
            tf.nodes.assign(fa.nodes, fa.nodes + fa.nodes_count);
        }

        // Copy tt
        if (fa.tt && fa.tt_count > 0) {
            tf.tt.assign(fa.tt, fa.tt + fa.tt_count);
        }

        // Copy input_phase
        if (fa.input_phase && fa.input_phase_count > 0) {
            tf.input_phase.assign(fa.input_phase, fa.input_phase + fa.input_phase_count);
        }

        // Copy output_phase
        if (fa.output_phase && fa.output_phase_count > 0) {
            tf.output_phase.assign(fa.output_phase, fa.output_phase + fa.output_phase_count);
        }

        dst.push_back(std::move(tf));
    }
}

// Convert GiaXorRemain array to std::vector<talisman::XorRemain>
inline void convert_xor_remains(
    const GiaXorRemain* src,
    int count,
    std::vector<talisman::XorRemain>& dst
) {
    dst.clear();
    dst.reserve(count);

    for (int i = 0; i < count; i++) {
        const GiaXorRemain& xr = src[i];
        talisman::XorRemain tx;

        tx.id = xr.id;

        // Copy nodes
        if (xr.nodes && xr.nodes_count > 0) {
            tx.nodes.assign(xr.nodes, xr.nodes + xr.nodes_count);
        }

        // Copy tt
        if (xr.tt && xr.tt_count > 0) {
            tx.tt.assign(xr.tt, xr.tt + xr.tt_count);
        }

        // Copy input_phase
        if (xr.input_phase && xr.input_phase_count > 0) {
            tx.input_phase.assign(xr.input_phase, xr.input_phase + xr.input_phase_count);
        }

        // Copy output_phase
        if (xr.output_phase && xr.output_phase_count > 0) {
            tx.output_phase.assign(xr.output_phase, xr.output_phase + xr.output_phase_count);
        }

        dst.push_back(std::move(tx));
    }
}

// ======================================================================
// Main detection function with structured output
// ======================================================================

// Detect HA/FA/XOR and return structured results directly
// This replaces JSON file reading with direct memory passing
// Also caches vAdds for reuse in get_kcover_nodes_optimized
inline bool detect_ha_fa_xor(
    Gia_Man_t* pGia,
    std::vector<talisman::HalfAdder>& half_adders,
    std::vector<talisman::FullAdder>& full_adders,
    std::vector<talisman::XorRemain>& xor_list
) {
    half_adders.clear();
    full_adders.clear();
    xor_list.clear();

    if (!pGia) {
        std::cerr << "Error: GIA is null" << std::endl;
        return false;
    }

    // Clear any previous cache
    DetectionCache& cache = get_detection_cache();
    cache.clear();

    // Use extended version that returns vAdds and vXors for reuse
    Vec_Int_t* vAdds = nullptr;
    Vec_Int_t* vXors = nullptr;
    GiaDetectingResults* pRes = Gia_detecting_ha_fa_ext(pGia, &vAdds, &vXors);

    if (!pRes) {
        std::cerr << "HA/FA detection failed!" << std::endl;
        if (vAdds) Vec_IntFree(vAdds);
        if (vXors) Vec_IntFree(vXors);
        return false;
    }

    // Cache vAdds for later use in get_kcover_nodes_optimized
    cache.vAdds = vAdds;
    cache.vXors = vXors;

    // Convert C structures to C++ vectors
    if (pRes->half_adders && pRes->half_adders_count > 0) {
        convert_half_adders(pRes->half_adders, pRes->half_adders_count, half_adders);
    }

    if (pRes->full_adders && pRes->full_adders_count > 0) {
        convert_full_adders(pRes->full_adders, pRes->full_adders_count, full_adders);
    }

    if (pRes->xor_remains && pRes->xor_remains_count > 0) {
        convert_xor_remains(pRes->xor_remains, pRes->xor_remains_count, xor_list);
    }

    // Free the C results
    Gia_free_detecting_results(pRes);

    if (verbose >= 1) std::cout << "Detected: half_adders=" << half_adders.size()
              << ", full_adders=" << full_adders.size()
              << ", xor_remaining=" << xor_list.size() << std::endl;

    return true;
}

// Detect affine (linear) single-output cuts with k <= 6.
inline bool detect_linear_cuts(
    Gia_Man_t* pGia,
    std::vector<talisman::LinearCut>& cuts,
    bool verbose = false
) {
    cuts.clear();
    if (!pGia) {
        std::cerr << "Error: GIA is null" << std::endl;
        return false;
    }

    Ree_LinearRes_t* pRes = Ree_ManComputeLinearCuts(pGia, 1);
    if (!pRes) {
        std::cerr << "Affine cut detection failed!" << std::endl;
        return false;
    }

    int n = Vec_IntSize(pRes->vRoots);
    cuts.reserve(n);
    for (int i = 0; i < n; i++) {
        talisman::LinearCut cut;
        cut.root = Vec_IntEntry(pRes->vRoots, i);
        Vec_Int_t* vLeaves = Vec_WecEntry(pRes->vLeaves, i);
        cut.leaves.assign(Vec_IntArray(vLeaves), Vec_IntArray(vLeaves) + Vec_IntSize(vLeaves));
        cut.truth = Vec_WrdEntry(pRes->vTruths, i);
        // Extract internal nodes between leaves and root
        Vec_Int_t* vNodes = Vec_WecEntry(pRes->vNodes, i);
        cut.nodes.assign(Vec_IntArray(vNodes), Vec_IntArray(vNodes) + Vec_IntSize(vNodes));

        // Compute ANF polynomial coefficients using Fast Möbius Transform
        compute_linear_cut_polynomial(cut);

        cuts.push_back(std::move(cut));
    }

    Ree_LinearResFree(pRes);

    if (verbose >= 1) std::cout << "Detected affine cuts (k<=6): " << cuts.size() << std::endl;
    return true;
}

// ======================================================================
// Print functions for debugging/display
// ======================================================================

// Print a single HalfAdder
inline void print_half_adder(const talisman::HalfAdder& ha) {
    if (verbose < 1) return;
    std::cout << "  HalfAdder[" << ha.id << "]:" << std::endl;
    std::cout << "    nodes: [";
    for (size_t i = 0; i < ha.nodes.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << ha.nodes[i];
    }
    std::cout << "]" << std::endl;

    std::cout << "    tt: [";
    for (size_t i = 0; i < ha.tt.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << "0x" << std::hex << (int)ha.tt[i] << std::dec;
    }
    std::cout << "]" << std::endl;

    std::cout << "    input_phase: [";
    for (size_t i = 0; i < ha.input_phase.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << ha.input_phase[i];
    }
    std::cout << "]" << std::endl;

    std::cout << "    output_phase: [";
    for (size_t i = 0; i < ha.output_phase.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << ha.output_phase[i];
    }
    std::cout << "]" << std::endl;
}

// Print a single FullAdder
inline void print_full_adder(const talisman::FullAdder& fa) {
    if (verbose < 1) return;
    std::cout << "  FullAdder[" << fa.id << "]:" << std::endl;
    std::cout << "    nodes: [";
    for (size_t i = 0; i < fa.nodes.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << fa.nodes[i];
    }
    std::cout << "]" << std::endl;

    std::cout << "    tt: [";
    for (size_t i = 0; i < fa.tt.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << "0x" << std::hex << (int)fa.tt[i] << std::dec;
    }
    std::cout << "]" << std::endl;

    std::cout << "    input_phase: [";
    for (size_t i = 0; i < fa.input_phase.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << fa.input_phase[i];
    }
    std::cout << "]" << std::endl;

    std::cout << "    output_phase: [";
    for (size_t i = 0; i < fa.output_phase.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << fa.output_phase[i];
    }
    std::cout << "]" << std::endl;
}

// Print a single XorRemain
inline void print_xor_remain(const talisman::XorRemain& xr) {
    if (verbose < 1) return;
    std::cout << "  XorRemain[" << xr.id << "]:" << std::endl;
    std::cout << "    nodes: [";
    for (size_t i = 0; i < xr.nodes.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << xr.nodes[i];
    }
    std::cout << "]" << std::endl;

    std::cout << "    tt: [";
    for (size_t i = 0; i < xr.tt.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << "0x" << std::hex << (int)xr.tt[i] << std::dec;
    }
    std::cout << "]" << std::endl;

    std::cout << "    input_phase: [";
    for (size_t i = 0; i < xr.input_phase.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << xr.input_phase[i];
    }
    std::cout << "]" << std::endl;

    std::cout << "    output_phase: [";
    for (size_t i = 0; i < xr.output_phase.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << xr.output_phase[i];
    }
    std::cout << "]" << std::endl;
}

// Print a single LinearCut
inline void print_linear_cut(const talisman::LinearCut& lc) {
    if (verbose < 1) return;
    std::cout << "  LinearCut[root=" << lc.root << "]:" << std::endl;
    std::cout << "    leaves: [";
    for (size_t i = 0; i < lc.leaves.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << lc.leaves[i];
    }
    std::cout << "]" << std::endl;

    std::cout << "    nodes: [";
    for (size_t i = 0; i < lc.nodes.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << lc.nodes[i];
    }
    std::cout << "]" << std::endl;

    std::cout << "    truth: 0x" << std::hex << lc.truth << std::dec << std::endl;

    // Print ANF polynomial: root = const ⊕ Σ(coeff_i * leaf_i)
    std::cout << "    ANF poly: n" << lc.root << " = ";
    bool first = true;
    if (lc.constant_coeff) {
        std::cout << "1";
        first = false;
    }
    for (size_t i = 0; i < lc.linear_coeffs.size() && i < lc.leaves.size(); i++) {
        if (lc.linear_coeffs[i]) {
            if (!first) std::cout << " ⊕ ";
            std::cout << "n" << lc.leaves[i];
            first = false;
        }
    }
    if (first) std::cout << "0";  // All coefficients are 0
    std::cout << std::endl;

    // Also print in linear algebra form: root - (const + Σ coeff_i * leaf_i) = 0
    std::cout << "    LinPoly: n" << lc.root;
    for (size_t i = 0; i < lc.linear_coeffs.size() && i < lc.leaves.size(); i++) {
        if (lc.linear_coeffs[i]) {
            std::cout << " - n" << lc.leaves[i];
        }
    }
    if (lc.constant_coeff) {
        std::cout << " - 1";
    }
    std::cout << " = 0" << std::endl;
}

// Print all LinearCuts
inline void print_linear_cuts(
    const std::vector<talisman::LinearCut>& cuts,
    int max_print = 10
) {
    if (verbose < 1) return;
    std::cout << "=== Linear Cuts Detection Results ===" << std::endl;
    std::cout << "Total Linear Cuts: " << cuts.size() << std::endl;

    for (size_t i = 0; i < cuts.size() && (int)i < max_print; i++) {
        print_linear_cut(cuts[i]);
    }
    if ((int)cuts.size() > max_print) {
        std::cout << "  ... (" << (cuts.size() - max_print) << " more)" << std::endl;
    }

    std::cout << "=====================================" << std::endl;
}

// Print all detection results
inline void print_detection_results(
    const std::vector<talisman::HalfAdder>& half_adders,
    const std::vector<talisman::FullAdder>& full_adders,
    const std::vector<talisman::XorRemain>& xor_list,
    int max_print = 10
) {
    if (verbose < 1) return;
    std::cout << "=== Detection Results ===" << std::endl;

    std::cout << "\nHalf Adders (" << half_adders.size() << " total):" << std::endl;
    for (size_t i = 0; i < half_adders.size() && (int)i < max_print; i++) {
        print_half_adder(half_adders[i]);
    }
    if ((int)half_adders.size() > max_print) {
        std::cout << "  ... (" << (half_adders.size() - max_print) << " more)" << std::endl;
    }

    std::cout << "\nFull Adders (" << full_adders.size() << " total):" << std::endl;
    for (size_t i = 0; i < full_adders.size() && (int)i < max_print; i++) {
        print_full_adder(full_adders[i]);
    }
    if ((int)full_adders.size() > max_print) {
        std::cout << "  ... (" << (full_adders.size() - max_print) << " more)" << std::endl;
    }

    std::cout << "\nXOR Remaining (" << xor_list.size() << " total):" << std::endl;
    for (size_t i = 0; i < xor_list.size() && (int)i < max_print; i++) {
        print_xor_remain(xor_list[i]);
    }
    if ((int)xor_list.size() > max_print) {
        std::cout << "  ... (" << (xor_list.size() - max_print) << " more)" << std::endl;
    }

    std::cout << "=========================" << std::endl;
}

// ======================================================================
// Legacy functions (kept for backward compatibility)
// ======================================================================

// Call Gia_detecting_ha_fa and print results (legacy version)
inline void detect_and_print_ha_fa(Gia_Man_t* pGia) {
    if (!pGia) {
        std::cerr << "Error: GIA is null" << std::endl;
        return;
    }

    if (verbose >= 1) std::cout << "=== Running HA/FA Detection ===" << std::endl;

    std::vector<talisman::HalfAdder> half_adders;
    std::vector<talisman::FullAdder> full_adders;
    std::vector<talisman::XorRemain> xor_list;

    if (detect_ha_fa_xor(pGia, half_adders, full_adders, xor_list)) {
        print_detection_results(half_adders, full_adders, xor_list);
    }

    if (verbose >= 1) std::cout << "=== HA/FA Detection Complete ===" << std::endl;
}

// Call Gia_PoLocalKCover_reverse and print results
inline void compute_and_print_kcover(Gia_Man_t* pGia, int k = 0) {
    if (!pGia) {
        std::cerr << "Error: GIA is null" << std::endl;
        return;
    }

    if (verbose >= 1) std::cout << "=== Running PO Local K-Cover Reverse (k=" << k << ") ===" << std::endl;

    Vec_Int_t* vBoundaryNodes = NULL;
    Vec_Int_t* vResult = Gia_PoLocalKCover_reverse(pGia, k, &vBoundaryNodes);

    if (!vResult) {
        std::cerr << "K-Cover computation returned null!" << std::endl;
        return;
    }

    int size = Vec_IntSize(vResult);
    if (verbose >= 1) std::cout << "Reverse FSA nodes count: " << size << std::endl;

    if (size > 0 && verbose >= 1) {
        std::cout << "Reverse FSA nodes: ";
        int entry, i;
        Vec_IntForEachEntry(vResult, entry, i) {
            std::cout << entry << " ";
            // Only print first 50 to avoid long output
            if (i >= 49) {
                std::cout << "... (truncated)";
                break;
            }
        }
        std::cout << std::endl;
    }

    // Print boundary nodes
    if (vBoundaryNodes) {
        int bsize = Vec_IntSize(vBoundaryNodes);
        if (verbose >= 1) std::cout << "Boundary nodes count: " << bsize << std::endl;

        if (bsize > 0 && verbose >= 1) {
            std::cout << "Boundary nodes: ";
            int entry, i;
            Vec_IntForEachEntry(vBoundaryNodes, entry, i) {
                std::cout << entry << " ";
                if (i >= 49) {
                    std::cout << "... (truncated)";
                    break;
                }
            }
            std::cout << std::endl;
        }
        Vec_IntFree(vBoundaryNodes);
    }

    // Free result
    Vec_IntFree(vResult);

    if (verbose >= 1) std::cout << "=== K-Cover Complete ===" << std::endl;
}

// Get K-Cover result as std::vector (returns both reverse_fsa_nodes and boundary_nodes)
// This version automatically uses cached vAdds from detect_ha_fa_xor if available
inline bool get_kcover_nodes(
    Gia_Man_t* pGia,
    int k,
    std::vector<int>& reverse_fsa_nodes,
    std::vector<int>& boundary_nodes
) {
    reverse_fsa_nodes.clear();
    boundary_nodes.clear();

    if (!pGia) {
        std::cerr << "Error: GIA is null" << std::endl;
        return false;
    }

    Vec_Int_t* vBoundaryNodes = NULL;
    Vec_Int_t* vResult = NULL;

    // Check if we have cached vAdds from detect_ha_fa_xor
    DetectionCache& cache = get_detection_cache();
    if (cache.vAdds) {
        // Use optimized version with cached vAdds
        if (verbose >= 1) std::cout << "Using cached vAdds for K-Cover computation" << std::endl;
        vResult = Gia_PoLocalKCover_reverse_with_adds(pGia, k, &vBoundaryNodes, cache.vAdds);
    } else {
        // Fall back to original version
        vResult = Gia_PoLocalKCover_reverse(pGia, k, &vBoundaryNodes);
    }

    if (!vResult) {
        std::cerr << "K-Cover computation returned null!" << std::endl;
        return false;
    }

    // Copy reverse FSA nodes
    int size = Vec_IntSize(vResult);
    reverse_fsa_nodes.reserve(size);
    int entry, i;
    Vec_IntForEachEntry(vResult, entry, i) {
        reverse_fsa_nodes.push_back(entry);
    }
    Vec_IntFree(vResult);

    // Copy boundary nodes
    if (vBoundaryNodes) {
        int bsize = Vec_IntSize(vBoundaryNodes);
        boundary_nodes.reserve(bsize);
        Vec_IntForEachEntry(vBoundaryNodes, entry, i) {
            boundary_nodes.push_back(entry);
        }
        Vec_IntFree(vBoundaryNodes);
    }

    if (verbose >= 1) std::cout << "K-Cover: reverse_fsa_nodes=" << reverse_fsa_nodes.size()
              << ", boundary_nodes=" << boundary_nodes.size() << std::endl;

    return true;
}

// Clear the detection cache (call this when done with the current AIG)
inline void clear_detection_cache() {
    get_detection_cache().clear();
}

// Print GIA basic statistics
inline void print_gia_stats(Gia_Man_t* pGia) {
    if (!pGia) {
        std::cerr << "Error: GIA is null" << std::endl;
        return;
    }

    if (verbose < 1) return;
    std::cout << "=== GIA Statistics ===" << std::endl;
    std::cout << "  Primary Inputs (PIs):  " << Gia_ManPiNum(pGia) << std::endl;
    std::cout << "  Primary Outputs (POs): " << Gia_ManPoNum(pGia) << std::endl;
    std::cout << "  AND gates:             " << Gia_ManAndNum(pGia) << std::endl;
    std::cout << "  Total objects:         " << Gia_ManObjNum(pGia) << std::endl;
    std::cout << "========================" << std::endl;
}

} // namespace abc_interface

#endif // ABC_SIMPLE_INTERFACE_HPP
