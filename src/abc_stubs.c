/**
 * @file abc_stubs.c
 * @brief Stub implementations for ABC functions that are disabled
 *
 * When linking ABC without kissat and eSLIM modules, we need to provide
 * stub implementations for functions that are still referenced in abc.c
 */

#include <stdio.h>

// Stub for kissat_solve_cnf (used in Abc_CommandAbc9Kissat)
int kissat_solve_cnf(void *solver, void *callback) {
    (void)solver;
    (void)callback;
    fprintf(stderr, "Warning: kissat_solve_cnf is not available (kissat disabled)\n");
    return -1;  // UNKNOWN
}

// Stubs for eSLIM (used in Abc_CommandAbc9eSLIM)
void seteSLIMParams(int a, int b, int c, int d, int e, int f) {
    (void)a; (void)b; (void)c; (void)d; (void)e; (void)f;
    fprintf(stderr, "Warning: seteSLIMParams is not available (eSLIM disabled)\n");
}

void* applyeSLIM(void *pGia, int a, int b, int c) {
    (void)pGia; (void)a; (void)b; (void)c;
    fprintf(stderr, "Warning: applyeSLIM is not available (eSLIM disabled)\n");
    return NULL;
}
