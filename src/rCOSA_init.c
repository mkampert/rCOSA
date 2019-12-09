// RegisteringDynamic Symbols + tracing function fortran
#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
//#include <stdlib.h> // for NULL
// #include <R_ext/Rdynload.h


/* Fortran Calls Used for Tracing in the COSA algorithms */
void F77_SUB(printdouble)(double *a) {
    Rprintf("Double: %lf\n", *a);
}

void F77_SUB(printint)(int *k) {
    Rprintf("Int: %3d\n", *k);
}

void F77_SUB(print2ints)(int *k1, int *k2) {
    Rprintf("i1=%3d, i2=%3d\n", *k1, *k2);
}

void F77_SUB(print3ints)(int *k1, int *k2, int *k3) {
    Rprintf("i1=%3d, i2=%3d, i3=%3d\n", *k1, *k2, *k3);
}

void F77_SUB(checkpoint)() {
    Rprintf("Still running ... \n");
}

void F77_SUB(traceprint)(double *a1, double *a2, double *a3, double *a4) {
    Rprintf("d1=%lf, d2=%lf, d3=%lf, d4=%lf\n", *a1, *a2, *a3, *a4);
}


void F77_SUB(trialmessage)() {
    Rprintf("\n Trial period for this version of COSA has expired.\n A new version can be obtained from: https://github.com/mkampert/rCOSA \n When in R, run the following: devtools::install_github('mkampert/rCOSA') \n require(rCOSA)");
}

void F77_SUB(runningcosa)() {
    Rprintf("\n\n COSA executing (press ESC or ctrl+c to terminate). \n\n");
}

void F77_SUB(cosaheadlog)() {
    Rprintf("Wchange #iit #oit #it Eta CritMd CritMn\n");
}

void F77_SUB(cosa2017headlog)() {
    Rprintf("Wchange #iit #oit #it Eta CritW CritC CritH\n");
}

void F77_SUB(cosatracelog)(double *wtchg, int *iit, int *oit, int *it, double *eta, double *critw, double *critc, double *crith) {
    Rprintf("%1.5f %3d %3d %3d %1.3f %3.4f %3.4f %3.4f\n", *wtchg, *iit, *oit, *it, *eta, *critw, *critc, *crith);
}

void F77_SUB(oldcosatracelog)(double *wtchg, int *iit, int *oit, int *it, double *eta, double *critw, double *critc) {
  Rprintf("%1.5f %3d %3d %3d %1.3f %3.4f %3.4f\n", *wtchg, *iit, *oit, *it, *eta, *critw, *critc);
}

void F77_SUB(allcosatracelog)(double *wtchg, int *iit, int *oit, int *it, double *eta, double *critw, double *critc, double *crithw, double *crithc) {
  Rprintf("%1.5f %3d %3d %3d %1.3f %3.4f %3.4f %3.4f\n %3.4f\n", *wtchg, *iit, *oit, *it, *eta, *critw, *critc, *crithw, *crithc);
}



/* */
void R_init_rCOSAl(DllInfo *dll)
{
  // idea from https://github.com/swihart/rmutil/issues/1
  // dll, CEntries , ??,  Fortran, ??
  R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
