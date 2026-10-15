#include "matio.h"

extern "C" {
mat_t *Mat_CreateVer(const char *, const char *, int) { return 0; }
matvar_t *Mat_VarCreate(const char *, int, int, int, size_t *, void *, int) { return 0; }
int Mat_VarWrite(mat_t *, matvar_t *, int) { return 0; }
void Mat_VarFree(matvar_t *) {}
int Mat_Close(mat_t *) { return 0; }
}
