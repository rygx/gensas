// Minimal stand-in for <matio.h>, for scripts/local_solver_test.sh only.
// Declares just what sas/SasComputation.cpp uses to write MAT files; the stub implementations
// (matio_stub.cpp) do nothing, so writeMatFile() produces no file in a local test build.
#ifndef GENSAS_MATIO_STUB_H
#define GENSAS_MATIO_STUB_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _mat_t mat_t;
typedef struct _matvar_t matvar_t;

enum { MAT_FT_DEFAULT = 0x0100 };
enum { MAT_C_DOUBLE = 6 };
enum { MAT_T_DOUBLE = 9 };
enum { MAT_COMPRESSION_NONE = 0 };

mat_t *Mat_CreateVer(const char *matname, const char *hdr_str, int mat_file_ver);
matvar_t *Mat_VarCreate(const char *name, int class_type, int data_type, int rank, size_t *dims, void *data, int opt);
int Mat_VarWrite(mat_t *mat, matvar_t *matvar, int compress);
void Mat_VarFree(matvar_t *matvar);
int Mat_Close(mat_t *mat);

#ifdef __cplusplus
}
#endif

#endif
