/* Test-only LD_PRELOAD shim: capture the actual MMG call prepared by DES.
 * Compile once per dimension. It leaves the production MMG adapter untouched.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef THREED
#include "mmg/mmg3d/libmmg3d.h"
#define ADAPT MMG3D_mmg3dlib
#define SAVE_MESH MMG3D_saveMesh
#define SAVE_SOL MMG3D_saveSol
#define ADAPT_NAME "MMG3D_mmg3dlib"
#define LINK_WRAPPER __wrap_MMG3D_mmg3dlib
#define LINK_ORIGINAL __real_MMG3D_mmg3dlib
#else
#include "mmg/mmg2d/libmmg2d.h"
#define ADAPT MMG2D_mmg2dlib
#define SAVE_MESH MMG2D_saveMesh
#define SAVE_SOL MMG2D_saveSol
#define ADAPT_NAME "MMG2D_mmg2dlib"
#define LINK_WRAPPER __wrap_MMG2D_mmg2dlib
#define LINK_ORIGINAL __real_MMG2D_mmg2dlib
#endif

#ifdef DES_CAPTURE_LINK_WRAP
extern int LINK_ORIGINAL(MMG5_pMesh, MMG5_pSol);
#undef ADAPT
#define ADAPT LINK_WRAPPER
#endif

int ADAPT(MMG5_pMesh mesh, MMG5_pSol metric) {
#ifdef DES_CAPTURE_LINK_WRAP
    int (*original)(MMG5_pMesh, MMG5_pSol) = LINK_ORIGINAL;
#else
    int (*original)(MMG5_pMesh, MMG5_pSol) = (int (*)(MMG5_pMesh, MMG5_pSol))dlsym(RTLD_NEXT, ADAPT_NAME);
#endif
    if (!original) { fprintf(stderr, "MMG capture: missing original adapter\n"); return MMG5_STRONGFAILURE; }
    const char *parent = getenv("DES_MMG_CAPTURE");
    if (!parent) return original(mesh, metric);
    char directory[4096], filename[4200];
    if (snprintf(directory, sizeof(directory), "%s/call-XXXXXX", parent) >= (int)sizeof(directory) || !mkdtemp(directory)) {
        fprintf(stderr, "MMG capture: cannot create capture directory\n"); return MMG5_STRONGFAILURE;
    }
    fprintf(stderr, "[mmg-capture] %s\n", directory);
    snprintf(filename, sizeof(filename), "%s/input.mesh", directory);
    if (!SAVE_MESH(mesh, filename)) return MMG5_STRONGFAILURE;
    snprintf(filename, sizeof(filename), "%s/input.sol", directory);
    if (!SAVE_SOL(mesh, metric, filename)) return MMG5_STRONGFAILURE;
    snprintf(filename, sizeof(filename), "%s/parameters.txt", directory);
    FILE *parameters = fopen(filename, "w");
    if (!parameters) return MMG5_STRONGFAILURE;
    fprintf(parameters, "-hmin %.17g\n-hmax %.17g\n-hausd %.17g\n-hgrad %.17g\n-hgradreq %.17g\n",
            mesh->info.hmin, mesh->info.hmax, mesh->info.hausd,
            mesh->info.hgrad < 0 ? -1 : exp(mesh->info.hgrad),
            mesh->info.hgradreq < 0 ? -1 : exp(mesh->info.hgradreq));
    if (mesh->info.nosizreq) fprintf(parameters, "-nosizreq\n");
    if (mesh->info.nosurf) fprintf(parameters, "-nosurf\n");
    if (fclose(parameters)) return MMG5_STRONGFAILURE;
    const int status = original(mesh, metric);
    if (status == MMG5_SUCCESS) {
        snprintf(filename, sizeof(filename), "%s/output.mesh", directory);
        if (!SAVE_MESH(mesh, filename)) return MMG5_STRONGFAILURE;
        snprintf(filename, sizeof(filename), "%s/output.sol", directory);
        if (!SAVE_SOL(mesh, metric, filename)) return MMG5_STRONGFAILURE;
    }
    return status;
}
