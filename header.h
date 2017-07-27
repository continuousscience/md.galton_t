#include <sil_ext.h>
#include <proto_prim.h>
#include <math.h>

#define IP2(x,y) ((x)[0]*(y)[0] + (x)[1]*(y)[1])

// Data structure holding problem setup for periodic 2D Lorentz
// gas.  Input and output data are always in the 'orig' frame,
// but internally all coordinates (L, iL, x, v, and sph_crds)
// are in the 'rotated' frame.  This is so that g is always along
// +x in the rotated frame.
typedef struct {
    double L[2][2];
    double iL[2][2];
    //double x[2], v[3];
    double frame[2][2]; // x(orig) = frame . x(rotated)
    double g; // accel. constant
    double gamma; // damping constant - should be bM_2 * sigma**2
    double sigma; // noise constant - should be sqrt(gam/bM_2), gam << bM_2
    int co_linear, nsph;
    double sph_crds[0]; // NSPH x 3 (x,y,R)
} Galton;
// number of doubles in Galton
#define GALTON_DOUBLES (12+3)

static const unsigned char galton_hash[HASH_SIZE+1] = /*!hash!*/;

// Handy 2x2 matrix multiplies.
static void d2mv(double y[2], const double A[2][2], const double x[2]) {
    y[0] = A[0][0]*x[0] + A[0][1]*x[1];
    y[1] = A[1][0]*x[0] + A[1][1]*x[1];
}

static void d2mvT(double y[2], const double A[2][2], const double x[2]) {
    y[0] = A[0][0]*x[0] + A[1][0]*x[1];
    y[1] = A[0][1]*x[0] + A[1][1]*x[1];
}

