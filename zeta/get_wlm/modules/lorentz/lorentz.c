#include <math.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <lorentz.h>

#ifndef M_PI 
#define M_PI 3.141592653589793238462643383279502884197169
#endif

void lorentz_gamma_inv(double const *const a, double const *const v, double *g)
{
	memcpy(g,a,3*sizeof(double));

	double vv = pow(v[0],2)+pow(v[1],2)+pow(v[2],2);

	/*return g=a if |v|=0 within precision*/
	if (vv <= DBL_EPSILON){
		return;
	}

	double ginvfac=sqrt(1-vv);
	double av=a[0]*v[0] + a[1]*v[1] + a[2]*v[2];
	for(int i=0; i<3; i++){
		g[i]=(ginvfac - 1)*v[i]*av/vv + a[i];
	}
}

void lorentz_gamma(double const *const a, double const *const v, double *g)
{
	memcpy(g,a,3*sizeof(double));

	double vv = pow(v[0],2)+pow(v[1],2)+pow(v[2],2);

	/*return g=a if |v|=0 within precision*/
	if (vv <= DBL_EPSILON){
		return;
	}

	double gfac=1/sqrt(1-vv);
	double av=a[0]*v[0] + a[1]*v[1] + a[2]*v[2];
	for(int i=0; i<3; i++){
		g[i]=(gfac - 1)*v[i]*av/vv + a[i];
	}
}

double lorentz_gamma_scalar(double const E, int const *const d, int const Nx)
{
	double twopiL = 2.0*M_PI/((double) Nx);
	double d2 = (pow(d[0],2)+pow(d[1],2)+pow(d[2],2))*pow(twopiL,2);
	double E2 = pow(E,2);
    double gamma = E / sqrt(E2 - d2);
    return gamma;
}
