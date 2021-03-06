#include <complex>
#include "fft.h"
using namespace std;
void mem_alloc2D(int nx, int ny, double ***pt);
void show_msg(char *fname);

//	********************************************
//
//		DOMAIN CLASS
//
//	********************************************
//-------------Circ Class--------------------
class Circ{
	public:
		double xc[2];	// center 
		double radi;	// radius
		bool isin(double* x);
	private:
};
//----------- Rectangle Class -------------
class RECT{
	public:
	double x[4],y[4];
	double xmin,xmax,ymin, ymax;
	void setup(
		double xc[2], //	center
		double b,  //	width
		double h,  //	height
		double th  //	counterclock-wise angle in degree
	);
	void setup(
		double x1[2], //endpoint 1
		double x2[2], //endpoint 2
		double t,	// thickness
		int endcap	// endcap (0:none, 1:add)
	);
	void out(char *fname);
	int isin(double xf[2]);
	private:
};

//-------------Dom2D Class ------------------
class Dom2D{
	public:
		double Xa[2],Xb[2],dx[2];// Computational Domain
		int Nx[2],Ndiv[2];
		int iprd[2];	// periodic B.C. (1:yes,2:No)
		int **kcell;
		int perfo(char *fn);
		int perfo_rec(char *fn);
		void perfo0(Circ cdat, int iphs);
		int mscs(char *fn);
		void out_kcell(char *fname);
		Dom2D();
		Dom2D(int ndiv1, int ndiv2);
		Dom2D(char *fname);
		void set_dx();
		int count(int i);
		double time;	// time (ps) in DEM simulation
		int draw_line(double x1[2], double x2[2],int iphs, int lw);
		void set_val(int val);
		void xy2ij(double x, double y, int indx[2]);
		void show_size();
		void clear_kcell();
		complex<double> **Kdat;
		DFT_prms dftx;
		DFT_prms dfty;
		void FFT2D();
		void out_Kdat(char *fname);
		void XRD(char *fnout);
	private:
		void mem_alloc();
		void mem_calloc();
};

double dist2d(double *x1, double *x2);

