//double dist(double *x, double *y);
//void show_msg(char *fname);

class Vec2{
	public:
	double x[2];
	void set(double,double);	
	void set(double*);	
	Vec2 times(double);	
	Vec2 div(double);	
	double len();
	void print();
	private:
};
double iprod(Vec2 a,Vec2 b);
double iprod(double a[2], double b[2]);
Vec2 vsum(Vec2 a, Vec2 b);
Vec2 vdiff(Vec2 a, Vec2 b);

class MatTriDiag{	// triple diagonal matrix
	public:
		double *A0;	// A[i-1,i]
		double *A1;	// A[i,i]
		double *A2;	// A[i,i+1]
		double *y;	// solution vector (used when LinSolve(b) is called) 
		int n;		// matrix size (n x n)
		void init(int N); // initialize class instance (must be used firstly)
		void print(); // print matrix components
		// Linear Equation Solvers (Gaussian Elimination)
		// Note that the coeeficnet matrices are destroyed during the solution process
		void LinSolve(double *b); // take one r.h.s. vector --> solution =y
		void LinSolve(double *b1, double *b2); // take two r.h.s. vectors --> solutions b1,b2
	private:
	protected:
};

class Curve2D{	// plane (2D) curve 
	public:
		int np; // number of points
		double *x,*y; // coordinate vectors x,y
		double *bx,*cx,*dx; // cubic spline coefficients for x
		double *by,*cy,*dy; // cubic spline coefficients for y
		void init(int N); // initialize instance (must be called firstly)
		void print(); // print data 
		void print(char *fname); // print data 
		void spline(); // obtain spline coefficients
		double intplx(double s); // interpolate x(s) for s in (0, np-1)
		double intply(double s); // interpolate y(s) for s in (0, np-1)
		double dxds(double s); // interpolate dx/ds 
		double dyds(double s); // interpolate dy/ds
	private:
	protected:
};

