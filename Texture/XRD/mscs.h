class Mscs{
	public :
		double x1[2], x2[2];   // centers of solid cylinder 
		double rd;		// radius of the cylinder
		double x0p[2], x0m[2]; // centers of meniscus perimeter
		double eh[2];	// unit vector 
		double nh[2];	// unit vector 
		double d12;
		double rho;	//radius of meniscus perimeter
		double alpha,theta;
		Mscs();
		void setup( double x1[2], double x2[2], double rd, double alpha, double theta);	
		int IsIn(double xf[2]);
		void show();
	private: 
};
