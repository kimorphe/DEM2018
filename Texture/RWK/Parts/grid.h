class Node{
	public:
		int iad;// data count (address in the Node class array NDs[]) 
		int id;	// grid number in the underlying regular grid
		int cnct[4];	// connected  nodes
		Node *cnds[4];	// connected nodes 
		Node();	//default constructor
		int nc;
		double x0,y0;
		bool sld;	// solid phase grid (true/false)
	private:
};
class Walker{	// Random Walker
	public:
		Node *nd0;
		double x0,y0;
		double xn,yn;
		int ofx,ofy;
	private:
};
class Grid{
	public:
		int Ng;		// total number of grids
		int Nx,Ny;	// number of underling regular grid (Ng != Nx x Ny)
		Node *NDs;	// grid points
		bool ready;
		double Xa[2],Wd[2],Xb[2],dx[2];
		int Ndiv[2];
		Grid();	// default constructor
		Grid(int n);
		void set_grid_params(double xa[2], double xb[2], int nx, int ny);
		void setup(Tree4 tr4);
		int find(int  k);
		void l2cod(int l,double *x, double *y);
		void indx2cod(int i,int j, double *x, double *y);
		void l2ij(int l, int *i, int *j);
		void grid_cod(int inod,double *xcod, double *ycod);
		void connect();
		int nwk;	// number of random walkers
		Walker *wks;	// random walkers
		void setup_walkers(int n, int sd);
		std::mt19937_64 engine;
		std::uniform_int_distribution<> irnd3;
		void init_rand(int seed);
		void rwk();
		double tolx,toly;
		void write_wks(char fname[128],int istp);
		double mean_u2();
		void mean_u(double *Ux, double *Uy);
		void dbg_connectivity();	
	private:
		void mem_alloc();
};

