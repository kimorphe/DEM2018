//---------------------------------------------------------------------------------
class Material{
	public:
		double gmm[3][3];	// interfacial energy (0:gas, 1:fluid, 2: solid)
		double thE;		// contact angle
		void load(double th); 	// set gmm values 
		double gmm_max;		// maximum value in gmm
		void  normalize();	// normalize gmm
		void  print_gmm(); // print interfacial energy
		double erg;	// cell energy
};
class Cell{
	public:
		QPatch *qp0;
		int cnct[8];//connected ?
		Cell *cncl[8]; //pointer to connected cells
		int cnct4[4];//connected ?
		Cell *cncl4[4]; //pointer to connected cells
		int nc;	// number  of connected cells
		Cell();	// constructuor
		bool bnd; // boundary cell (T/F)
		int ID;		// linear index for 2D Grid
		int iad;	// address in PoreCells[ncell]
		int phs;	// 0=gas, 1=fluid, 2=solid 
		int phs_bff;
		double erg;
		double erg_bff;
	private:
};
class cWalker{
	public:
		Cell *cl0;
		double x0,y0;
		double xn,yn;
		int ix,iy;
		int ix0,iy0;
		int ofx,ofy;
	private:
};
class PoreCells:public Tree4{
	public:
		int ncell; 	// number of cells
		Cell *cells;	// pointer(array) to cell class instances
		PoreCells();	// constructor
		void setup(Ellip *els,int nelp,bool set_opr, int Lev_Max,Bbox bx); // generate cells 
		void isetup(Ellip *els,int nelp,bool set_opr, int Lev_Max,int Lev_Exact,Bbox bx); // generate cells 
		int find(int id);	// find cell having a given linear grid number(id).
		void connect();	// establish neghboring cell connection 
		void connect4();	// establish neghboring cell connection 
		void l2ij(int l, int *i, int *j); // index transform ( linear to 2D index)
		double Sr;	// degree of water saturation
		int n_void,n_water; // number of gas and fluid cells, resp.
		int init(double sr); // initialize phase distribution
		int setup_phs_list(); // setup phase lists
		int *indx_w;	// index set of fluid phase
		int *indx_v;	// index set of gas phase
		double cell_energy(int iad); // evaluate interfacial energy/cell
		double cell_energy_diff(int iad); // evaluate interfacial energy/cell
		double total_energy(); // evaluate total interfacial energy
		Material mtrl;	// material constants (interfaceial energy)
		void load_gmm(double th); // load gmm data (th = fluid/solid contact angle)
		double Etot;	// total interfacial energy
		double swap(int id, int jd); // swap cell id & jd tempralily
		void reject_swap(int id, int jd); // apply swap 
		double MC_stepping(Temp_Hist TH, int *nswap, int seed);	// Monte Carlo stepping 
		void write_phs();
		int count_grids();
		int grid_type(int i,int j);
		int grid_type_verb(int i,int j);
		bool is_bnd_cell_grid(int i, int j);
		int grid_loc(int i,int j);
		void grid_connect(int i0, int j0, int cnct[4]);
		void fwrite_cells(char fn[128]);
		void load_cell_data(char fn[128]);
		void write_leaves();
		int ngap;		// inter-particle gap = 1:closed , 2:open
		int nwk;	// number of random walkers
		cWalker *wks;	// random walker array
		void setup_walkers(int n, int mseed);
		void rwk(int seed);
		double mean_u2();
		void mean_u(double *ux, double *uy);
		void write_wks(char fname[128],int istp);
		void rwk2(int seed);
		double mean_u2_new();
		void mean_u_new(double *ux, double *uy);
		void write_wks_new(char fname[128],int istp);
	private:
};
void copy_PoreCell_data(PoreCells *pc_From, PoreCells *pc_To);
void refine_PoreCell_data(PoreCells *pc_From, PoreCells *pc_To);
//---------------------------------------------------------------------------------
