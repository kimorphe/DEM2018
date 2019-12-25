#include<stdio.h>
#include<stdlib.h>
#include<math.h> 
#include<random>


class Cell{
	public:
		int cnct[4];//connected ?
		Cell *cncl[4]; //pointer to connected cells
		int nc;	// number  of connected cells
		Cell();	// constructuor
		bool bnd; // boundary cell (T/F)
		int ID;		// linear index for 2D Grid
		int iad;	// address in PoreCells[ncell]
		int phs;	// 0=gas, 1=fluid, 2=solid 
		int phs_bff;
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
//-------------Dom2D Class ------------------
class Dom2D{
	public:
		double Xa[2],Xb[2],dx[2];// Computational Domain
		int Nx[2],Ndiv[2];
		int iprd[2];	// periodic B.C. (1:yes,2:No)
		int **kcell;
		Dom2D();
		Dom2D(int ndiv1, int ndiv2);
		Dom2D(char *fname);

		void out_kcell(char *fname);
		void clear_kcell();
		void set_dx();
		int count(int i);
		double time;	// time (ps) in DEM simulation
		void set_val(int val);
		void xy2ij(double x, double y, int indx[2]);

		int ncell;
		Cell *cl;
		void setup_cells();
		void connect_cells();
		int find_cell(int ID);

		int nwk;
		cWalker *wks;
		void init_walkers(int n);
	private:
		void mem_alloc();
};
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void show_msg(char *fname){
	printf("Can't find '%s'\n",fname);
	printf(" --> process terminated.");
	exit(-1);
}

Cell::Cell(){};


//-------------Dom2D Class ------------------
void Dom2D:: init_walkers(int n){

	std::mt19937_64 eng(-2);
	std::uniform_int_distribution<int>irnd(0,ncell-1);
	nwk=n;	// number of random walkers
	wks=(cWalker *)malloc(sizeof(cWalker)*nwk);

	int ID;
	for(int i=0;i<nwk;i++){
	       	wks[i].cl0=cl+irnd(eng);
		ID=wks[i].cl0->ID;
		printf("%d %d\n",int(ID/Ndiv[1]),ID%Ndiv[1]);
	}

};
void Dom2D:: setup_cells(){
	int i,j,k,id;
	ncell=0;
	for(i=0; i<Ndiv[0]; i++){
	for(j=0; j<Ndiv[0]; j++){
		if(kcell[i][j]==0) ncell++;
	}
	}

	cl=(Cell *)malloc(sizeof(Cell)*ncell);

	id=0;
	k=0;
	for(i=0; i<Ndiv[0]; i++){
	for(j=0; j<Ndiv[0]; j++){
		if(kcell[i][j]==0){
			cl[id].iad=id++;
			cl[id].ID=k;
			cl[id].phs=0;	// fluid phase 
		}	
		k++;
	}
	}
}
void Dom2D::connect_cells(){
	int i,j,k,ic;
	int id,jd;
	int iad;
	int ofsti[4]={ 0,1,0,-1};
	int ofstj[4]={-1,0,1, 0};
	Cell cll; 

	for(ic=0;ic<ncell;ic++){	// Cells
		cll=cl[ic];
		i=floor(cll.ID/Ndiv[1]);
		j=cll.ID%Ndiv[1];
		cl[ic].nc=0;
		for(k=0;k<4;k++){
			id=i+ofsti[k];
			jd=j+ofstj[k];
			cl[ic].cncl[k]=cl+ic;
			while(id <0) id+=Ndiv[0];
			while(id >= Ndiv[0]) id-=Ndiv[0];
			while(jd <0) jd+=Ndiv[1];
			while(jd >= Ndiv[1]) jd-=Ndiv[1];

			if(kcell[id][jd]!=0) continue;
			iad=find_cell(cll.ID);
			if(iad==-1) continue;
			cl[ic].cncl[k]=cl+iad;
			cl[ic].nc++;
		}
		
	};

};

int Dom2D::find_cell(int ID){

	int i1,i2,im;
	i1=0;
	i2=ncell-1;
	if(ID<cl[i1].ID) return(-1);
	if(ID>cl[i2].ID) return(-1);


	if(ID == cl[i1].ID) return(i1);
	if(ID == cl[i2].ID) return(i2);
	while(i2-i1>1){
		im=int(0.5*(i1+i2));
		if(ID ==cl[im].ID) return(im);
		if(ID < cl[im].ID){
			i2=im;
		}else{
			i1=im;
		}	
	};
	return(-1);
};
int Dom2D :: count(int iphs){
	int i,j,isum=0;
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		if(kcell[i][j]-1 == iphs) isum++;
	}
	}
	return(isum);
}
void Dom2D :: out_kcell(char *fname){

	FILE *fp=fopen(fname,"w");
	int i,j,kdat;

	fprintf(fp,"# time (ps) in DEM simulation\n");
	fprintf(fp,"%lf\n",time);

	fprintf(fp,"# Xa[0], Xa[1]\n");
	fprintf(fp," %lf %lf\n",Xa[0],Xa[1]);
	fprintf(fp,"# Xb[0], Xb[1]\n");
	fprintf(fp," %lf %lf\n",Xb[0],Xb[1]);
	fprintf(fp,"# Ndiv[0], Ndiv[1]\n");
	fprintf(fp,"%d %d\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"# kcell[i][j]\n");
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		kdat=kcell[i][j];
		fprintf(fp, "%d\n",kdat);
	}
	}
	fclose(fp);
};
void Dom2D::clear_kcell(){
	int i,j;
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		kcell[i][j]=0;
	}
	}
};
Dom2D::Dom2D(int Nx,int Ny){
	Ndiv[0]=Nx; 
	Ndiv[1]=Ny; 
	mem_alloc();
}
void Dom2D::set_dx(){
	int i,ndim=2;
	for(i=0;i<ndim;i++){
		dx[i]=(Xb[i]-Xa[i])/Ndiv[i];
	}
	iprd[0]=0;
	iprd[1]=0;
};
void Dom2D::set_val(int val){
	int i,j;
	for(i=0; i<Ndiv[0]; i++){
	for(j=0; j<Ndiv[1]; j++){
		kcell[i][j]=val;
	}
	}
};
void Dom2D::xy2ij(double x, double y, int indx[2]){
	int ix,iy; 
	ix=int((x-Xa[0])/dx[0]);
	iy=int((y-Xa[1])/dx[1]);
	while(ix < 0) ix+=Ndiv[0];
	while(iy < 0) iy+=Ndiv[1];
	while(ix >=Ndiv[0]) ix-=Ndiv[0];
	while(iy >=Ndiv[1]) iy-=Ndiv[1];
	indx[0]=ix;
	indx[1]=iy;
};
Dom2D::Dom2D(){		// Constructor 1
	int i,ndim=2;
	for(i=0;i<ndim;i++){
		Ndiv[i]=1; 
		Xa[i]=0.0; 
		Xb[i]=1.0; 
		dx[i]=(Xb[i]-Xa[i])/Ndiv[i];
		Nx[i]=Ndiv[i]+1;
	}
	mem_alloc();
};

Dom2D::Dom2D(char *fname){ //Contructor 2

	int i,j,ndim=2;
	FILE *fp;
	char cbff[128];
	fp=fopen(fname,"r");	
	if(fp==NULL) show_msg(fname);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&time);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Xa,Xa+1);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Xb,Xb+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",iprd,iprd+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",Ndiv,Ndiv+1);

	for(i=0;i<ndim;i++){
		dx[i]=(Xb[i]-Xa[i])/Ndiv[i];
		Nx[i]=Ndiv[i]+1;
	}
	mem_alloc();

	int itmp;
	fgets(cbff,128,fp);
	for(i=0; i<Ndiv[0]; i++){
	for(j=0; j<Ndiv[1]; j++){
		fscanf(fp,"%d\n",&itmp);
		kcell[i][j]=itmp;
	}
	}

	fclose(fp);

};
void Dom2D::mem_alloc(){
	int i,j,*ptmp;
	ptmp=(int *)malloc(sizeof(int)*Ndiv[0]*Ndiv[1]);
	kcell=(int **)malloc(sizeof(int*)*Ndiv[0]);

	for(i=0; i<Ndiv[0];i++){
		kcell[i]=(ptmp+(i*Ndiv[1]));
	}

	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		kcell[i][j]=0;
	}}
};
/*
int cod2indx(double x, double Xa,double dx){
	return(floor((x-Xa)/dx));
}
int cod2indx(double x, double Xa,double dx, int periodic){
	int indx=floor((x-Xa)/dx);
	if(periodic>0) indx=(indx%periodic);
};
double indx2cod(int indx,double Xa, double dx){
	return((indx+0.5)*dx+Xa);
};
*/

//----------------------------------------------------------
int main(){
	char fname[128]="k240.dat";
	FILE *fp;
	Dom2D dom(fname);
	dom.setup_cells();
	printf("ncell=%d\n",dom.ncell);
	dom.connect_cells();
	dom.init_walkers(100);

	return(0);
};

