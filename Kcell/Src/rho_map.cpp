#include<stdio.h>
#include<stdlib.h>
#include<string.h>

class IMG{
	public:
		int Ndat[2];
		double Xa[2],Xb[2],Wd[2];
		int **kcell;
		double dx[2];
		char fname[128];
		void load(char fn[128]);
		void print();
		double smooth(double x, double y, double wx, double wy);
	private:
		int **mem_alloc(int nx,int ny);
};
class FIELD{
	public:
		int Ndat[2];
		double Xa[2],Xb[2],Wd[2];
		double **dat;
		double dx[2];
		//char fname[128];
		///void load(char fn[128]);
		void init(int nx,int ny);
		void print();
	private:
		double **mem_alloc(int nx,int ny);
};
void FIELD::init(int nx, int ny){
	FIELD::mem_alloc(nx,ny);
	Xa[0]=0.0; Xa[1]=0.0;
	Xb[0]=1.0; Xb[1]=1.0;
	Ndat[0]=nx; Ndat[1]=ny;
};
double **FIELD::mem_alloc(int nx, int ny){
	int ndat=nx*ny;
	double *pt=(double *)malloc(sizeof(double)*ndat);
	double **pt2=(double **)malloc(sizeof(double *)*nx);
	for(int i=0;i<nx;i++) pt2[i]=pt+i*ny;
	for(int i=0;i<ndat;i++) pt[i]=0.0;
	return(pt2);
};


int **IMG::mem_alloc(int nx, int ny){
	int ndat=nx*ny;
	int *pt=(int *)malloc(sizeof(int)*ndat);
	int **pt2=(int **)malloc(sizeof(int *)*nx);
	for(int i=0;i<nx;i++) pt2[i]=pt+i*ny;
	for(int i=0;i<ndat;i++) pt[i]=0;
	return(pt2);
};
void IMG::load(char fn[128]){
	strcpy(fname,fn);

	FILE *fp=fopen(fname,"r");
	char cbff[128];	

	fgets(cbff,128,fp);
	double time;
	fscanf(fp,"%lf\n",&time);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&time);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Xb,Xb+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",Ndat,Ndat+1);
	fgets(cbff,128,fp);
	puts(cbff);

	int i;
	for( i=0;i<2;i++){
		Wd[i]=Xb[i]-Xa[i];
		dx[i]=Wd[i]/(Ndat[i]-1);
	};

	//int *pt=(int *)malloc(sizeof(int)*Ndat[0]*Ndat[1]);
	//kcell=(int **)malloc(sizeof(int *)*Ndat[0]);
	//for(i=0;i<Ndat[0];i++) kcell[i]=pt+i*Ndat[1]; 
	kcell=IMG::mem_alloc(Ndat[0],Ndat[1]);
	for(i=0;i<Ndat[0]*Ndat[1];i++) fscanf(fp,"%d\n",kcell[0]+i);
	fclose(fp);
};
void IMG::print(){
	int i,j,inc=10;
	for(i=0;i<Ndat[0];i+=inc){
	for(j=0;j<Ndat[1];j+=inc){
		printf("%d ",kcell[i][j]);
	}
	printf("\n");
	}
};

double IMG::smooth(double x, double y, double wx, double wy){
	double rho=0.0;
	
	int ix,iy;
	ix=(x-Xa[0])/dx[0];
	iy=(y-Xa[1])/dx[1];


	return(rho);
};
//--------------------------------------------------------


int main(){

	IMG im;
	FIELD fld;

	char fname[]="k240.dat";
	im.load(fname);
	im.print();

	fld.init(100,100);
	return(0);
};
