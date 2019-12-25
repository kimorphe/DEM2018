#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "domain.h"
#include "spline.h"
//#include "fft.h"
#include <complex>
using namespace std;
/*
void show_msg(char *fname){
	printf("Can't find '%s'\n",fname);
	printf(" --> process terminated.");
	exit(-1);
}
*/
class SHEET{
	public:
		int Np;	// number of particles 
		double K1,K2;	// spring constants
		double r1,r2;	// natural lengths
		double sig_w,sig_s; // characteristic lengths
		void init(int N); // initialization
		int *list;	// particle numbers
		int wall;	// is_wall ? --> 1:yes, 0:no		
		double xa[2],xb[2],xc[2],rd; // bounding box;
	private:
};
class PRTCL{
	public:
		int mobile;	// 1 for mobile, 0 for immobile particle
		double x[2];	// position
		double v[2];	// velocity
		double F[2];	// force vector
		double m;	// mass
		int irev[2];	// cell index 
		void init();	// initialize class instance
		void setX(double x1, double x2); // set position vector
		double sigs[2];
		double n[2];	// normal vector;
	private:
 };

class DEM_DATA{
	public:
		SHEET *st;
		PRTCL *PT;
		Curve2D *crvs,*hgts;
		void load_dem_inp(char fname[128]);
		void load_sheet_data(char fname[128]);
		void load_ptc_data(char fname[128],bool init);
		void spline_fit(bool init);
		void paint(Dom2D &dom);
		int iprd[2];
		double Xa[2],Xb[2],Wd[2];
		double time;

		int nst; 
		int npt;
		void write_normal(char fname[128]);
	private:
};
class FIELD{
	public:
		int Ndiv[2];
		double Xa[2],Xb[2],Wd[2],dx[2];
		double **dat;
		void init(int nx, int ny);
		void set_dx();
		void set_lims(double xa[2], double xb[2]);
		void print();
		void fwrite(char *fname);
		void eval_rhox(PRTCL *PT, int npt);
	private:
		void mem_alloc();
};
void FIELD::eval_rhox(PRTCL *PT,int npt){
	int ipt,i,j, ix,iy;
	double xx,yy;
	int i1,i2,j1,j2;
	double sig3=4.0;
	double sig=sig3/3.0;
	double argx,argy,rx,ry,wgt;
	for(ipt=0; ipt<npt; ipt++){
		xx=PT[ipt].x[0];
		yy=PT[ipt].x[1];
		i1=int((xx-sig3-Xa[0])/dx[0]);
		i2=int((xx+sig3-Xa[0])/dx[0]);
		j1=int((yy-sig3-Xa[1])/dx[1]);
		j2=int((yy+sig3-Xa[1])/dx[1]);
		for(i=i1; i<=i2; i++){
			rx=Xa[0]+dx[0]*i-xx;
			argx=rx/sig;
			argx*=argx;
			ix=i;
			while(ix<0) ix+=Ndiv[0];
			while(ix>=Ndiv[0]) ix-=Ndiv[0];
		for(j=j1; j<=j2; j++){
			ry=Xa[1]+dx[1]*j-yy;
			argy=ry/sig;
			argy*=argy;
			wgt=exp(-0.5*(argx+argy));
			iy=j;
			while(iy<0) iy+=Ndiv[1];
			while(iy>=Ndiv[1]) iy-=Ndiv[1];
			dat[ix][iy]+=wgt;
		}
		}

		//ix=floor(xx-Xa[0])/dx[0];
		//iy=floor(yy-Xa[1])/dx[1];
		//dat[ix][iy]+=1.0;
	}
};
void FIELD::init(int nx,int ny){
	Ndiv[0]=nx;
	Ndiv[1]=ny;
	mem_alloc();
};
void FIELD::set_dx(){
	int i;
	for(i=0;i<2;i++){
		Wd[i]=Xb[i]-Xa[i];
		dx[i]=Wd[i]/Ndiv[i];
	}
};
void FIELD::fwrite(char *fname){
	FILE *fp=fopen(fname,"w");
	int i,j;
	fprintf(fp,"# time\n");
	fprintf(fp," 0.0\n");
	fprintf(fp,"# Xa[0], Xa[1]\n");
	fprintf(fp,"%lf %lf\n",Xa[0],Xa[1]);
	fprintf(fp,"# Xb[0], Xb[1]\n");
	fprintf(fp,"%lf %lf\n",Xb[0],Xb[1]);
	fprintf(fp,"# Ndiv[0], Ndiv[1]\n");
	fprintf(fp,"%d %d\n",Ndiv[0],Ndiv[1]);
	fprintf(fp,"# dat[i][j]\n");
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		fprintf(fp,"%lf\n",dat[i][j]);
	}
	}
	fclose(fp);

};
void FIELD::mem_alloc(){
	int ndat=Ndiv[0]*Ndiv[1];
	double *pt=(double *)malloc(sizeof(double)*ndat);
	dat=(double **)malloc(sizeof(double *)*Ndiv[0]);
	for(int i=0;i<Ndiv[0];i++) dat[i]=pt+i*Ndiv[1];
	for(int i=0;i<ndat;i++) pt[i]=0.0;
};
void FIELD::set_lims(double xa[2], double xb[2]){
	for( int i=0;i<2;i++){
		Xa[i]=xa[i];
		Xb[i]=xb[i];
	}
	FIELD::set_dx();
};

void SHEET:: init(int N){

	Np=N;
	list=(int *)malloc(sizeof(int)*Np);
};
void PRTCL::init(){
	x[0]=0.0; x[1]=0.0;
	v[0]=0.0; v[1]=0.0;
	irev[0]=0; irev[1]=0;
	F[0]=0.0; F[1]=0.0;
	m=1.0;	// mass
	mobile=1;	// set to "mobile" mode
}; 

void PRTCL::setX(double x1, double x2){
	x[0]=x1;
	x[1]=x2;
};


void DEM_DATA::load_dem_inp(char fname[128]){
	char cbff[128];
	int nhead=28; // number of header lines (to be skipped)

	FILE *fp=fopen(fname,"r"); // "dem.inp"
	if(fp==NULL) show_msg(fname); 

	for(int i=0;i<nhead;i++) fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",iprd,iprd+1); // periodic B.C. 
	fclose(fp);

	printf("iprd=%d %d\n",iprd[0],iprd[1]);
};
//	----------DEM SHEET DATA --------------
void DEM_DATA::load_sheet_data(char fname[128]){
	char cbff[128];
	FILE *fp=fopen(fname,"r"); // "sheet.dat"
	if(fp==NULL) show_msg(fname); 

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nst);
		st=(SHEET *)malloc(nst*sizeof(SHEET));
	int nsum=0;
	int i,j,N;
	for(i=0;i<nst;i++){
		fgets(cbff,128,fp);
		fscanf(fp,"%d\n",&N);
		nsum+=(N-1);
		st[i].init(N);
		for(j=0;j<N;j++) fscanf(fp,"%d",st[i].list+j);
		for(j=0;j<N;j++) st[i].list[j]--;
		fgets(cbff,2,fp);
	}
	fclose(fp);
};
void DEM_DATA::load_ptc_data(char fname[128],bool init){
	char cbff[128];
	FILE *fp=fopen(fname,"r"); // "x***.dat"
	if(fp==NULL) show_msg(fname); 

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&time);
	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	double exy,eyx;
	fscanf(fp,"%lf %lf %lf %lf\n",Wd,Wd+1,&exy,&eyx);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&npt);
	Xb[0]=Xa[0]+Wd[0];
	Xb[1]=Xa[1]+Wd[1];
	//;dom.set_dx();

	if(init) PT=(PRTCL *)malloc(sizeof(PRTCL)*npt);

	fgets(cbff,128,fp);
	int i,ir0,ir1;
	double xx,yy,vx,vy;
	double sigs[2];
	for(i=0;i<npt;i++){
		fscanf(fp,"%d %d %le %le %le %le %le %le\n",&ir0,&ir1,&xx,&yy,&vx,&vy,sigs,sigs+1);
		PT[i].setX(xx,yy);
		PT[i].irev[0]=ir0;
		PT[i].irev[1]=ir1;
		PT[i].sigs[0]=sigs[0];
		PT[i].sigs[1]=sigs[1];
	}
	fclose(fp);
};
void DEM_DATA::spline_fit(bool init){
	if(init){
		crvs=(Curve2D *)malloc(sizeof(Curve2D)*nst);
		hgts=(Curve2D *)malloc(sizeof(Curve2D)*nst);
	}
	int i,j,k,ip1;
	double x1[2];
	for(i=0;i<nst;i++){
		crvs[i].init(st[i].Np);
		hgts[i].init(st[i].Np);
	for(j=0;j<st[i].Np;j++){
		ip1=st[i].list[j];
		for(k=0;k<2;k++){
			x1[k]=PT[ip1].x[k]+PT[ip1].irev[k]*Wd[k];
		}
		//printf("%lf %lf\n",x1[0],x1[1]);
		crvs[i].x[j]=x1[0]; // parctcle x-cooridnate
		crvs[i].y[j]=x1[1]; // particle y-coordinate
		hgts[i].x[j]=PT[ip1].sigs[0]; // hydrated water thickness 
		hgts[i].y[j]=PT[ip1].sigs[1]; // 
	}
		crvs[i].spline(); // generate spline curves
		hgts[i].spline(); 
	}


	double dxds,dyds,ss,tt;
	Vec2 tb;
	for(i=0;i<nst;i++){
	for(j=0;j<st[i].Np;j++){
		ip1=st[i].list[j];
		ss=j;
		dxds=crvs[i].dxds(ss);
		dyds=crvs[i].dyds(ss);
		tb.set(dxds,dyds);	// tangential vector
		tt=tb.len();
		tb.div(tt);	// unit tangential vector
		PT[ip1].n[0]=-tb.x[1];
		PT[ip1].n[1]= tb.x[0];
	}
	}
};

void DEM_DATA::write_normal(char fname[128]){
	FILE *fp=fopen(fname,"w");
	int i;
	for(i=0;i<npt;i++){
		fprintf(fp,"%lf %lf\n",PT[i].n[0],PT[i].n[1]);
	};
	fclose(fp);
};
void DEM_DATA::paint(Dom2D &dom){
	int i,j;
	double x1[2],x2[2];
	dom.Xa[0]=Xa[0]; 
	dom.Xa[1]=Xa[1];
	dom.Xb[0]=Xb[0];
	dom.Xb[1]=Xb[1];
	dom.time=time;
	dom.set_dx();
	dom.iprd[0]=1;
	dom.iprd[1]=1;


	int Ns;
	double ds,ss,tt,xx,yy,dxds,dyds;
	Vec2 tb,nb;
	for(i=0;i<nst;i++){
		Ns=crvs[i].np*100;
		ds=double(crvs[i].np-1)/(Ns-1);
		for(j=0;j<Ns;j++){
			ss=ds*j;
			xx=crvs[i].intplx(ss); 
			yy=crvs[i].intply(ss);
			dxds=crvs[i].dxds(ss);
			dyds=crvs[i].dyds(ss);
			tb.set(dxds,dyds);	// tangential vector
			tt=tb.len();
			tb.div(tt);	// unit tangential vector
			nb.set(-tb.x[1],tb.x[0]);	// unit normal vector
/*
			sigp=hgts[i].intplx(ss)*0.5;
			sigm=hgts[i].intply(ss)*0.5;
			x1[0]=xx+nb.x[0]*sigp;
			x1[1]=yy+nb.x[1]*sigp;
			x2[0]=xx-nb.x[0]*sigm;
			x2[1]=yy-nb.x[1]*sigm;
			dom.draw_line(x1,x2,1,1); // paint pore water (fluid phase)
*/

			x1[0]=xx+nb.x[0]*0.40;
			x1[1]=yy+nb.x[1]*0.40;
			x2[0]=xx-nb.x[0]*0.40;
			x2[1]=yy-nb.x[1]*0.40;
			dom.draw_line(x1,x2,2,1); // paint clay sheet (solid phase)
		};
	};


	for(i=0;i<nst;i++){
	for(j=0;j<st[i].Np-1;j++){
		x1[0]=crvs[i].x[j];
		x1[1]=crvs[i].y[j];
		x2[0]=crvs[i].x[j+1];
		x2[1]=crvs[i].y[j+1];	
		dom.draw_line(x1,x2,2,1); // paint clay sheet (solid phase)
	}
	}	
};

double xrd_sum(PRTCL *PT, int npt){
	double pi=4.0*atan(1.0);
	int i,j,k,Nth=361;
	double kv[2];
	double lmb=1.5418e-01; // CuKa [nm]
	double kin=2.*pi/lmb;
	double th,dth=0.125*pi/(Nth-1);
	complex<double> zi=complex<double>(0.0,1.0);
	complex<double> Ith,dIth;
	int Nin=181;
	double phi,dphi=pi/(Nin-1);

	for(j=0;j<Nth;j++){
		th=dth*j;
		Ith=complex<double>(0.0,0.0);
	for(k=0;k<Nin;k++){
		phi=dphi*k;
		kv[0]=2.*kin*sin(th)*cos(0.5*pi+th+phi);
		kv[1]=2.*kin*sin(th)*sin(0.5*pi+th+phi);
		dIth=complex<double>(0.0,0.0);
		for(i=0;i<npt;i++){
			dIth+=exp(zi*(kv[0]*PT[i].x[0]+kv[1]*PT[i].x[1]));
		}
		Ith+=abs(dIth);
	}
		printf("%lf %lf\n",2.*th/pi*180.0,abs(Ith)/npt);
	}
};
double two_body_cor(PRTCL *PT, int npt, double *Wd, char *fname, bool iso){
	double pi=4.0*atan(1.0);
	int i,j,I,J;
	double rmax=Wd[0];
	if(rmax > Wd[1]) rmax=Wd[1];
	double dr=0.1;
	double r0=dr;	// radius of exclusion volume

	int Nr=int(rmax/dr);
	rmax=dr*Nr;
	double *hist;
	hist=(double *)malloc(sizeof(double)*Nr);
	for(i=0;i<Nr;i++) hist[i]=0.0;

	double x1[2],x2[2],xx,yy,rij;
	double n1[2],n2[2],beta;

	for(i=0; i<npt; i++){
		x1[0]=PT[i].x[0];
		x1[1]=PT[i].x[1];
		n1[0]=PT[i].n[0];
		n1[1]=PT[i].n[1];
	for(j=0; j<npt; j++){
		n2[0]=PT[j].n[0];
		n2[1]=PT[j].n[1];
		for(I=-1;I<=1;I++){
			x2[0]=PT[j].x[0]+I*Wd[0];
		for(J=-1;J<=1;J++){
			x2[1]=PT[j].x[1]+J*Wd[1];
			xx=x2[0]-x1[0];
			yy=x2[1]-x1[1];
			rij=sqrt(xx*xx+yy*yy);
			if(rij>=rmax) continue;
			if(rij < r0) continue;
			xx/=rij;
			yy/=rij;
			beta=(n1[0]*xx+n1[1]*yy)*(n2[0]*xx+n2[1]*yy);
			if(iso) beta=1.;	// particle alignment not considerd
			hist[int(rij/dr)]+=(beta*beta);
		}
		}
	}
	}

	FILE *fp=fopen(fname,"w");
	double ri;

	for(i=0;i<Nr;i++){
		ri=dr*(i+0.5);
	       	fprintf(fp,"%lf %lf\n",ri,hist[i]/(2.*pi*ri));
	};
	fclose(fp);

};

int main(int argc, char *argv[] ){

	double pi=atan(1.0)*4.0;
	char fname[128];
	char fnout1[128]; // output file (pixel image data)
	char fnout2[128]; // output file (2D FFT) 
	char fnout3[128]; // output file (XRD pattern)
	char fnout4[128]; // output file (Radial distribution function)
	char fnout5[128]; // output file (normal vector)
	char fnout6[128]; // output file (local number density)
	char fndem[128]; // "folder/dem.inp" (DEM main input)
	char fnsht[128]; // "folder/sheet.dat" (Clay sheet data)
	char fndat[128]; // "folder/x***.dat" (particle data file)
	FILE *fp;
	double Wd[2];
	int Ndiv[2];
	char cbff[128],dir[128],dir_out[128],head[128],tail[128];

//   ----------INPUT DATA (control data)-----------
	strcpy(fname,"xrd.inp");
	fp=fopen(fname,"r");
	if(fp==NULL) show_msg(fname);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fndem); // DEM data file
	fscanf(fp,"%s\n",fnsht); // sheet data file
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",dir);	// data directory
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",dir_out); // output data directory

	int nf1,nf2,nf_inc;
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d %d\n",&nf1,&nf2,&nf_inc); // File No.s to be read
	printf("%d %d %d\n",nf1,nf2,nf_inc);

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",head); // data file name prefix 
	fscanf(fp,"%s\n",tail); // data file name extension

	fgets(cbff,128,fp);
	int px,py;
	fscanf(fp,"%d %d\n",&px,&py); // Number of pixels
	Ndiv[0]=pow(2,px);
	Ndiv[1]=pow(2,py);
	printf("Ndiv=%d %d\n",Ndiv[0],Ndiv[1]);
	fclose(fp);

//		-------- MAIN ROUTINE --------
	DEM_DATA DM;		// DEM data class
	DM.load_dem_inp(fndem);	// load DEM parameters
	DM.load_sheet_data(fnsht); //load DEM sheet data 
	Dom2D dom(Ndiv[0],Ndiv[1]); // domain class to manipulate image data
	FIELD rho;
	rho.init(128,128);
	bool init=true,iso=false;
	for(int nf=nf1;nf<=nf2;nf+=nf_inc){
		//	FILE NAMES
		sprintf(fndat,"%s/x%d.dat",dir,nf);	// particle data file
		sprintf(fnout1,"%s%s%d.%s",dir_out,head,nf,tail); // image data file
		sprintf(fnout2,"%s%s%d.%s",dir_out,head,nf,"fft"); // FFT wave number spectrum 
		sprintf(fnout3,"%s%s%d.%s",dir_out,head,nf,"xrd"); // XRD Intensity plot 
		sprintf(fnout4, "%sx%d.%s",dir_out,nf,"rad"); // XRD Intensity plot 
		sprintf(fnout5, "%sx%d.%s",dir_out,nf,"nml"); // normal vector 
		sprintf(fnout6, "%sx%d.%s",dir_out,nf,"rho"); // local number density 
		printf("%s --> %s\n",fndat,fnout1); // Input,Output data files

		DM.load_ptc_data(fndat,init);	// load particle snapshot data
		DM.spline_fit(init);		// generate spline curves
		DM.write_normal(fnout5);
		init=false;

		rho.set_lims(DM.Xa,DM.Xb);
		rho.eval_rhox(DM.PT,DM.npt);
		rho.fwrite(fnout6);

		DM.paint(dom);	// convert particle to image data
		dom.FFT2D();	// Perform 2D FFT
		dom.out_kcell(fnout1); // write image to file
		dom.out_Kdat(fnout2); // write FFT data to file
		dom.XRD(fnout3);// synthesize XRD pattern by sampling FFT wave number spectrum 

		//xrd_sum(DM.PT,DM.npt);
		two_body_cor(DM.PT,DM.npt,DM.Wd,fnout4,iso); // radial distribution functioin
		
		dom.clear_kcell();	// clear image data
	}
	return(0);
}
