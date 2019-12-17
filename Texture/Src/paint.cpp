#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "domain.h"
#include "spline.h"

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
	private:
 };


int main(int argc, char *argv[] ){

	double pi=atan(1.0)*4.0;
	char fname[128];
	char fnout[128]; // output file 
	char fnout2[128]; // output file 
	char fndem[128]; // "folder/dem.inp" (DEM main input)
	char fnsht[128]; // "folder/sheet.dat" (Clay sheet data)
	char fndat[128]; // "folder/x***.dat" (particle data file)
	FILE *fp,*fo;

	int npt,nst,ir0,ir1,ip1,ip2;
	int nsum; 
	int i,j,k,l,N;
	int Ndiv[2],irev[2],iprd[2];
		//,nsig;

	double xx,yy,sig0,sig,tt;
	double x1[2],x2[2],Wd[2],Xa[2],Xb[2],Xc[2];
	char cbff[128],dir[128],dir_out[128],head[128],tail[128];
	SHEET *st;
	PRTCL *PT;

//   ----------INPUT DATA (control data)-----------
	strcpy(fname,"paint.inp");
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
	fscanf(fp,"%d %d\n",Ndiv,Ndiv+1); // Number of pixels
	fclose(fp);

//	----------DEM PARAMETERS--------------
	sig0=1.5;
	fp=fopen(fndem,"r"); //	"dem.inp"
	int nhead=26; // number of header lines (to be skipped)
	if(fp==NULL) show_msg(fndem); 
	for(i=0;i<nhead;i++){
		fgets(cbff,128,fp);
	}
	fscanf(fp,"%d %d\n",iprd,iprd+1); // periodic B.C. 

	fclose(fp);

//	----------DEM SHEET DATA --------------
	fp=fopen(fnsht,"r");
	if(fp==NULL) show_msg(fnsht); 
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nst);
		st=(SHEET *)malloc(nst*sizeof(SHEET));
	nsum=0;
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

	
	for(int nf=nf1;nf<=nf2;nf+=nf_inc){
		sprintf(fndat,"%s/x%d.dat",dir,nf);
		sprintf(fnout,"%s%s%d.%s",dir_out,head,nf,tail);
		sprintf(fnout2,"%s%sn%d.%s",dir_out,head,nf,tail);
		printf("%s --> %s\n",fndat,fnout); // Input,Output data files

//	-----------LOAD PARTICLE MOTION DATA -------------
	fp=fopen(fndat,"r");
	if(fp==NULL) show_msg(fndat); 

	fo=fopen(fnout,"w");
	if(fo==NULL) show_msg(fnout);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&tt);

	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	double exy,eyx;
	fscanf(fp,"%lf %lf %lf %lf\n",Wd,Wd+1,&exy,&eyx);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&npt);

	PT=(PRTCL *)malloc(sizeof(PRTCL)*npt);

	fgets(cbff,128,fp);
	double vx,vy;
	double sigs[2];
	for(i=0;i<npt;i++){
		fscanf(fp,"%d %d %le %le %le %le %le %le\n",&ir0,&ir1,&xx,&yy,&vx,&vy,sigs,sigs+1);
		PT[i].setX(xx,yy);
		PT[i].irev[0]=ir0;
		PT[i].irev[1]=ir1;
		PT[i].sigs[0]=sigs[0];
		PT[i].sigs[1]=sigs[1];
	}

	Dom2D dom(Ndiv[0],Ndiv[1]);
	dom.Xa[0]=Xa[0]; 
	dom.Xa[1]=Xa[1];
	dom.Xb[0]=Xa[0]+Wd[0];
	dom.Xb[1]=Xa[1]+Wd[1];
	dom.time=tt;
	dom.set_dx();

	Dom2D Th(Ndiv[0],Ndiv[1]);
	Th.Xa[0]=Xa[0]; 
	Th.Xa[1]=Xa[1];
	Th.Xb[0]=Xa[0]+Wd[0];
	Th.Xb[1]=Xa[1]+Wd[1];
	Th.time=tt;
	Th.set_dx();
	Th.set_val(-90);



	Curve2D *crvs,*hgts;
	crvs=(Curve2D *)malloc(sizeof(Curve2D)*nst);
	hgts=(Curve2D *)malloc(sizeof(Curve2D)*nst);
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

	double ss,ds;
	int Ns;
	double dxds,dyds,sigp,sigm;
	Vec2 tb,nb;
	double th_n; 
	int ix,iy,indx[2];
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
			//nb.times(sig);
			//printf("%lf %lf %lf\n",xx,yy,ss);
			sigp=hgts[i].intplx(ss)*0.5;
			sigm=hgts[i].intply(ss)*0.5;
			x1[0]=xx+nb.x[0]*sigp;
			x1[1]=yy+nb.x[1]*sigp;
			x2[0]=xx-nb.x[0]*sigm;
			x2[1]=yy-nb.x[1]*sigm;
			dom.draw_line(x1,x2,1,1); // paint pore water (fluid phase)

			x1[0]=xx+nb.x[0]*0.45;
			x1[1]=yy+nb.x[1]*0.45;
			x2[0]=xx-nb.x[0]*0.45;
			x2[1]=yy-nb.x[1]*0.45;
			dom.draw_line(x1,x2,2,1); // paint clay sheet (solid phase)

			if(nb.x[0] > 1.0) nb.x[0]=1.0;
			if(nb.x[0] <-1.0) nb.x[0]=-1.0;
			th_n=acos(nb.x[0])/pi*180.0;
			if(nb.x[1]<0.0) th_n=360.0-th_n;
			Th.xy2ij(xx,yy,indx);
			/*
			ix=int((xx-Xa[0])/Th.dx[0]);
			iy=int((yy-Xa[1])/Th.dx[1]);
			if(ix < 0) ix+=Ndiv[0];
			if(iy < 0) iy+=Ndiv[1];
			if(ix >=Ndiv[0]) ix-=Ndiv[0];
			if(iy >=Ndiv[1]) iy-=Ndiv[1];
			*/
			Th.kcell[indx[0]][indx[1]]=(int(th_n)%180);
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
	dom.out_kcell(fnout);
	Th.out_kcell(fnout2);

	free(crvs);
	free(hgts);

	}
	return(0);
}

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

