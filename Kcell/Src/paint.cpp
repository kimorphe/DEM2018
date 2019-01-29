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

	char fname[128];
	char fnout[128]; // output file 
	char fndem[128]; // "folder/dem.inp"
	char fnsht[128]; // "folder/sheet.dat"
	char fndat[128]; // "folder/x***.dat" (particle data file)
	FILE *fp,*fo;

	int npt,nst,ir0,ir1,ip1,ip2,*iphs;
	int i,j,k,l,N;
	int icavi, iline,Ndiv[2],nsig;
	int nsum,endcap;
	int irev[2],iprd[2];

	double xx,yy,sig0,sig,tt;
	double x1[2],x2[2],Wd[2],Xa[2],Xb[2],Xc[2];
	double tmp1,tmp2,exx,eyy,*alph;
	char cbff[128];
	SHEET *st;
	PRTCL *PT;

//	----------INPUT DATA -----------
	strcpy(fname,"mkinp.inp");
	fp=fopen(fname,"r");
	if(fp==NULL) show_msg(fname);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fndem);
	fscanf(fp,"%s\n",fnsht);
	fscanf(fp,"%s\n",fndat);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fnout);
	puts(fndat);

	sprintf(fnout,"kcell.dat");
	puts(fnout);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nsig);
		alph=(double *)malloc(sizeof(double)*nsig);
		iphs=(int *)malloc(sizeof(int)*nsig);
	for(i=0;i<nsig;i++){ 	
		fscanf(fp,"%lf %d\n",alph+i,iphs+i);
		printf("alph=%lf iphs=%d\n",alph[i],iphs[i]);
	}
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&endcap);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&icavi);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&iline);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",Ndiv,Ndiv+1);

	fclose(fp);


//	----------DEM PARAMETERS--------------
	sig0=1.5;
	fp=fopen(fndem,"r"); //	"dem.inp"
	//int nhead=17;
	int nhead=26;
	if(fp==NULL) show_msg(fname); 
	for(i=0;i<nhead;i++){
		fgets(cbff,128,fp);
	}
	fscanf(fp,"%d %d\n",iprd,iprd+1);

	fclose(fp);

//	----------DEM SHEET DATA --------------
//	strcpy(fname,"sheet.dat");
	fp=fopen(fnsht,"r");
	if(fp==NULL) show_msg(fname); 
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

//	------------------------
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
	printf("Xa=%lf %lf\n",Xa[0],Xa[1]);
	printf("Wd=%lf %lf\n",Wd[0],Wd[1]);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&npt);
	printf("npt=%d\n",npt);

		PT=(PRTCL *)malloc(sizeof(PRTCL)*npt);
		printf("npt=%d\n",npt);

	fgets(cbff,128,fp);
	double vx,vy,sigi,sigj;
	double sigs[2];
	for(i=0;i<npt;i++){
		//fscanf(fp,"%lf %lf %d %d\n",&xx,&yy,&ir0,&ir1);
		fscanf(fp,"%d %d %le %le %le %le %le %le\n",&ir0,&ir1,&xx,&yy,&vx,&vy,sigs,sigs+1);
		PT[i].setX(xx,yy);
		PT[i].irev[0]=ir0;
		PT[i].irev[1]=ir1;
		PT[i].sigs[0]=sigs[0];
		PT[i].sigs[1]=sigs[1];
	}

	/*
	fprintf(fo,"# time (ps) in DEM simulation\n");
	fprintf(fo,"%lf\n",tt);
	fprintf(fo,"# Xa[2], Xb[2]\n");
	fprintf(fo," %lf %lf\n",Xa[0],Xa[1]);
	fprintf(fo," %lf %lf\n",Xa[0]+Wd[0],Xa[1]+Wd[1]);
	fprintf(fo,"# Periodic Boundary (1:yes, 2:No) x,y-directions\n");
	fprintf(fo,"%d %d\n",iprd[0],iprd[1]);
	fprintf(fo,"# Ndiv[0], Ndiv[1]\n");
	fprintf(fo," %d %d\n",Ndiv[0],Ndiv[1]);
	*/
	Dom2D dom(Ndiv[0],Ndiv[1]);
	dom.Xa[0]=Xa[0]; dom.Xa[1]=Xa[1];
	dom.Xb[0]=Xa[0]+Wd[0];
	dom.Xb[1]=Xa[1]+Wd[1];
	dom.set_dx();


/*
	if(icavi==1){
		fprintf(fo,"##Circ\n");
		fprintf(fo,"%d \n",npt*nsig);
		for(l=0;l<nsig;l++){
			sig=sig0*alph[l];
			for(i=0;i<npt;i++) fprintf(fo,"%lf %lf %lf %d\n",PT[i].x[0],PT[i].x[1],sig,iphs[l]);
		}
		fprintf(fo,"##END\n");
	}

	double t0=0.9*1.0;
	if(iline==1){
		fprintf(fo,"##RECT\n");
		fprintf(fo,"%d\n",nsum*nsig);
		for(l=0;l<nsig;l++){
//			sig=sig0*alph[l];
		for(i=0;i<nst;i++){
		for(j=0;j<st[i].Np-1;j++){
			ip1=st[i].list[j];
			ip2=st[i].list[j+1];
			for(k=0;k<2;k++){
				x1[k]=PT[ip1].x[k];
				x2[k]=PT[ip2].x[k];
				if(iprd[k]==1){
					irev[k]=PT[ip2].irev[k]-PT[ip1].irev[k];
					x2[k]+=(irev[k]*Wd[k]);
				}	
			}
			sig=t0+(sigs[0]-t0)*alph[l];
			if(l==1) sig=t0;
			fprintf(fo,"%lf %lf  %lf %lf  %lf  %d %d\n",x1[0],x1[1],x2[0],x2[1],sig,endcap,iphs[l]);

		}
		}
		}
	fprintf(fo,"##END\n");
	}
	fclose(fp);
	fclose(fo);
*/

	Curve2D *crvs,crv,*hgts;
	crvs=(Curve2D *)malloc(sizeof(Curve2D)*nst);
	hgts=(Curve2D *)malloc(sizeof(Curve2D)*nst);
		for(i=0;i<nst;i++){
			crvs[i].init(st[i].Np);
			hgts[i].init(st[i].Np);
		for(j=0;j<st[i].Np;j++){
			ip1=st[i].list[j];
			//ip2=st[i].list[j+1];
			for(k=0;k<2;k++){
				x1[k]=PT[ip1].x[k]+PT[ip1].irev[k]*Wd[k];
				//x2[k]=PT[ip2].x[k]+PT[ip2].irev[k]*Wd[k];
			}
			//printf("%lf %lf\n",x1[0],x1[1]);
			crvs[i].x[j]=x1[0];
			crvs[i].y[j]=x1[1];
			hgts[i].x[j]=PT[ip1].sigs[0];
			hgts[i].y[j]=PT[ip1].sigs[1];
		}
			crvs[i].spline();
			hgts[i].spline();
		//puts("");
		}

	double ss,ds;
	int Ns;
	double dxds,dyds,sigp,sigm;
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

			tb.set(dxds,dyds);
			tt=tb.len();
			tb.div(tt);
			nb.set(-tb.x[1],tb.x[0]);
			//nb.times(sig);
			//printf("%lf %lf %lf\n",xx,yy,ss);
			sigp=hgts[i].intplx(ss)*0.5;
			sigm=hgts[i].intply(ss)*0.5;
			x1[0]=xx+nb.x[0]*sigp;
			x1[1]=yy+nb.x[1]*sigp;
			x2[0]=xx-nb.x[0]*sigm;
			x2[1]=yy-nb.x[1]*sigm;
			dom.draw_line(x1,x2,1,1);

			x1[0]=xx+nb.x[0]*0.45;
			x1[1]=yy+nb.x[1]*0.45;
			x2[0]=xx-nb.x[0]*0.45;
			x2[1]=yy-nb.x[1]*0.45;
			dom.draw_line(x1,x2,2,1);
		};
		//puts("");
	};
	for(i=0;i<nst;i++){
		//crv=crvs[i];
	for(j=0;j<st[i].Np-1;j++){
		x1[0]=crvs[i].x[j];
		x1[1]=crvs[i].y[j];
		x2[0]=crvs[i].x[j+1];
		x2[1]=crvs[i].y[j+1];	
		dom.draw_line(x1,x2,2,1);
	}
	}	

	dom.out_kcell(fnout);
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

