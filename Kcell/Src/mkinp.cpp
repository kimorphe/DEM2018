#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void show_msg(char *fname){
	printf("Can't find '%s'\n",fname);
	printf(" --> process terminated.");
	exit(-1);
}
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
	}

	fprintf(fo,"# time (ps) in DEM simulation\n");
	fprintf(fo,"%lf\n",tt);
	fprintf(fo,"# Xa[2], Xb[2]\n");
	fprintf(fo," %lf %lf\n",Xa[0],Xa[1]);
	fprintf(fo," %lf %lf\n",Xa[0]+Wd[0],Xa[1]+Wd[1]);
	fprintf(fo,"# Periodic Boundary (1:yes, 2:No) x,y-directions\n");
	fprintf(fo,"%d %d\n",iprd[0],iprd[1]);
	fprintf(fo,"# Ndiv[0], Ndiv[1]\n");
	fprintf(fo," %d %d\n",Ndiv[0],Ndiv[1]);

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

