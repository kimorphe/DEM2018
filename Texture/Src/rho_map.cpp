#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

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
		void set_lim();
		void print();
		void fwrite(char *fname);
	private:
		double **mem_alloc(int nx,int ny);
};
void FIELD::init(int nx, int ny){
	dat=FIELD::mem_alloc(nx,ny);
	Xa[0]=0.0; Xa[1]=0.0;
	Xb[0]=1.0; Xb[1]=1.0;
	Ndat[0]=nx; Ndat[1]=ny;
};
void FIELD::set_lim(){
	int i;
	for(i=0;i<2;i++){
		Wd[i]=Xb[i]-Xa[i];
		dx[i]=Wd[i]/(Ndat[i]-1);
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
	fprintf(fp,"# Ndat[0], Ndat[1]\n");
	fprintf(fp,"%d %d\n",Ndat[0],Ndat[1]);
	fprintf(fp,"# dat[i][j]\n");
	for(i=0;i<Ndat[0];i++){
	for(j=0;j<Ndat[1];j++){
		fprintf(fp,"%lf\n",dat[i][j]);
	}
	}
	fclose(fp);

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

	while(x<Xa[0]) x+=Wd[0];
	while(y<Xa[1]) x+=Wd[1];
	while(x>Xb[0]) x-=Wd[0];
	while(y>Xb[1]) x-=Wd[1];

	ix=int((x-Xa[0])/dx[0]);
	iy=int((y-Xa[1])/dx[1]);

	int nwx=int(wx/dx[0]*0.5);
	int nwy=int(wy/dx[1]*0.5);

	int i,j,id,jd;
	double sum=0.0;
	int nblk=(nwx+1)*(nwy+1);
	double rx,ry;
	double arg,wgt,sig=wx/4.0;
	sig*=sig;
	int phs;
	for(i=ix-nwx; i<=ix+nwx; i++){
		id=i;
		if(id<0) id+=Ndat[0];
		if(id>=Ndat[0]) id-=Ndat[0];
		rx=(i-ix)*dx[0];
	for(j=iy-nwy; j<=iy+nwy; j++){
		ry=(j-iy)*dx[1];
		arg=(rx*rx+ry*ry)/sig;
		arg*=0.5;
		wgt=exp(-arg);

		jd=j;
		if(jd<0) jd+=Ndat[1];
		if(jd>=Ndat[1]) jd-=Ndat[1];
		phs=kcell[id][jd];
		if(phs >1){
		       	phs=1;
		}else{
			phs=0;
		}
		sum+=(phs*wgt);
	}
	}
	return(sum);
};
//--------------------------------------------------------


int main(){

	IMG im;
	FIELD fld;

	char fname[]="k240.dat";
	char fnout[]="k240s.dat";
	im.load(fname);
	//im.print();

	int ndat[2]={50,50};
	fld.init(ndat[0],ndat[1]);
	fld.Xa[0]=im.Xa[0];
	fld.Xa[1]=im.Xa[1];
	fld.Xb[0]=im.Xb[0];
	fld.Xb[1]=im.Xb[1];

	fld.set_lim();

	int i,j;
	double xx,yy;
	double wx=5.,wy=5.;
	for(i=0;i<ndat[0];i++){
		xx=fld.Xa[0]+fld.dx[0]*i;
		printf("i=%d/%d\n",i,ndat[0]);
	for(j=0;j<ndat[1];j++){
		yy=fld.Xa[1]+fld.dx[1]*j;
		fld.dat[i][j]=im.smooth(xx,yy,wx,wy);
	}
	}
	fld.fwrite(fnout);

	return(0);
};
