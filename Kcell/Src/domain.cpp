#include<iostream> 
#include<stdio.h>
#include<math.h> 
#include<string.h>
#include<stdlib.h>
#include "domain.h"
#include "mscs.h"
using namespace std;

void show_msg(char *fname){
	printf("Can't find '%s'\n",fname);
	printf(" --> process terminated.");
	exit(-1);
}
bool Circ :: isin(double *x){
	double dist;

	dist=(x[0]-xc[0])*(x[0]-xc[0]);
	dist+=(x[1]-xc[1])*(x[1]-xc[1]);
	dist=sqrt(dist);

	if(dist <= radi){
		return true;
	}{
		return false;
	}
};

//-------------Dom2D Class ------------------
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
//		fprintf(fp, "%d\n",kcell[i][j]-1);
		kdat=kcell[i][j];
//		if(kdat >0) kdat=1;
		fprintf(fp, "%d\n",kdat);
	}
	}
	fclose(fp);
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

	int i,ndim=2;
	FILE *fp;
	char cbff[128];
	fp=fopen(fname,"r");	
	if(fp==NULL) show_msg(fname);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&time);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf",Xa,Xa+1);
	fscanf(fp,"%lf %lf\n",Xb,Xb+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",iprd,iprd+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d %d\n",Ndiv,Ndiv+1);

	for(i=0;i<ndim;i++){
		dx[i]=(Xb[i]-Xa[i])/Ndiv[i];
		Nx[i]=Ndiv[i]+1;
	}

	fclose(fp);

	mem_alloc();
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

int Dom2D::perfo(char *fname){
	int i,j,k,ncirc=0,iphs;
	double xcod[2];
	bool io;
	Circ cdat,ctmp;
	int ix,iy,nrx,nry;
	int i1,i2,j1,j2;
	double Wd[2];


	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
	printf("Wd=%lf %lf\n",Wd[0],Wd[1]);

	FILE *fp;
	char cbff[128];
	fp=fopen(fname,"r");
	if(fp==NULL) puts("File open failed");
	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Circ")==0){
			fscanf(fp,"%d\n",&ncirc);
			for(k=0;k<ncirc;k++){
				fscanf(fp,"%lf %lf %lf %d\n",cdat.xc,cdat.xc+1, &cdat.radi,&iphs);
				perfo0(cdat,iphs);
			}
		}
	};
	fclose(fp);
	return(ncirc);

};
void Dom2D :: perfo0(Circ cdat, int iphs){

	int i,j,id,jd;
	int ix,iy,nrx,nry;
	int i1,i2,j1,j2;
	double xcod[2];
	bool io;

	ix=floor( (cdat.xc[0]-Xa[0])/dx[0] );
	iy=floor( (cdat.xc[1]-Xa[1])/dx[1] );
	nrx=ceil(cdat.radi/dx[0]);
	nry=ceil(cdat.radi/dx[1]);
	i1=ix-nrx;
	i2=ix+nrx;
	j1=iy-nry;
	j2=iy+nry;

	if(iprd[0]!=1){
		if (i1 < 0) i1=0;
		if (i2 >= Ndiv[0]) i2=Ndiv[0];
	}

	if(iprd[1]!=1){
		if (j1 < 0) j1=0;
		if (j2 >= Ndiv[1]) j2=Ndiv[1];
	}

	for(i=i1; i<i2; i++){
		xcod[0]=Xa[0]+dx[0]*(i+0.5);	
	for(j=j1; j<j2; j++){
		xcod[1]=Xa[1]+dx[1]*(j+0.5);	
		io=cdat.isin(xcod);	
		id=i%Ndiv[0];
		jd=j%Ndiv[1];
		if(id < 0 ) id+=Ndiv[0];
		if(jd < 0 ) jd+=Ndiv[1];
		if(io==true && kcell[id][jd]<iphs){
		 	kcell[id][jd]=iphs;
		}
	}
	}
	
};
int Dom2D::perfo_rec(char *fname){
	int i,j,k,nrec=0;
	double xcod[2],x1[2],x2[2],tt;
	int endcap,io,iphs;
	RECT rec;
	int i1,i2,j1,j2;
	int id,jd; 

	FILE *fp;
	char cbff[128];
	fp=fopen(fname,"r");
	if(fp==NULL) show_msg(fname);
	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##RECT")==0){
			fscanf(fp,"%d\n",&nrec);
			for(k=0;k<nrec;k++){
				fscanf(fp,"%lf %lf %lf %lf %lf %d %d\n",x1,x1+1,x2,x2+1,&tt,&endcap,&iphs);

				rec.setup(x1,x2,tt,endcap);

				i1=floor((rec.xmin-Xa[0])/dx[0]);
				i2=ceil((rec.xmax-Xa[0])/dx[0]);
				j1=floor((rec.ymin-Xa[1])/dx[1]);
				j2=ceil((rec.ymax-Xa[1])/dx[1]);

				if(iprd[0]!=1){	// non-periodic B.C. in x-direction
					if(i1< 0) i1=0;
					if(i2>= Ndiv[0]) i2=Ndiv[0];
				}
				if(iprd[1]!=1){	// non-periodic B.C. in y-direction
					if(j1< 0) j1=0;
					if(j2>= Ndiv[1]) j2=Ndiv[1];
				}

				for(i=i1; i<i2; i++){
					xcod[0]=Xa[0]+dx[0]*(i+0.5);	
				for(j=j1; j<j2; j++){
					xcod[1]=Xa[1]+dx[1]*(j+0.5);	
					io=rec.isin(xcod);	
						id=i%Ndiv[0];
						jd=j%Ndiv[1];
						if(id < 0) id+=Ndiv[0];
						if(jd < 0) jd+=Ndiv[1];
					if(io==1 && kcell[id][jd] < iphs){
						kcell[id][jd]=iphs;
					}
				}	// END_j
				}	// END_i
			}	// END_k
		}
	};	// END_while
	fclose(fp);
	return(nrec);
};
int Dom2D::mscs(char *fname){
	int i,j,k,ncirc=0;
	int ix,jy,iin,iphs;
	double xcod[2];
	bool io;
	Circ *cdat;
	Mscs mscs;
	
	double alp;	// meniscus aperture angle 
	double th;	// contact angle
	double rmax;	// maximum separation between cylinders

	FILE *fp;
	char cbff[128];

	fp=fopen(fname,"r");
	if(fp==NULL) puts("File open failed");
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Mscs")==0){
			fgets(cbff,128,fp);
			fscanf(fp,"%lf %lf %lf\n",&alp,&th,&rmax);
			printf("alp=%lf th=%lf rmax=%lf\n",alp,th,rmax);
		}
	};
	fclose(fp);

	fp=fopen(fname,"r");
	if(fp==NULL) puts("File open failed");
	i=0;
	while(fgets(cbff,7,fp) !=NULL){
		if(strcmp(cbff,"##Circ")==0){
			fscanf(fp,"%d\n",&ncirc);
			cdat=(Circ *)malloc(sizeof(Circ)*ncirc);	
			for(k=0;k<ncirc;k++){
				fscanf(fp,"%lf %lf %lf\n",cdat[k].xc,cdat[k].xc+1, &cdat[k].radi);
				printf("%lf %lf %lf\n",cdat[k].xc[0],cdat[k].xc[1], cdat[k].radi);
			}

			for(i=0;i<ncirc;i++){
			for(j=0;j<ncirc;j++){
				if(i == j) continue;
				if(dist2d(cdat[j].xc,cdat[i].xc)-2.*cdat[i].radi > rmax) continue; 
				mscs.setup(cdat[i].xc,cdat[j].xc,cdat[i].radi,alp,th);
				for(ix=0;ix<Ndiv[0];ix++){
					xcod[0]=Xa[0]+dx[0]*(ix+0.5);	
				for(jy=0;jy<Ndiv[1];jy++){
					xcod[1]=Xa[1]+dx[1]*(jy+0.5);	
					iin=mscs.IsIn(xcod);
					iphs=kcell[ix][jy];
					if(iin==0 && iphs ==0) kcell[ix][jy]=2;
					}
				}
				}
			}	
			}
	}
	return(ncirc);
};

int cod2indx(double x, double Xa,double dx, int imax){
	int indx=floor((x-Xa)/dx);
	if (indx<0) indx=0;
	if (indx> imax) indx=imax;
	return(indx);
};
double indx2cod(int indx,double Xa, double dx){
	return(indx*dx+Xa);
};
int Dom2D::draw_line(double x1[2], double x2[2],int iphs, int lw){

	double r12[2];
	int i1,i2,j1,j2,id,jd;
	double xmin,xmax;
	double ymin,ymax;
	double xcod,ycod;

	xmin=x1[0]; ymin=x1[1];
	if(xmin > x2[0]) xmin=x2[0];
	if(ymin > x2[1]) ymin=x2[1];

	xmax=x1[0]; ymax=x1[1];
	if(xmax < x2[0]) xmax=x2[0];
	if(ymax < x2[1]) ymax=x2[1];


	i1=cod2indx(x1[0],Xa[0],dx[0],Ndiv[0]-1);
	i2=cod2indx(x2[0],Xa[0],dx[0],Ndiv[0]-1);

	j1=cod2indx(x1[1],Xa[1],dx[1],Ndiv[1]-1);
	j2=cod2indx(x2[1],Xa[1],dx[1],Ndiv[1]-1);

	r12[0]=x2[0]-x1[0];
	r12[1]=x2[1]-x1[1];

	int ncell=0,itmp;
	if(fabs(r12[0]) > fabs(r12[1])){
		itmp=i1;
		if(i1>i2){
		       i1=i2;;
		       i2=itmp;
		}
		for(id=i1; id<=i2; id++){
			xcod=indx2cod(id,Xa[0],dx[0]);
			ycod=(xcod-x1[0])/r12[0]*r12[1]+x1[1];
			//if(ycod < ymin) continue;
			//if(ycod > ymax) continue;
			if(ycod > Xb[1]) continue;
			if(ycod < Xa[1]) continue;
			jd=cod2indx(ycod,Xa[1],dx[1],Ndiv[1]-1);
			kcell[id][jd]=iphs;
			for(int l=1;l<lw;l++){
				if(jd-l >0) kcell[id][jd-l]=iphs;
				if(jd+l <Ndiv[1]) kcell[id][jd+l]=iphs;
			}
			ncell++;
		};
	}else{
		itmp=j1;
		if(j1>j2){
		       j1=j2;;
		       j2=itmp;
		}
		for(jd=j1; jd<=j2; jd++){
			ycod=indx2cod(jd,Xa[1],dx[1]);
			xcod=(ycod-x1[1])/r12[1]*r12[0]+x1[0];
			//if(xcod < xmin) continue;
			//if(xcod > xmax) continue;
			if(xcod > Xb[0]) continue;
			if(xcod < Xa[0]) continue;
			id=cod2indx(xcod,Xa[0],dx[0],Ndiv[0]-1);
			kcell[id][jd]=iphs;
			for(int l=1;l<lw;l++){
				if(id-l >0) kcell[id-l][jd]=iphs;
				if(id+l <Ndiv[0]) kcell[id+l][jd]=iphs;
			}
			ncell++;
		}
	}

	return(ncell);
}


