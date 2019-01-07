#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dem.h"
#include "MT.h"

void SHEET:: init(int N_particle){

	Np=N_particle;
	list=(int *)malloc(sizeof(int)*Np);
};
void SHEET :: set_prms(
	CNTRL prms,	// control parameters
	PRTCL *PT 	// particle class data
){
	int  i1,i2,i3,imb;

	K1=prms.K1;
	K2=prms.K2;

	r1=0.0; r2=0.0;
	if(Np >1){
		i1=list[0];
		i2=list[1];	
		 r1=dist(PT[i1].x, PT[i2].x);
	}
	if(Np >2){
		i3=list[2];
		r2=dist(PT[i1].x, PT[i3].x);
	}

	sig_w=prms.sig;
	sig_s=pow(2,0.16666667)*r1;
	wall=0;
	if(Np == 1){
		imb=PT[list[0]].mobile;
		switch(imb){
			case  1:  // regular particle
				wall=0;
				break;
			case -1: // horizontal wall (bottom)
			case -3: // horizontal wall (top)
				wall=1;
				break;
			case -2: // vertical wall
			case -4: // vertical wall
				wall=2;
				break;
		}
	}
};


void SHEET :: set_prms(
	double k1,	// spring const 1 
	double k2, 	// spring const 2
	PRTCL *PT, 	// particle class data
	double sig	// characteristic length
){

	int  i1,i2,i3,imb;

	K1=k1; // spring constant 1  
	K2=k2; // spring constant 2

	r1=0.0; r2=0.0;
	if(Np >1){
		i1=list[0];
		i2=list[1];	
		 r1=dist(PT[i1].x, PT[i2].x);
	}
	if(Np >2){
		i3=list[2];
		r2=dist(PT[i1].x, PT[i3].x);
	}

	sig_w=sig;
	sig_s=pow(2,0.16666667)*r1;

	wall=0;
	if(Np == 1){
		imb=PT[list[0]].mobile;
		switch(imb){
			case  1:  // regular particle
				wall=0;
				break;
			case -1: // horizontal wall (bottom)
			case -3: // horizontal wall (top)
				wall=1;
				break;
			case -2: // vertical wall
			case -4: // vertical wall
				wall=2;
				break;
		}
	}	
};


void SHEET :: bbox(PRTCL *PT, CNTRL prms){
	int i;
	double xx,yy;
	
	xa[0]=PT[list[0]].x[0]+PT[list[0]].irev[0]*prms.Wd[0];
	xa[1]=PT[list[0]].x[1]+PT[list[0]].irev[1]*prms.Wd[1];
	xb[0]=xa[0];
	xb[1]=xa[1];

	x0[0]=0.0;
	x0[1]=0.0;
	for(i=1;i<Np;i++){
		xx=PT[list[i]].x[0]+PT[list[i]].irev[0]*prms.Wd[0];
		yy=PT[list[i]].x[1]+PT[list[i]].irev[1]*prms.Wd[1];
		if(xx < xa[0]) xa[0]=xx;
		if(yy < xa[1]) xa[1]=yy;
		if(xx > xb[0]) xb[0]=xx;
		if(yy > xb[1]) xb[1]=yy;

		x0[0]+=xx;
		x0[1]+=yy;
	}

	x0[0]/=Np;
	x0[1]/=Np;

	xc[0]=.5*(xa[0]+xb[0]);
	xc[1]=.5*(xa[1]+xb[1]);

	rd =(xb[0]-xc[0])*(xb[0]-xc[0]);
	rd+=(xb[1]-xc[1])*(xb[1]-xc[1]);
	rd=sqrt(rd);
}

double SHEET ::incN(PRTCL *PTC, REV rev, double Sab[2][2]){

	int jp,jpt1,jpt2; 
	double dFn[2],UE=0.0;
	double dSab[2][2];

	Sab[0][0]=0.0;
	Sab[1][0]=0.0;
	Sab[1][1]=0.0;
	for(jp=0;jp<Np-1;jp++){
		jpt1=list[jp];
		jpt2=list[jp+1];
		UE+=STF1(PTC[jpt1],PTC[jpt2],dFn,r1,K1,rev,dSab);
		PTC[jpt1].F[0]+=dFn[0];
		PTC[jpt1].F[1]+=dFn[1];
		PTC[jpt2].F[0]-=dFn[0];
		PTC[jpt2].F[1]-=dFn[1];

		Sab[0][0]+=dSab[0][0];
		Sab[1][0]+=dSab[1][0];
		Sab[1][1]+=dSab[1][1];
	}
	return(UE);
}
double SHEET :: incQ(PRTCL *PTC, REV rev, double Sab[2][2]){
	int jp,jpt1,jpt2;
	double dFn[2],UE=0.0;
	double dSab[2][2];

	Sab[0][0]=0.0;
	Sab[1][0]=0.0;
	Sab[1][1]=0.0;
	for(jp=0;jp<Np-2;jp++){
		jpt1=list[jp];
		jpt2=list[jp+2];
		UE+=STF1(PTC[jpt1],PTC[jpt2],dFn,r2,K2,rev,dSab);
		PTC[jpt1].F[0]+=dFn[0];
		PTC[jpt1].F[1]+=dFn[1];
		PTC[jpt2].F[0]-=dFn[0];
		PTC[jpt2].F[1]-=dFn[1];

		Sab[0][0]+=dSab[0][0];
		Sab[1][0]+=dSab[1][0];
		Sab[1][1]+=dSab[1][1];
	}
	return(UE);
}

void SHEET :: set_vel(PRTCL *PT, double vx, double vy){

	int k,ipt;

	for(k=0;k<Np;k++){
		ipt=list[k];	
		PT[ipt].v[0]=vx;
		PT[ipt].v[1]=vy;
	}
}

//------------------------------------------------------------------
void Set_Vel(PRTCL *PT, SHEET *st,CNTRL prms){

	int i,k;
	double et[2],rr;
	double v0,vx0,vy0;
	double vmin=prms.vmin,vmax=prms.vmax;
	printf("vmin,vmax[nm/ps]=%lf %lf\n",vmin,vmax);
	init_genrand(10);
	//init_genrand((unsigned)time(NULL));

	for(i=0;i<prms.nst;i++){
		v0=genrand_real1()*(vmax-vmin)+vmin;
		st[i].bbox(PT,prms);
		et[0]=prms.Xc[0]-st[i].x0[0];
		et[1]=prms.Xc[1]-st[i].x0[1];
		rr=sqrt(et[0]*et[0]+et[1]*et[1]);
		et[0]/=rr;
		et[1]/=rr;
		vx0=v0*et[0];
		vy0=v0*et[1];

		st[i].set_vel(PT,vx0,vy0);
	}
	
	vx0=0.0;
	vy0=0.0;
	for(i=0;i<prms.np;i++){
		vx0+=PT[i].v[0];
		vy0+=PT[i].v[1];
	}
	vx0/=prms.np;
	vy0/=prms.np;
	for(i=0;i<prms.np;i++){
		PT[i].v[0]-=vx0;
		PT[i].v[1]-=vy0;
	}

}
//------------------------------------------------------------------
void SHEET::xy2crv(REV rev, PRTCL *PTC){
	int i,ipt;
	Vec2 xf;
	for(i=0;i<Np;i++){
		ipt=list[i];
		xf.set(PTC[ipt].x);
	//	crv.x[i]=PTC[ipt].x[0]+PTC[ipt].rev.;
		crv.y[i]=PTC[ipt].x[1];
	}
};
