#include <stdio.h>
#include <math.h>
#include "domain.h"

void RECT :: setup(
	double x1[2], //endpoint 1
	double x2[2], //endpoint 2
	double t,	// thickness
	int endcap	// endcap (0:none, 1:add)
){
	double e1[2];	// unit tangential vector ( x1 -> x2 )
	double e2[2];	// unit normal vector 
	double d2=0.0,t2=t*0.5;;
	int i;

	for(i=0;i<2;i++){
		e1[i]=x2[i]-x1[i];
		d2+=e1[i]*e1[i];
	}

	d2=sqrt(d2);
	e1[0]/=d2;	
	e1[1]/=d2;
	e2[0]=-e1[1];
	e2[1]= e1[0];

	switch(endcap){
	case 0 :
		x[0]=x1[0]-t2*e2[0];
		y[0]=x1[1]-t2*e2[1];

		x[1]=x2[0]-t2*e2[0];
		y[1]=x2[1]-t2*e2[1];

		x[2]=x2[0]+t2*e2[0];
		y[2]=x2[1]+t2*e2[1];

		x[3]=x1[0]+t2*e2[0];
		y[3]=x1[1]+t2*e2[1];

		break;
	case 1 :
		x[0]=x1[0]-t2*(e1[0]+e2[0]);
		y[0]=x1[1]-t2*(e1[1]+e2[1]);

		x[1]=x2[0]+t2*(e1[0]-e2[0]);
		y[1]=x2[1]+t2*(e1[1]-e2[1]);

		x[2]=x2[0]+t2*(e1[0]+e2[0]);
		y[2]=x2[1]+t2*(e1[1]+e2[1]);

		x[3]=x1[0]+t2*(-e1[0]+e2[0]);
		y[3]=x1[1]+t2*(-e1[1]+e2[1]);

		break;
	};

	xmin=x[0];
	xmax=x[0];
	ymin=y[0];
	ymax=y[0];

	for( i=1;i<4;i++){
		if(xmin > x[i]) xmin=x[i];
		if(xmax < x[i]) xmax=x[i];
		if(ymin > y[i]) ymin=y[i];
		if(ymax < y[i]) ymax=y[i];
	}

};

void RECT :: setup(
	double xc[2], //	center
	double b,  //	width
	double h,  //	height
	double th  //	counterclock-wise angle in degree
){

	int i;
	double xd[4],yd[4];
	double R[2][2];
	double pi=4.0*atan(1.0);

	th=th/180.0*pi;
	R[0][0]= cos(th);
	R[0][1]=-sin(th);
	R[1][0]= sin(th);
	R[1][1]= cos(th);

	xd[0]=-0.5*b; yd[0]=-0.5*h;
	xd[1]= 0.5*b; yd[1]=-0.5*h;
	xd[2]= 0.5*b; yd[2]= 0.5*h;
	xd[3]=-0.5*b; yd[3]= 0.5*h;

	for(i=0;i<4;i++){
		x[i]=R[0][0]*xd[i]+R[0][1]*yd[i]+xc[0];
		y[i]=R[1][0]*xd[i]+R[1][1]*yd[i]+xc[1];
	}

};

void RECT ::out(char *fname){
	FILE *fp=fopen(fname,"w");

	for(int i=0;i<4;i++){
		fprintf(fp,"%lf %lf\n",x[i],y[i]);
	}
	fprintf(fp,"%lf %lf\n",x[0],y[0]);
	fprintf(fp,"\n");

	fclose(fp);
};

int RECT :: isin(double xf[2]){

	double e1[2],e2[2],d2;
	int i,j,i2,in;
	double sgn;
	

	in=1;
	for(i=0;i<4;i++){
		i2=((i+1)%4);	
		
		e1[0]=x[i2]-x[i];
		e1[1]=y[i2]-y[i];
		d2=e1[0]*e1[0]+e1[1]*e1[1];
		d2=sqrt(d2);
		e1[0]/=d2;
		e1[1]/=d2;


		e2[0]=-e1[1];
		e2[1]= e1[0];

		sgn=e2[0]*(xf[0]-x[i])+e2[1]*(xf[1]-y[i]);

		if(sgn < 0.0 ){
			in=0;
			break;
		}
		
	}
	return(in);
};

