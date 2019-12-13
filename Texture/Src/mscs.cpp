#include <stdio.h>
#include <math.h>
#include "mscs.h"



double dist2d(double* x, double* y){

	double dx=x[0]-y[0];
	double dy=x[1]-y[1];

	return( sqrt(dx*dx+dy*dy));
}

Mscs :: Mscs(){};
void Mscs :: setup(
	double X1[2],	// center 1
	double X2[2],	// center 2
	double Rd,	// radius of the circles
	double alp,	// meniscus width [deg] 
	double th	// contact angle [deg]

){
	int i;
	double EPS=1.e-07;
	double PI=4.0*atan(1.0);

	double Ex,Ey;
	double dnm;

	for( i=0;i<2;i++){
		x1[i]=X1[i]; x2[i]=X2[i];
	}
	rd=Rd;

	alpha=alp/180.*PI;
	theta=th/180.*PI;

	d12=dist2d(x1,x2);

	eh[0]=(x2[0]-x1[0])/d12;
	eh[1]=(x2[1]-x1[1])/d12;

	nh[0]=-eh[1];
	nh[1]= eh[0];

	rho=0.5*d12-rd*cos(alpha);
	dnm=cos(alpha+theta);
	if(fabs(dnm) < EPS) dnm=EPS;
	rho/=dnm;

	printf("rho=%lf\n",rho);

	Ex=rd*cos(alpha)+rho*cos(alpha+theta);
	Ey=rd*sin(alpha)+rho*sin(alpha+theta);

	x0p[0]=x1[0]+Ex*eh[0]+Ey*nh[0];
	x0p[1]=x1[1]+Ex*eh[1]+Ey*nh[1];
	
	Ey=rd*sin(-alpha)+rho*sin(-alpha-theta);
	x0m[0]=x1[0]+Ex*eh[0]+Ey*nh[0];
	x0m[1]=x1[1]+Ex*eh[1]+Ey*nh[1];
};

int Mscs :: IsIn(double xf[2]){

	double dp=dist2d(xf,x0p);
	double dm=dist2d(xf,x0m);
	double xfe,xfn;

	double sgn=1.0;
	if(rho < 0.0) sgn=-1.0;
	int iflg=0;
	if( (dp - fabs(rho))*sgn < 0.0) iflg=1;
	if( (dm - fabs(rho))*sgn < 0.0) iflg=1;

	xfe=(xf[0]-x1[0])*eh[0]+(xf[1]-x1[1])*eh[1];
	xfn=(xf[0]-x1[0])*nh[0]+(xf[1]-x1[1])*nh[1];
	if(xfe < rd*cos(alpha)) iflg=1;
	if(xfe > d12-rd*cos(alpha)) iflg=1;

	if(rho > 0.0){
		if(fabs(xfn) >= rd*sin(alpha)) iflg=1;
	}

	if(dist2d(xf,x1) <= rd) iflg=1;
	if(dist2d(xf,x2) <= rd) iflg=1;

	return(iflg);

};

void Mscs::show(){
	printf("X1=(%lf, %lf)\n",x1[0],x1[1]);
	printf("X2=(%lf, %lf)\n",x2[0],x2[1]);
	printf("rho=%lf \n",rho);
	printf("x0p=(%lf, %lf)\n",x0p[0],x0p[1]);
	printf("x0m=(%lf, %lf)\n",x0m[0],x0m[1]);
};

