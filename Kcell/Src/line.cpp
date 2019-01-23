#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "domain.h"
#include "spline.h"

int main(int argc, char *argv[]){

	int Ncell;
	char fname[128]="kcell.inp";
	if(argc > 1) strcpy(fname,argv[1]);

	char fnout[128]="kcell.dat";
	if(argc > 2) strcpy(fnout,argv[2]);

	Dom2D dom(100,120);

	dom.Xa[0]=0.0;
	dom.Xa[1]=0.0;
	dom.Xb[0]=10.0;
	dom.Xb[1]=12.0;

	dom.set_dx();

	double x1[2],x2[2];

	x1[0]=2.0; x1[1]=3.0;
	x2[0]=7.0; x2[1]=8.0;

	int lw=3;
	dom.draw_line(x1,x2,1,lw);

	x1[0]=2.0; x1[1]=8.0;
	x2[0]=7.0; x2[1]=3.0;
	dom.draw_line(x1,x2,1,lw);

	/*
	x1[0]=-2.0; x1[1]=5.0;
	x2[0]=14.0; x2[1]=6.0;
	dom.draw_line(x1,x2,1,lw);

	x1[0]=-5.0; x1[1]=6.0;
	x2[0]=10.0; x2[1]=13.0;
	dom.draw_line(x1,x2,1,lw);

	x1[0]=-5.0; x1[1]=2.0;
	x2[0]=12.0; x2[1]=2.0;
	dom.draw_line(x1,x2,1,lw);

	x1[0]=8.0; x1[1]=13.0;
	x2[0]=8.0; x2[1]=-2.0;
	dom.draw_line(x1,x2,1,lw);
	*/

	double th,dth;
	double pi=4.0*atan(1.0);
	int N=361;
	dth=pi*2.0/(N-1);
	double r1=5.0;
	double r2=6.0;
	for(int i=0;i<N;i++){
		th=i*dth;
		x1[0]=r1*cos(th);
		x1[1]=r1*sin(th);
		x2[0]=r2*cos(th);
		x2[1]=r2*sin(th);
		dom.draw_line(x1,x2,2,1);
	}

	dom.out_kcell(fnout);

	int k0=dom.count(-1);
	int ks=dom.count(0);
	int kf=dom.count(1);

	printf("solid cells: %d\n",ks);
	printf("fluid cells: %d\n",kf);
	printf(" void cells: %d\n",k0);
	printf("      total: %d\n",dom.Ndiv[0]*dom.Ndiv[1]);
	printf("------------------------\n");

	return(0);
}
