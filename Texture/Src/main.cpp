#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "domain.h"

int main(int argc, char *argv[]){

	int Ncell;
	char fname[128]="kcell.inp";
	if(argc > 1) strcpy(fname,argv[1]);

	char fnout[128]="kcell.dat";
	if(argc > 2) strcpy(fnout,argv[2]);

	Dom2D dom(fname);


	dom.perfo(fname);
	dom.perfo_rec(fname);
//	dom.mscs(fname);
	dom.out_kcell(fnout);

	int k0=dom.count(-1);
	int ks=dom.count(0);
	int kf=dom.count(1);

	printf("solid cells: %d\n",ks);
	printf("fluid cells: %d\n",kf);
	printf(" void cells: %d\n",k0);
	printf("      total: %d\n",dom.Ndiv[0]*dom.Ndiv[1]);
	printf("------------------------\n");
	FILE *fp=fopen("ncell.dat","a");
		fprintf(fp,"%d %d %d %d\n",ks,kf,k0,ks+kf+k0);

	double rho_s=2.4;	// mass density (solid phase) [g/cm^3]
	double rho_w=1.0;	// mass density (fulid phase) [g/cm^3]

	Ncell=dom.Ndiv[0]*dom.Ndiv[1];
	double rho=(ks*rho_s+kf*rho_w)/Ncell; // dry density
	double rho_dry=(ks*rho_s)/Ncell;	// wet density (apparent mass density)

	printf("wet density  =%lf [g/cm^3]\n",rho);
	printf("dry density  =%lf [g/cm^3]\n",rho_dry);
	printf("   porosity n=%lf [percent] \n",(k0+kf)*100./Ncell);
	printf(" void ratio e=%lf\n",(k0+kf)/(double)ks);
	printf("moisture content w=%lf\n",(rho_w*kf)/(rho_s*ks));
	printf("degree of saturation Sr=%lf\n",kf/(double)(kf+k0)*100.);
	fclose(fp);
}
