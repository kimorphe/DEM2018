#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "dem.h"

using namespace std;

double vfac(double TK, double T0, double T1, double t_now, double t_start, double t_end){
	double Tset;
	if(t_now > t_end) return(1.0);
	if(t_now < t_start){
		Tset=T1;	
	}else{
		Tset=(T1-T0)*(t_now-t_start)/(t_end-t_start)+T0;
	}
	double Vfac=sqrt(Tset/TK);

	if(Vfac < 0.90) Vfac=0.90;
	if(Vfac > 1.10) Vfac=1.10;
	return(Vfac);
};
void join_chars( char *str1, char *str2 , char *str_out);
void join_chars( char *str1, char *str2 );

int main(){

	double PI=4.0*atan(1.0);
	char fninp[128]="dem.inp";	// General DEM parameters



	char fnptc[128]="ptc.dat";	// Particle Data
	char fnsht[128]="sheet.dat";// Clay Sheet Data
	char fnump[128]="ptc_nums.dat"; // particle number/sheet 
	char fname[128],cbff[128];

	char fnerg[128]="energy.out";
	char fnstr[128]="stress.out";
	char fnptcl[128]="ptcl.out";

//	------------- READ DEM PARAMETERS -------------
	CNTRL prms;
	prms.load(fninp);
	join_chars(prms.Dir,fnerg);
	join_chars(prms.Dir,fnstr);
	join_chars(prms.Dir,fnptcl);
	join_chars(prms.Dir,fnump);

	puts(fnerg);
	puts(fnstr);
	puts(fnptcl);
	puts(fnump);

	FILE *fp,*ftmp;
	FILE *ferg=fopen(fnerg,"w");
	FILE *fstr=fopen(fnstr,"w");

	double Sab[2][2],dS1[2][2],dS2[2][2];
	double m0;
	double kb=1.381e-23; //[J/K] Boltzmann constant
	double Na=6.022e+23; // Abogadro Number
	double TK;	// temperature [K]
	double KE,UE;	// Kinetic & potential energy

	int nsmp=300,ismp=0;
	double Sab_smp[2][2];
	double dsxx,dsxy,dsyy;
	double dsxx0,dsxy0,dsyy0;
	double T0,T1;
	double Tfac=1.0,Vfac;

	int nst;	// number of clay molecules
	int i,j,k,np,imb,ipt,ir0,ir1,i0,i1;
	double Xmin,Xmax,Ymin,Ymax,Wmax[2],tt,Mtot,rho_d,pr;
	double dt; // time increment
	int Nt;	  // number of time steps	
	int Nout,ninc,isum; // number of output times 
	double x1,x2,v1,v2,mj,r1,r2;
	double sig;
	int i0_T_start;

	SHEET *st;
	PRTCL *PTC,pt;
	REV rev;
	WALL wll;

	FILE *fpl_out=fopen(fnptcl,"w");

//	------------  READ WALL DISPALACEMENT DATA -------------
	wll.setup();
//	------------  READ PARTICLE DATA -------------
	fp=fopen(fnptc,"r");
	if(fp == NULL) show_msg(fnptc);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&np);
	printf("np=%d\n",np);

	PTC= (PRTCL *)malloc(sizeof(PRTCL)*np);
	for(i=0;i<np;i++) PTC[i].init();

	fgets(cbff,128,fp);
	int npl=0;
	for(i=0;i<np;i++){
		fscanf(fp,"%le %le %le %d %d %d %le\n",&x1,&x2,&mj,&imb,&i0,&i1,&sig);
		if(i==0){
			Xmin=x1; Xmax=x1;
			Ymin=x2; Ymax=x2;
		}
		if(x1>Xmax) Xmax=x1;
		if(x2>Ymax) Ymax=x2;
		if(x1<Xmin) Xmin=x1;
		if(x2<Ymin) Ymin=x2; 
		PTC[i].setX(x1,x2);
		PTC[i].mobile=imb;
		PTC[i].m=mj*(sig*sig/1.5/1.5);	
		PTC[i].irev[0]=i0;
		PTC[i].irev[1]=i1;
		if( sig > 10.0) npl++;
		//PTC[i].sig=sig;
		PTC[i].sig=1.0;
		PTC[i].sigs[0]=sig;
		PTC[i].sigs[1]=sig;
	}

	int *indx=(int *)malloc(sizeof(int)*npl);
	int jsum=0;
	for(i=0;i<np;i++){
		if(PTC[i].sig >10.0) indx[i]=jsum++;
	}
	printf("npl=%d\n",npl);
	for(i=0;i<npl;i++) printf("indx=%d\n",indx[i]);
	Wmax[0]=Xmax-Xmin;
	Wmax[1]=Ymax-Ymin;
	fclose(fp);

//	------------------ DEM PARAMETERS -------------
	//prms.load(fninp);
	ninc=prms.Nt/prms.Nout;
	prms.np=np;
	if(ninc ==0) ninc=1;
	i0_T_start=int(prms.time_T_start/prms.dt);
	if(i0_T_start < 1) i0_T_start=1;
//	------------  READ AGGREGATE DATA -------------
	puts(fnump);
	ftmp=fopen(fnump,"w");
	fp=fopen(fnsht,"r");
	if(fp == NULL) show_msg(fnsht);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&nst);
	prms.nst=nst;
	st=(SHEET *)malloc(sizeof(SHEET)*nst);
	for(i=0;i<nst;i++){
		fgets(cbff,128,fp);
		fscanf(fp,"%d\n",&j);
 		st[i].init(j);
		st[i].crv.init(j);
		for(k=0;k<st[i].Np;k++){
			fscanf(fp,"%d",st[i].list+k);
			st[i].list[k]--;

			fgets(cbff,2,fp);
			PTC[st[i].list[k]].ist=i;
			PTC[st[i].list[k]].ipt=k;

		}
		st[i].set_prms(prms,PTC);
		fprintf(ftmp,"%d\n",st[i].Np);
	}
	fclose(ftmp);
	fclose(fp);

	Set_Vel(PTC,st,prms);

	if(prms.rstat==1) restart(&prms,np,PTC,nst,st);

	rev.setup(prms,wll);	// setup unit cell
	rev.update(0,prms.dt,wll);
	rev.print();		// print REV parameters
	for(i=0;i<nst;i++) st[i].xy2crv(rev,PTC);


//	----------------------------------------------
	
	Mtot=np*prms.m0;
	m0=prms.m0*1.e24;

	rho_d=Mtot/rev.Vol*1.e27*1.e-03;
	pr=1.0-np*1.0*1.0/rev.Vol;

	save_ptc_data(0,0.0,rho_d,pr,rev,np,PTC,nst,st,prms.Dir);
	printf("\n");
	fprintf(fpl_out,"%d\n",npl);
	for(j=0;j<npl;j++) fprintf(fpl_out,"%lf\n",PTC[indx[j]].sig);
	for(j=0;j<npl;j++) fprintf(fpl_out,"%lf %lf\n",PTC[indx[j]].x[0],PTC[indx[j]].x[1]);

	fprintf(ferg,"# time [ps], kinetic & potential energy [kJ/mol], Cell Size [nm], Temperature [K]\n");
	fprintf(fstr,"# time [ps], Crystal consts:Wd[0],Wd[1][nm], Wd[2],Wd[3][rad], Stress [GPa]: Sxx, Syx, Syy \n");

	double ds11,ds12,ds22;

	Sab_smp[0][0]=0.0;
	Sab_smp[0][1]=0.0;
	Sab_smp[1][1]=0.0;
	isum=0;

	ftmp=fopen("erg.out","w");


	SUBCELL *sbcll;
	for(i=1;i<=prms.Nt;i++){	// Time Step (START) 

		//printf("i=%d\n",i);
		for(j=0;j<np;j++) PTC[j].scale(i,prms.dt,rev);
		wll.disp(prms,i);
		rev.update(i,prms.dt,wll);
		prms.Wd[0]=rev.Wd[0];
		prms.Wd[1]=rev.Wd[1];
		prms.Xa[0]=rev.Xa[0];
		prms.Xa[1]=rev.Xa[1];

		sbcll=(SUBCELL *)malloc(sizeof(SUBCELL)*rev.nsub);
		for(j=0;j<rev.nsub;j++) sbcll[j].setup(rev,prms.sig); 
		regist_ptc(rev,PTC,sbcll,np);

		rho_d=Mtot/rev.Vol*1.e27*1.e-03;
		pr=1.0-np*1.0*1.0/rev.Vol;

		for(j=0;j<np;j++){	// Clear Force Vectors
			PTC[j].F[0]=0.0; 
			PTC[j].F[1]=0.0;
		}
		
		Sab[0][0]=0.0;
		Sab[1][0]=0.0;
		Sab[1][1]=0.0;
		UE=0.0;
		UE+=(VDF(rev,sbcll,PTC,prms.iprd,prms.sig,prms.Eps,dS1)*2.*prms.Eps0);
		//UE+=(VDF_L(rev,PTC,prms.iprd,prms.Eps,dS1,np,npl,indx)*2.*prms.Eps0);
		Sab[0][0]+=dS1[0][0]; 
		Sab[1][0]+=dS1[1][0];
		Sab[1][1]+=dS1[1][1];
		for(j=0;j<nst;j++){ // START_SHEET_j
			UE+=st[j].incN(PTC,rev,dS1)*prms.m0*1.e06; // axial force
			UE+=st[j].incQ(PTC,rev,dS2)*prms.m0*1.e06; // shearing force
			Sab[0][0]+=(dS1[0][0]+dS2[0][0]);
			Sab[1][0]+=(dS1[1][0]+dS2[1][0]);
			Sab[1][1]+=(dS1[1][1]+dS2[1][1]);
		} // END_SHEET_j

		Sab[0][0]*=0.5;
		Sab[1][0]*=0.5;
		Sab[1][1]*=0.5;
		for(j=0;j<np;j++){
			pt=PTC[j];

			Sab[0][0]+=(pt.m*pt.v[0]*pt.v[0]);
			Sab[1][0]+=(pt.m*pt.v[1]*pt.v[0]);
			Sab[1][1]+=(pt.m*pt.v[1]*pt.v[1]);
		}
		Sab[0][0]/=rev.Vol;
		Sab[1][0]/=rev.Vol;
		Sab[1][1]/=rev.Vol;


//			--------- Velocity & Position Vectors Updated  ---------
		for(j=0;j<np;j++){	// PARTICLE j (START)
			PTC[j].incV(prms);
			PTC[j].incX(prms,rev);
		}

//			---------  OUTPUT  ---------
		if(i%ninc ==0){ 
			isum++;
			tt=i*prms.dt;
			save_ptc_data(isum,tt,rho_d,pr,rev,np,PTC,nst,st,prms.Dir);
			printf("( %le[g/cm3] )\n",rho_d);

			for(j=0;j<npl;j++){
				fprintf(fpl_out,"%lf %lf\n",PTC[indx[j]].x[0],PTC[indx[j]].x[1]);
			}
		}

		KE=0.0;
		for(j=0;j<np;j++){
			KE+=PTC[j].KE(prms.m0);
		}
		TK=KE/(1.5*prms.np*kb);


		rev.smooth(Sab,TK);
		KE=KE/np*Na*1.e-03;
		UE=UE/np*Na*1.e-03; // [kJ/mol]
		fprintf(ferg,"%le %le %le %le %le %le %le\n",i*prms.dt,KE,UE,rev.Wd[0],rev.Wd[1],TK,rev.Tb);
		fprintf(fstr,"%le %le %le %le %le ",i*prms.dt,rev.Wd[0],rev.Wd[1],rev.Wd[2],rev.Wd[3]);
		fprintf(fstr,"%le %le %le ",Sab[0][0]*m0,Sab[1][0]*m0,Sab[1][1]*m0);
		fprintf(ftmp,"%le %le",PTC[0].UE[0],PTC[0].UE[1]);
		fprintf(ftmp," %le %le",PTC[20].UE[0],PTC[20].UE[1]);
		fprintf(ftmp," %le %le",PTC[100].UE[0],PTC[100].UE[1]);

		double dUE=VarUE(rev,sbcll,PTC,prms.iprd,prms.sig,prms.Eps,20,0,-0.5); //*2.*prms.Eps0;
		//printf("dUE=%lf\n",)
		fprintf(ftmp," %le\n",dUE);




		if((ismp%nsmp)==0){
			dsxx=rev.sxx-rev.sxxb*m0;
			dsxy=rev.sxy-rev.sxyb*m0;
			dsyy=rev.syy-rev.syyb*m0;
			dsxx0=rev.sxx-Sab_smp[0][0];
			dsxy0=rev.sxy-Sab_smp[0][1];
			dsyy0=rev.syy-Sab_smp[1][1];
			if(ismp > 0) {
				if(fabs(dsxx/(dsxx0+1.e-07)) < 0.95) rev.kxx*=1.1; 
				if(fabs(dsxy/(dsxy0+1.e-07)) < 0.95) rev.kxy*=1.1; 
				if(fabs(dsyy/(dsyy0+1.e-07)) < 0.95) rev.kyy*=1.1; 
				if(fabs(dsxx/(dsxx0+1.e-07)) >= 1.1) rev.kxx*=0.9; 
				if(fabs(dsxy/(dsxy0+1.e-07)) >= 1.1) rev.kxy*=0.9; 
				if(fabs(dsyy/(dsyy0+1.e-07)) >= 1.1) rev.kyy*=0.9; 
			}

		
			Sab_smp[0][0]=rev.sxxb*m0;
			Sab_smp[0][1]=rev.sxyb*m0;
			Sab_smp[1][1]=rev.syyb*m0;
		}
		ismp++;
		rev.sxxb*=m0;
		rev.syyb*=m0;
		rev.sxyb*=m0;
		fprintf(fstr,"%le %le %le \n",rev.sxxb,rev.sxyb,rev.syyb);

		if(prms.T_cntrl==1){
			if(i==i0_T_start){
				T0=rev.Tb;
				printf("T0=%lf\n",T0);
			}

			if(i%100==0){

			if(i >= i0_T_start){
				tt=i*prms.dt;
				Vfac=vfac(rev.Tb, T0, prms.Tset,tt,prms.time_T_start, prms.time_T_end); 
				for(j=0;j<prms.np;j++){
					 PTC[j].v[0]*=Vfac;
					 PTC[j].v[1]*=Vfac;
				}
			}
			}
		}

		//wll.disp(prms,i+1);
		for(j=0;j<rev.nsub;j++) free(sbcll[j].list); // to prevent memory leak
		free(sbcll);

	}	// Time Step (END)

	isum++;
	save_ptc_data(isum,prms.dt*prms.Nt,rho_d,pr,rev,np,PTC,nst,st,prms.Dir);

	return(0);
} 

double dist(double *x, double *y){
	double dx,dy;
	dx=x[0]-y[0];
	dy=x[1]-y[1];
	return sqrt(dx*dx+dy*dy);
}

void print_ptc(
	PRTCL *PTC,	// particle class array 
	int np,		// number of particles
	int isum,	// file number
	double tt,	// time 
	REV rev,	// unit cell 
	double rho_d,	// dry density
	double pr,	// porosity	
	char *Dir	// Output data directory
){
	int j,ir0,ir1;
	char fname[128];
	double x,y;
	FILE *fp;


	sprintf(fname,"x%d.dat",isum);
	join_chars(Dir,fname);
	printf(" %s\n",fname);
	fp=fopen(fname,"w");


	fprintf(fp,"# time [pico sec]\n");
	fprintf(fp,"%lf\n",tt);
	fprintf(fp,"# dry density [g/cm^3], porosity [%%]\n");
	fprintf(fp,"%lf %lf\n",rho_d,pr*100.);
	fprintf(fp,"# boundig box (Xa[2],& size Wd[2])\n");
	fprintf(fp,"%lf %lf\n",rev.Xa[0],rev.Xa[1]);
	fprintf(fp,"%lf %lf\n",rev.Wd[0],rev.Wd[1]);
	fprintf(fp,"# position: x, y, irev[0:1] \n");
	for(j=0;j<np;j++){	// PARTICLE j (START)
		x=PTC[j].x[0]; 
		y=PTC[j].x[1]; 
		if(dist(rev.Xa,PTC[j].x) > 5*(rev.Wd0[0]+rev.Wd0[1]))
		{
			puts("Too large displacement detected!");
			puts(" --> process terminated");
			printf("Wmax=(%lf,%lf)\n",rev.Wd0[0],rev.Wd0[1]);
			printf("PTC[%d].x=(%lf,%lf)\n",j,PTC[j].x[0],PTC[j].x[1]);
			exit (-1); 
		}
		ir0=PTC[j].irev[0];
		ir1=PTC[j].irev[1];
		fprintf(fp,"%lf %lf %d %d\n",x,y,ir0,ir1);
	}
	fclose(fp);
};

void regist_ptc(
	REV rev,
	PRTCL *PTC,
	SUBCELL *sbcll,
	int np
){
	int i,j,ic,jc,l,ipt;
	Vec2 xf;
	for(i=0;i<np;i++){ // Particle 

		if(PTC[i].sig>10.0) continue;

		xf.set(PTC[i].x[0]-rev.Xa[0], PTC[i].x[1]-rev.Xa[1]);
		ic=floor(iprod(xf,rev.ah)*rev.Nh[0]);
		jc=floor(iprod(xf,rev.bh)*rev.Nh[1]);
		l=ic*rev.Nh[1]+jc;	
		if(l>=rev.nsub || l<0){
			puts("no subcell exsits for given l!");
			printf("ic,jc,l=%d %d %d\n",ic,jc,l);
			printf("x,Xa,hd=%lf %lf %lf\n",PTC[i].x[0],rev.Xa[0],rev.hd[0]);
			exit(-1);
		}
		PTC[i].sbcll=l;
		j=sbcll[l].np;	// sbcll[l].np is supposed to be cleared by SUBCELL.setup(...);
		if(j>=sbcll[l].np_max || j<0){
			puts("too many particles for one cell!");
			printf("np=%d,np_max=%d\n",j,sbcll[l].np_max);
			printf("size=%lf %lf\n",rev.hd[0],rev.hd[1]);
			exit(-1);
		}
		sbcll[l].list[j]=i;
		sbcll[l].np++;
	} // Particle
};

void save_ptc_data(
	int isum,	// file number
	double tt, 	// time
	double rho_d, 	// dry density  
	double pr,	// porosity
	REV rev,	// REV data
	int np,		// number of particles
	PRTCL *PTC,	// particle data
	int nst,	// number of sheets
	SHEET *st,	// sheet
	char *Dir	// Output data directory
){
	int j;
	FILE *fp;

//	----------
	char fname[128];
	sprintf(fname,"x%d.dat",isum);
	join_chars(Dir,fname);
	printf("%s ",fname);
	fp=fopen(fname,"w");
	fprintf(fp,"# time [pico sec]\n");
	fprintf(fp,"%lf\n",tt);
	fprintf(fp,"# dry density [g/cm^3], porosity [%%]\n");
	fprintf(fp,"%lf %lf\n",rho_d,pr*100.);
//	----------

	fprintf(fp,"#Xa[1:2]\n");
	fprintf(fp,"%le %le\n",rev.Xa[0],rev.Xa[1]);
	fprintf(fp,"#Wd[1:4]\n");
	for(j=0;j<4;j++) fprintf(fp,"%22.16le ",rev.Wd[j]);
	fprintf(fp,"\n");
	fprintf(fp,"#np\n");
	fprintf(fp,"%d\n",np);
	fprintf(fp,"# irev[1:2] x1   x2    v1   v2\n");
	for(j=0;j<np;j++){	
		fprintf(fp,"%d %d %22.16le %22.16le %22.16le %22.16le\n",
		PTC[j].irev[0],PTC[j].irev[1],PTC[j].x[0], PTC[j].x[1],PTC[j].v[0],PTC[j].v[1]);
	}
	fprintf(fp,"#number of sheets\n");
	fprintf(fp,"%d\n",nst);
	fprintf(fp,"#natural spring lengths");
	for(j=0;j<nst;j++){	
		fprintf(fp,"%22.16le %22.16le\n",st[j].r1,st[j].r2);
	}	
	fclose(fp);

};

void restart(
	CNTRL *prms,
	int np,	// number of particles
	PRTCL *PTC,
	int nst,
	SHEET *st
){

	FILE *fp=fopen((*prms).fnrst,"r");
	char cbff[128];
	int i0,i1,i;
	double x1,x2,v1,v2,r1,r2;
	double Xa[2],Wd[4],Xc[2];


	printf(" reading data from %s\n",(*prms).fnrst);
	if(fp==NULL) show_msg((*prms).fnrst);

	double tt,rho_d,pr;	// not used later
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&tt);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf\n",&rho_d,&pr);

	fgets(cbff,128,fp);
	fscanf(fp,"%le %le\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%le %le %le %le\n",Wd,Wd+1,Wd+2,Wd+3);
		(*prms).Xa[0]=Xa[0]; (*prms).Xa[1]=Xa[1];

		for(i=0;i<4;i++){
			 (*prms).Wd[i]=Wd[i]; 
			 (*prms).Wd0[i]=Wd[i]; 
		}
	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	for(i=0;i<np;i++){
		fscanf(fp,"%d %d %le %le %le %le\n",&i0,&i1,&x1,&x2,&v1,&v2);
		PTC[i].x[0]=x1; PTC[i].x[1]=x2;
		PTC[i].v[0]=v1; PTC[i].v[1]=v2;
		PTC[i].irev[0]=i0; PTC[i].irev[1]=i1;
	}

	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	fgets(cbff,128,fp);
	for(i=0;i<nst;i++){
		fscanf(fp,"%le %le",&r1,&r2);
		st[i].r1=r1; 
		st[i].r2=r2; 
	}
	fclose(fp);
};

void join_chars( char *str1, char *str2 ){
	int nchar=strlen(str1);
	if( str1[nchar-1]=='\n') str1[nchar-1]=int(NULL);
	char tmp[128];

	strcpy(tmp,str1);
	strcat(tmp,str2);
	strcpy(str2,tmp);

};

void join_chars( char *str1, char *str2 , char *str_out){

	int nchar=strlen(str1);
	if( str1[nchar-1]=='\n') str1[nchar-1]=int(NULL);

	strcpy(str_out,str1);
	strcat(str_out,str2);
};
