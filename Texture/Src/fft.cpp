#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fft.h"
#include <complex>

using namespace std;


int dft(complex<double> *Amp, int N, int isgn){
	
	int i,j;
	double PI=4.0*atan(1.0);
	complex<double> zi(0.0,1.0);
	complex<double> Wij;
	complex<double> *dat=(complex<double> *)malloc(sizeof(complex<double>)*N);

	for(i=0;i<N;i++) dat[i]=Amp[i];

	zi*=isgn;

//		Discrete Fourier Transform
	for(i=0;i<N;i++){
		Amp[i]=complex<double>(0.0,0.0);
	for(j=0;j<N;j++){
		Wij=exp((2.*PI/N*i*j)*zi);
		Amp[i]+=(dat[j]*Wij);
	}
		if(isgn==1) Amp[i]/=N;
	}

	return( ceil(log2(N)) );	
}; 

int fft(complex<double> *Amp, int N, int isgn){

	int p=ceil(log2(N));	
	if(N-pow(2,p)!=0){
		puts(" Error: Signal length must be a power of 2.");
		printf("N=%d, p=%lf\n",N,log2(N));
		puts(" --> process terminated.. ");
		exit(-1);
	}

	int i,k0,k_loc,nlen,irev;
	complex<double> tmp;

//	 Bit-reversal Scramble for FFT 
	for(i=0;i<N;i++){
		k_loc=i;
		k0=0;
		nlen=N;
		while(nlen>2){
			if(k_loc%2 != 0) k0+=(nlen/2); 
			k_loc/=2;
			nlen/=2;
		}
		irev=k0+k_loc;
		//printf("i=%d --> %d\n",i,irev);
		if(i < irev){
			tmp=Amp[i];
			Amp[i]=Amp[irev];
			Amp[irev]=tmp;
		};
	}


//		Butterfly Operation 
	int nseg;
	int j,i1,i2;
	double PI=4.0*atan(1.0);
	complex<double> zi(0.0,1.0);
	complex<double> X1,X2,W1,W2;

	zi*=isgn;

	nseg=N/2;
	nlen=1;
	while(nseg>0){	
		for(i=0;i<nseg;i++){
		for(j=0;j<nlen;j++){
			i1=i*(nlen*2)+j;
			i2=i1+nlen;
			X1=Amp[i1];
			X2=Amp[i2];
			W1=exp(zi*(PI/(double)nlen*j));
			W2=exp(zi*(PI/(double)nlen*(j+nlen)));
			Amp[i1]=X1+W1*X2;
			Amp[i2]=X1+W2*X2;
		}
		}
		nseg/=2;
		nlen*=2;
	}

	if(isgn==1) for(i=0;i<N;i++) Amp[i]/=N;

	return(p);
};

/*
int main(){

	double PI=4.0*atan(1.0);
	int i,q;
	int p=4;
	int N=pow(2,p);
	double *dat=(double *)malloc(sizeof(double)*N);
	complex<double> *Amp=(complex<double> *)malloc(sizeof(complex<double>)*N);

	FILE *f1=fopen("dft.out","w");
	FILE *f2=fopen("fft.out","w");
	FILE *f3=fopen("idft.out","w");
	FILE *f4=fopen("ifft.out","w");
	fprintf(f1,"# No. Re{F}, Im{F}\n");
	fprintf(f2,"# No. Re{F}, Im{F}\n");
	fprintf(f3,"# No. f(i) Re{f(i)}, Im{f(i)}\n");
	fprintf(f4,"# No. f(i) Re{f(i)}, Im{f(i)}\n");

//		Input Signal
	for(i=0;i<N;i++){
		dat[i]=sin(2*PI/N*3*i);
		Amp[i]=complex<double>(dat[i],0.0);
	}

//		DFT
	q=dft(Amp,N,1);
	for(i=0;i<N;i++) fprintf(f1,"%d %lf %lf\n",i,Amp[i].real(),Amp[i].imag());
	puts(" DFT  -->  dft.out");

//		IDFT

	q=dft(Amp,N,-1);
	for(i=0;i<N;i++) fprintf(f3,"%d %lf %lf %lf\n",i,dat[i],Amp[i].real(),Amp[i].imag());
	puts(" IDFT --> idft.out");

//		FFT
	q=fft(Amp,N,1);
	for(i=0;i<N;i++) fprintf(f2,"%d %lf %lf\n",i,Amp[i].real(),Amp[i].imag());
	puts(" FFT  -->  fft.out");

//		IFFT
	q=fft(Amp,N,-1);
	for(i=0;i<N;i++) fprintf(f4,"%d %lf %lf %lf\n",i,dat[i],Amp[i].real(),Amp[i].imag());
	puts(" IFFT --> ifft.out");

	fclose(f1);
	fclose(f2);
	fclose(f3);
	fclose(f4);

	return(0);
};
*/

//------------------------------------------------------------
double DFT_prms::t(int i){ 
	return(ts+dt*i);
};
double DFT_prms::f(int i){
	return(df*i);
};

void DFT_prms::sigma(complex<double> *Amp){
	int i;
	double sig,PI=4.0*atan(1.0);
	double arg;

	for(i=1;i<Nt/2;i++){
		arg=(PI*i)/(100);
		sig=sin(arg)/arg;
		sig=pow(sig,6);
		Amp[i]*=sig;
		Amp[Nt-i]*=sig;
	};	
};
void DFT_prms::diff(complex<double> *Amp){
	complex<double> zi(0.0,1.0);
	complex<double> zw;
	int i;

	for(i=1;i<Nt/2;i++){
		zw=-zi*(dw*i);
		Amp[i]*=zw;
		Amp[Nt-i]*=(-zw);
	};	
};
void DFT_prms::integ(complex<double> *Amp){
	complex<double> zi(0.0,1.0);
	complex<double> zw;
	int i;

	for(i=1;i<Nt/2;i++){
		zw=-zi*(dw*i);
		Amp[i]/=zw;
		Amp[Nt-i]/=(-zw);
	};	
};

void DFT_prms::set_time(double t1, double t2, int n){
	double PI=4.0*atan(1.0);
	ts=t1;
	te=t2;
	Td=te-ts;
	dt=(te-ts)/n;

	df=1./Td;
	dw=2.*PI*df;
	
};
void DFT_prms::set_time(double dtau, int n){
	double PI=4.0*atan(1.0);
	Nt=n;
	dt=dtau;	
	te=Nt*dt;
	Td=te;
	df=1./Td;
	dw=2.*PI*df;
	ts=0.0;
};
//------------------------------------------------------------

