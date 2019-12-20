#include <complex>
using namespace std;
int dft(complex<double> *Amp, int N, int isgn);
int fft(complex<double> *Amp, int N, int isgn);
class DFT_prms{
	public:
		int Nt;
		int p;	// p**2=Nt
		double dt;
		double dw;
		double df;
		double ts,te;
		double Td;
		double fmax,wmax;
		void set_time(double t1, double t2, int n);
		void set_time(double dtau, int n);
		double t(int i); // return time_i
		double f(int i); // return freq_i
		double w(int i); // return omega_i
		void diff(complex<double> *Amp);
		void integ(complex<double> *Amp);
		void sigma(complex<double> *Amp);
	private:
};
