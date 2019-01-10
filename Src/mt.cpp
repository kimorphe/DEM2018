#include <stdio.h>
#include <random>

int main(){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> nums(1,25);
	std::uniform_real_distribution<double> dbs(-1.0,1.0); 

	for(int i=0;i<1000;i++){
		printf("%d %lf\n",nums(mt),dbs(mt));
	}
	return(0);
};
