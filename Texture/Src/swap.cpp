#include <stdio.h>

void swap(int *i, int *j){
	int tmp;
	tmp=*i;
	*i=*j;
	*j=tmp;
}
int ring(int i, int N){
	if(i>=0){
	       return(i%N);
	}else{
		while(i<0) i+=N;
		return(i%N);
	}
}
int main(){

	int i,j;
	i=10;
	j=20;
	printf("(i,j)=%d,%d\n",i,j);
	swap(&i,&j);
	printf("(i,j)=%d,%d\n",i,j);

	int N=4;
	for(int i=-10;i<5;i++){
		printf("i=%d, i(mod %d)= %d, ring(%d)=%d\n",i,N,i%N,i,ring(i,N));
	}



	return(0);
};
