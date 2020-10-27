#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <x86intrin.h>
#include <stdint.h>


 void SB_Comba(int  *f, int *g, int *h)
{
    #define p 24

    int fg[p+p-1];
    int result;
    int i,j;
    
    for (i = 0;i < p;++i) {
        result = 0;
        for (j = 0;j <= i;++j) 
        result = result+f[j]*g[i-j];
         fg[i] = result;
    }
    
    
    
    for (i = p;i < p+p-1;++i) {
        result = 0;
        for (j = i-p+1;j < p;++j) 
        result = (result+f[j]*g[i-j]);
        fg[i] = result;
    }
    
  
    for (i = 0;i < p+p-1;++i) h[i] = fg[i]%3;
}

void sum(int *a, int *b, int *sum, int n){

    for(int i=0; i<n; ++i)
        sum[i] = a[i] + b[i];
}


void ka2_real(int a[], int b[], int c[], int n){
    
    #define KA2_REAL_MAX_AB_SIZE 384 // n/2
    #define KA2_REAL_MAX_P_SIZE 767 // n-1

    /* base case */
    if(n == 24){
        SB_Comba(a, b, c);
        return;
    }

    int n2 = n/2;

    int *a0 = a;
    int *a1 = a+n2;
    int *b0 = b;
    int *b1 = b+n2;

    int sum_a0_a1[KA2_REAL_MAX_AB_SIZE];
    int sum_b0_b1[KA2_REAL_MAX_AB_SIZE];

    sum(a0, a1, sum_a0_a1, n2);
    sum(b0, b1, sum_b0_b1, n2);

    int p1[KA2_REAL_MAX_P_SIZE];
    int p2[KA2_REAL_MAX_P_SIZE];
    int p3[KA2_REAL_MAX_P_SIZE];
    memset(p1, 0, sizeof(p1));
    memset(p2, 0, sizeof(p2));
    memset(p3, 0, sizeof(p3));

    ka2_real(a0, b0, p1, n2);
    ka2_real(sum_a0_a1, sum_b0_b1, p2, n2);
    ka2_real(a1, b1, p3, n2);


    for(int i=0; i<n-1; ++i){
        c[i] += p1[i];
        c[i+n2] -= p3[i];
    }

    for(int i=n+n2-2; i>=0 ; --i){
        c[i+n2] -= c[i];
    }

    for(int i=0; i<n-1; ++i){
        c[i+n2] += p2[i];
    }

}





  

   
int main(){


  int *a = malloc(768*sizeof(int));
    int *b = malloc(768*sizeof(int));
    int *c = calloc((2*768-1), sizeof(int));


     int i;
////////////////////////reading a and b from the file inp768test
	FILE *myFile;
    myFile = fopen("inp768test", "r");
    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (0);}

    for (i = 0; i < 768; i++){
        fscanf(myFile, "%d", &a[i] );
       
    }
    
    for (i = 0; i < 768; i++){
        fscanf(myFile, "%d", &b[i] );
       
    }
    fclose(myFile);
    

        unsigned long long   LastCycleCount = _rdtsc();
    	clock_t start = clock ();	     
   		for (i=0;i<99999;i++){
        ka2_real(a, b, c, 768);
        memset(c, 0, sizeof(int)*(2*768-1));
    	}
    	ka2_real(a, b, c, 768);
	    unsigned long long   EndCycleCount = _rdtsc();
        unsigned long long   CyclesElapsed = EndCycleCount - LastCycleCount;
        CyclesElapsed = CyclesElapsed/ 100000;
   
        double Timelapsed=(clock()-start)/(double) CLOCKS_PER_SEC;
        Timelapsed=Timelapsed/100000;
 /////////////////////////////////////////////write on c_a3 file	
	FILE *outFile;
    outFile = fopen("c_hybrid_1", "w");
    if (outFile == NULL)
	{
    	printf("Cannot Open File\n");
        exit (0);
    }
    
   for (i = 0; i < 2*768-1; i++)
        fprintf(outFile, "%d ", ((c[i])%3+3)%3 );
        
	fclose(outFile);
/////////////////////////////////////////////////////////////printing on console
	for( i=0; i<2*768-1; i++)
		printf("%d ", ((c[i])%3+3)%3);
    printf("\n\n\n");
	printf("Cycles: %llu\n", CyclesElapsed);
    printf("Multiplication Time:");
    printf( "%lf",Timelapsed);
    printf("sec");
    printf("\n \n ");

	free(a);
	free(b);
	free(c);

	return 0;
}
