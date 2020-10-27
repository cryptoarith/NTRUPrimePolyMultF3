#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <x86intrin.h>
#include <stdint.h>

struct complex{

	int re;
	int im;
};


void copy(struct complex *source, struct complex *destination, int from, int to){

	for(int i=from, j=0; i<to; ++i, ++j)
		destination[j] = source[i];
}

void sum(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re + b[i].re;
		target[i].im = a[i].im + b[i].im;
	}
}

void sum_neg_neg(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = - a[i].re - b[i].re;
		target[i].im = - a[i].im - b[i].im;
	}
}

void sum_multiply_w(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = - a[i].im - b[i].im;
		target[i].im = a[i].re + b[i].re;
	}
}

void difference(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re - b[i].re;
		target[i].im = a[i].im - b[i].im;
	}
}

void multiplication(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		int ac = a[i].re * b[i].re;
		int bd = a[i].im * b[i].im;
		target[i].re = ac - bd;
		target[i].im = (a[i].re + a[i].im)*(b[i].re + b[i].im) - ac - bd;
	}
}

void negation(struct complex *a, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = -a[i].re;
		target[i].im = -a[i].im;
	}
}


void sum_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re + b[i].re;

	}
}

void sum_neg_neg_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = - a[i].re - b[i].re;

	}
}

void sum_re_im_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re + b[i].im;

	}
}

void sum_re_neg_im_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re - b[i].im;

	}
}

void sum_im_im_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].im + b[i].im;

	}
}

void sum_im_neg_im_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].im - b[i].im;

	}
}

void difference_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re - b[i].re;

	}
}

void multiplication_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re * b[i].re;

	}
}

void negation_real(struct complex *a, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = -a[i].re;

	}
}

void create_term(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re - b[i].im;
		target[i].im = a[i].im + b[i].re;
	}
}

void create_term_neg(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re + b[i].im;
		target[i].im = a[i].im - b[i].re;
	}
}

void phase_change(struct complex *a, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re;
		target[i].im = -a[i].im;		
	}
}




void ka2_real(struct complex a[], struct complex b[], struct complex c[], int n);
void lt_real(struct complex a[], struct complex b[], struct complex c[], int n);


void sb_comba_real(struct complex a[], struct complex b[], struct complex c[], int n){

	int result;

	for(int i=0; i<n; ++i){
		result = 0;
		for(int j=0; j<=i; ++j)
			result += a[j].re*b[i-j].re;
		c[i].re = result;
	}

	for(int i=n; i<2*n-1; ++i){
		result = 0;
		for(int j=i-n+1; j<n; ++j)
			result += a[j].re*b[i-j].re;
		c[i].re = result;
	}

}

void ub_real(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define UB_REAL_MAX_AB_SIZE 10 // ceil(n/2)
	#define UB_REAL_MAX_P_SIZE 19 // n

	// removed base case

	int n0 = (n+1)/2;
	int n1 = n-n0;

	struct complex *a0 = a;
	struct complex *a1 = a+n0;
	struct complex *b0 = b;
	struct complex *b1 = b+n0;

	struct complex sum_a0_a1[UB_REAL_MAX_AB_SIZE];
	struct complex sum_b0_b1[UB_REAL_MAX_AB_SIZE];

	// assumed n is odd.
	sum(a0, a1, sum_a0_a1, n1);
	sum_a0_a1[n0-1] = a0[n0-1];
	sum(b0, b1, sum_b0_b1, n1);
	sum_b0_b1[n0-1] = b0[n0-1];


	struct complex p1[UB_REAL_MAX_P_SIZE];
	struct complex p2[UB_REAL_MAX_P_SIZE];
	struct complex p3[UB_REAL_MAX_P_SIZE];
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));

	ka2_real(a0, b0, p1, n0);
	ka2_real(sum_a0_a1, sum_b0_b1, p2, n0);
	lt_real(a1, b1, p3, n1);


	sum_real(c, p1, c, 2*n0-1);
	difference_real(c+n0, p3, c+n0, 2*n1-1);

	for(int i=n+n1-2; i>=0 ; --i){
		c[i+n0].re -= c[i].re;
		c[i+n0].im -= c[i].im;
	}

	sum_real(c+n0, p2, c+n0, 2*n0-1);
	
}


void a2_complex(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define A2_COMPLEX_MAX_AB_SIZE 38 // n
	#define A2_COMPLEX_MAX_P_SIZE 75 // 2*n-1

	// removed base case

	int p_size = 2*n-1;

	struct complex a0[A2_COMPLEX_MAX_AB_SIZE];
	struct complex a1[A2_COMPLEX_MAX_AB_SIZE];
	struct complex b0[A2_COMPLEX_MAX_AB_SIZE];
	struct complex b1[A2_COMPLEX_MAX_AB_SIZE];

	for(int i=0; i<n; ++i){
		a0[i].re = a[i].re;
		a0[i].im = 0;
		a1[i].re = a[i].im;
		a1[i].im = 0;
		b0[i].re = b[i].re;
		b0[i].im = 0;
		b1[i].re = b[i].im;
		b1[i].im = 0;
	}

	struct complex a0_p_a1[A2_COMPLEX_MAX_AB_SIZE];
	struct complex b0_p_b1[A2_COMPLEX_MAX_AB_SIZE];

	sum(a0, a1, a0_p_a1, n);
	sum(b0, b1, b0_p_b1, n);

	struct complex a0_b0[A2_COMPLEX_MAX_P_SIZE];
	struct complex a1_b1[A2_COMPLEX_MAX_P_SIZE];
	struct complex p[A2_COMPLEX_MAX_P_SIZE];
	memset(a0_b0, 0, sizeof(a0_b0));
	memset(a1_b1, 0, sizeof(a1_b1));
	memset(p, 0, sizeof(p));

	ka2_real(a0, b0, a0_b0, n);
	ka2_real(a1, b1, a1_b1, n);
	ka2_real(a0_p_a1, b0_p_b1, p, n);

	for(int i=0; i<p_size; ++i){
		c[i].re += a0_b0[i].re - a1_b1[i].re;
		c[i].im += p[i].re - a0_b0[i].re - a1_b1[i].re;
	}

}

void ka2_real(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define KA2_REAL_MAX_AB_SIZE 76 // n/2
	#define KA2_REAL_MAX_P_SIZE 151 // n-1

	/* base cases */
	if(n == 19){
		ub_real(a, b, c, n);
		return;
	}

	else if(n == 4 || n == 5){
		sb_comba_real(a, b, c, n);
		return;
	}

	int n2 = n/2;

	struct complex *a0 = a;
	struct complex *a1 = a+n2;
	struct complex *b0 = b;
	struct complex *b1 = b+n2;

	struct complex sum_a0_a1[KA2_REAL_MAX_AB_SIZE];
	struct complex sum_b0_b1[KA2_REAL_MAX_AB_SIZE];

	sum(a0, a1, sum_a0_a1, n2);
	sum(b0, b1, sum_b0_b1, n2);

	struct complex p1[KA2_REAL_MAX_P_SIZE];
	struct complex p2[KA2_REAL_MAX_P_SIZE];
	struct complex p3[KA2_REAL_MAX_P_SIZE];
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));

	if(n == 76){
		a2_complex(a0, b0, p1, n2);
		a2_complex(sum_a0_a1, sum_b0_b1, p2, n2);
		a2_complex(a1, b1, p3, n2);
	}

	else{
		ka2_real(a0, b0, p1, n2);
		ka2_real(sum_a0_a1, sum_b0_b1, p2, n2);
		ka2_real(a1, b1, p3, n2);
	}
	


	sum_real(c, p1, c, n-1);
	difference_real(c+n2, p3, c+n2, n-1);

	for(int i=n+n2-2; i>=0 ; --i){
		c[i+n2].re -= c[i].re;
	}

	sum_real(c+n2, p2, c+n2, n-1);

}

void n1_complex(struct complex a[], struct complex b[], struct complex c[], int n){

	#define N1_COMPLEX_MAX_AB_SIZE 38 // n/4
	#define N1_COMPLEX_MAX_P_SIZE 75 // 2*n/4-1

	// removed base case

	/* Assumed n = 4^k for some non-negative integer k */
	int n4 = n/4;
	int p_size = 2*n4-1;


	struct complex *a0 = a;
	struct complex *a1 = a+n4;
	struct complex *a2 = a+2*n4;
	struct complex *a3 = a+3*n4;
	struct complex *b0 = b;
	struct complex *b1 = b+n4;
	struct complex *b2 = b+2*n4;
	struct complex *b3 = b+3*n4;

	struct complex ra1[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra2[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra3[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra4[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra5[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra6[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra7[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra8[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra9[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra10[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra11[N1_COMPLEX_MAX_AB_SIZE];
	struct complex ra12[N1_COMPLEX_MAX_AB_SIZE];

	struct complex rb1[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb2[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb3[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb4[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb5[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb6[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb7[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb8[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb9[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb10[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb11[N1_COMPLEX_MAX_AB_SIZE];
	struct complex rb12[N1_COMPLEX_MAX_AB_SIZE];


	sum(a1, a3, ra1, n4);
	sum(ra1, a0, ra2, n4);
	sum(ra2, a2, ra3, n4);
	difference(a0, ra1, ra4, n4);
	difference(a1, a2, ra5, n4);
	difference(ra5, a3, ra6, n4);
	sum_neg_neg(a1, a2, ra7, n4);
	sum(ra7, a3, ra8, n4);
	create_term(ra2, ra6, ra9, n4);
	create_term(ra4, ra8, ra10, n4);
	create_term_neg(ra2, ra6, ra11, n4);
	create_term_neg(ra4, ra8, ra12, n4);

	sum(b1, b3, rb1, n4);
	sum(rb1, b0, rb2, n4);
	sum(rb2, b2, rb3, n4);
	difference(b0, rb1, rb4, n4);
	difference(b1, b2, rb5, n4);
	difference(rb5, b3, rb6, n4);
	sum_neg_neg(b1, b2, rb7, n4);
	sum(rb7, b3, rb8, n4);
	create_term(rb2, rb6, rb9, n4);
	create_term(rb4, rb8, rb10, n4);
	create_term_neg(rb2, rb6, rb11, n4);
	create_term_neg(rb4, rb8, rb12, n4);

	struct complex p0[N1_COMPLEX_MAX_P_SIZE];
	struct complex p1[N1_COMPLEX_MAX_P_SIZE];
	struct complex p2[N1_COMPLEX_MAX_P_SIZE];
	struct complex p3[N1_COMPLEX_MAX_P_SIZE];
	struct complex p4[N1_COMPLEX_MAX_P_SIZE];
	struct complex p5[N1_COMPLEX_MAX_P_SIZE];
	struct complex p6[N1_COMPLEX_MAX_P_SIZE];
	memset(p0, 0, sizeof(p0));
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));
	memset(p4, 0, sizeof(p4));
	memset(p5, 0, sizeof(p5));
	memset(p6, 0, sizeof(p6));

	a2_complex(a0, b0, p0, n4);
	a2_complex(ra3, rb3, p1, n4);
	a2_complex(ra9, rb9, p2, n4);
	a2_complex(ra11, rb11, p3, n4);	
	a2_complex(ra10, rb10, p4, n4);
	a2_complex(ra12, rb12, p5, n4);
	a2_complex(a3, b3, p6, n4);

	struct complex u1[N1_COMPLEX_MAX_P_SIZE];
	struct complex u2[N1_COMPLEX_MAX_P_SIZE];
	struct complex u3[N1_COMPLEX_MAX_P_SIZE];
	struct complex u4[N1_COMPLEX_MAX_P_SIZE];
	struct complex u5[N1_COMPLEX_MAX_P_SIZE];
	struct complex u6[N1_COMPLEX_MAX_P_SIZE];
	struct complex u7[N1_COMPLEX_MAX_P_SIZE];
	struct complex u8[N1_COMPLEX_MAX_P_SIZE];
	struct complex u9[N1_COMPLEX_MAX_P_SIZE];
	struct complex u10[N1_COMPLEX_MAX_P_SIZE];
	struct complex u11[N1_COMPLEX_MAX_P_SIZE];
	struct complex u12[N1_COMPLEX_MAX_P_SIZE];
	struct complex u13[N1_COMPLEX_MAX_P_SIZE];

	sum_neg_neg(p0, p1, u1, p_size);
	sum(p4, p5, u2, p_size);
	sum(p2, p3, u3, p_size);
	difference(p2, p3, u4, p_size);
	difference(p4, p5, u5, p_size);
	sum_neg_neg(u3, u2, u6, p_size);
	sum(u1, u6, u7, p_size);
	difference(u7, p6, u8, p_size);
	difference(u5, u4, u9, p_size);
	difference(u2, u3, u10, p_size);
	sum(u4, u5, u11, p_size);
	sum(u1, u2, u12, p_size);
	difference(u12, p6, u13, p_size);


	struct complex c1[N1_COMPLEX_MAX_P_SIZE];
	struct complex c2[N1_COMPLEX_MAX_P_SIZE];
	struct complex c3[N1_COMPLEX_MAX_P_SIZE];
	struct complex c4[N1_COMPLEX_MAX_P_SIZE];
	struct complex c5[N1_COMPLEX_MAX_P_SIZE];

	struct complex *c0 = p0;

	create_term_neg(u8, u4, c1, p_size);
	create_term(p6, u11, c2, p_size);
	create_term(u10, u9, c3, p_size);
	sum(p0, u6, c4, p_size);
	create_term(u13, u11, c5, p_size);

	struct complex *c6 = p6;

	sum(c, c0, c, p_size);
	sum(c+n4, c1, c+n4, p_size);
	sum(c+2*n4, c2, c+2*n4, p_size);
	sum(c+3*n4, c3, c+3*n4, p_size);
	sum(c+4*n4, c4, c+4*n4, p_size);
	sum(c+5*n4, c5, c+5*n4, p_size);
	sum(c+6*n4, c6, c+6*n4, p_size);


	/*for(int i=0; i<2*n-1; i++){
		c[i].re = ((c[i].re%3)+3)%3;
		c[i].im = ((c[i].im%3)+3)%3;
	}*/

}

void v1_real(struct complex a[], struct complex b[], struct complex c[], int n){

	#define V1_REAL_MAX_AB_SIZE 152 // n/5
	#define V1_REAL_MAX_P_SIZE 303 // 2*n/5-1


	// removed base case

	/* Assumed n = 5^k for some non-negative integer k */
	int n5 = n/5;
	int p_size = 2*n5-1;


	struct complex *a0 = a;
	struct complex *a1 = a+n5;
	struct complex *a2 = a+2*n5;
	struct complex *a3 = a+3*n5;
	struct complex *a4 = a+4*n5;
	struct complex *b0 = b;
	struct complex *b1 = b+n5;
	struct complex *b2 = b+2*n5;
	struct complex *b3 = b+3*n5;
	struct complex *b4 = b+4*n5;

	struct complex ra0[V1_REAL_MAX_AB_SIZE];
	struct complex ra1[V1_REAL_MAX_AB_SIZE];
	struct complex ra2[V1_REAL_MAX_AB_SIZE];
	struct complex ra3[V1_REAL_MAX_AB_SIZE];
	struct complex ra4[V1_REAL_MAX_AB_SIZE];
	struct complex ra5[V1_REAL_MAX_AB_SIZE];
	struct complex ra6[V1_REAL_MAX_AB_SIZE];
	struct complex ra7[V1_REAL_MAX_AB_SIZE];
	struct complex ra8[V1_REAL_MAX_AB_SIZE];
	struct complex ra9[V1_REAL_MAX_AB_SIZE];
	struct complex ra10[V1_REAL_MAX_AB_SIZE];
	struct complex ra11[V1_REAL_MAX_AB_SIZE];
	struct complex ra12[V1_REAL_MAX_AB_SIZE];
	struct complex ra13[V1_REAL_MAX_AB_SIZE];
	struct complex ra14[V1_REAL_MAX_AB_SIZE];
	struct complex ra15[V1_REAL_MAX_AB_SIZE];
	struct complex ra16[V1_REAL_MAX_AB_SIZE];
	struct complex ra17[V1_REAL_MAX_AB_SIZE];

	struct complex rb0[V1_REAL_MAX_AB_SIZE];
	struct complex rb1[V1_REAL_MAX_AB_SIZE];
	struct complex rb2[V1_REAL_MAX_AB_SIZE];
	struct complex rb3[V1_REAL_MAX_AB_SIZE];
	struct complex rb4[V1_REAL_MAX_AB_SIZE];
	struct complex rb5[V1_REAL_MAX_AB_SIZE];
	struct complex rb6[V1_REAL_MAX_AB_SIZE];
	struct complex rb7[V1_REAL_MAX_AB_SIZE];
	struct complex rb8[V1_REAL_MAX_AB_SIZE];
	struct complex rb9[V1_REAL_MAX_AB_SIZE];
	struct complex rb10[V1_REAL_MAX_AB_SIZE];
	struct complex rb11[V1_REAL_MAX_AB_SIZE];
	struct complex rb12[V1_REAL_MAX_AB_SIZE];
	struct complex rb13[V1_REAL_MAX_AB_SIZE];
	struct complex rb14[V1_REAL_MAX_AB_SIZE];
	struct complex rb15[V1_REAL_MAX_AB_SIZE];
	struct complex rb16[V1_REAL_MAX_AB_SIZE];
	struct complex rb17[V1_REAL_MAX_AB_SIZE];


	sum(a1, a3, ra1, n5);
	difference(a0, a4, ra2, n5);
	sum(a0, a4, ra3, n5);
	difference(a1, a3, ra4, n5);
	difference(ra4, a2, ra5, n5);
	sum_neg_neg(a2, ra4, ra6, n5);
	difference(ra3, a2, ra7, n5);
	sum(ra1, ra2, ra8, n5);
	difference(ra2, ra1, ra9, n5);
	sum(ra1, ra3, ra10, n5);
	sum(ra10, a2, ra11, n5);
	create_term(ra7, ra4, ra12, n5);
	create_term_neg(ra7, ra4, ra13, n5);
	create_term(ra8, ra5, ra14, n5);
	create_term_neg(ra8, ra5, ra15, n5);
	create_term(ra9, ra6, ra16, n5);
	create_term_neg(ra9, ra6, ra17, n5);

	sum(b1, b3, rb1, n5);
	difference(b0, b4, rb2, n5);
	sum(b0, b4, rb3, n5);
	difference(b1, b3, rb4, n5);
	difference(rb4, b2, rb5, n5);
	sum_neg_neg(b2, rb4, rb6, n5);
	difference(rb3, b2, rb7, n5);
	sum(rb1, rb2, rb8, n5);
	difference(rb2, rb1, rb9, n5);
	sum(rb1, rb3, rb10, n5);
	sum(rb10, b2, rb11, n5);
	create_term(rb7, rb4, rb12, n5);
	create_term_neg(rb7, rb4, rb13, n5);
	create_term(rb8, rb5, rb14, n5);
	create_term_neg(rb8, rb5, rb15, n5);
	create_term(rb9, rb6, rb16, n5);
	create_term_neg(rb9, rb6, rb17, n5);


	struct complex p0[V1_REAL_MAX_P_SIZE];
	struct complex p1[V1_REAL_MAX_P_SIZE];
	struct complex p2[V1_REAL_MAX_P_SIZE];
	struct complex p3[V1_REAL_MAX_P_SIZE];
	struct complex p4[V1_REAL_MAX_P_SIZE];
	struct complex p5[V1_REAL_MAX_P_SIZE];
	struct complex p6[V1_REAL_MAX_P_SIZE];
	struct complex p7[V1_REAL_MAX_P_SIZE];
	struct complex p8[V1_REAL_MAX_P_SIZE];
	memset(p0, 0, sizeof(p0));
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));
	memset(p4, 0, sizeof(p4));
	memset(p5, 0, sizeof(p5));
	memset(p6, 0, sizeof(p6));
	memset(p7, 0, sizeof(p7));
	memset(p8, 0, sizeof(p8));


	ka2_real(a0, b0, p0, n5);
	ka2_real(ra11, rb11, p1, n5);
	n1_complex(ra12, rb12, p2, n5);
	phase_change(p2, p3, p_size);
	n1_complex(ra14, rb14, p4, n5);
	phase_change(p4, p5, p_size);
	n1_complex(ra16, rb16, p6, n5);
	phase_change(p6, p7, p_size);
	ka2_real(a4, b4, p8, n5);

	struct complex u1[V1_REAL_MAX_P_SIZE];
	struct complex u2[V1_REAL_MAX_P_SIZE];
	struct complex u3[V1_REAL_MAX_P_SIZE];
	struct complex u4[V1_REAL_MAX_P_SIZE];
	struct complex u5[V1_REAL_MAX_P_SIZE];
	struct complex u6[V1_REAL_MAX_P_SIZE];
	struct complex u7[V1_REAL_MAX_P_SIZE];
	struct complex u8[V1_REAL_MAX_P_SIZE];
	struct complex u9[V1_REAL_MAX_P_SIZE];
	struct complex u10[V1_REAL_MAX_P_SIZE];
	struct complex u11[V1_REAL_MAX_P_SIZE];
	struct complex u12[V1_REAL_MAX_P_SIZE];
	struct complex u13[V1_REAL_MAX_P_SIZE];
	struct complex u14[V1_REAL_MAX_P_SIZE];

	sum_real(p0, p8, u1, p_size);
	difference_real(p1, u1, u2, p_size);
	difference_real(p2, p6, u3, p_size);
	difference_real(u3, p4, u4, p_size);
	difference_real(p2, p4, u5, p_size);
	sum_real(p6, p4, u6, p_size);
	sum_im_neg_im_real(p2, p4, u7, p_size);
	sum_re_im_real(u7, p6, u8, p_size);
	sum_im_im_real(p4, p6, u9, p_size);
	sum_im_im_real(p2, p4, u10, p_size);
	sum_re_neg_im_real(u10, p6, u11, p_size);
	sum_real(u2, u3, u12, p_size);
	sum_real(u1, u4, u13, p_size);
	sum_real(u2, u5, u14, p_size);


	struct complex c1[V1_REAL_MAX_P_SIZE];
	struct complex c2[V1_REAL_MAX_P_SIZE];
	struct complex c3[V1_REAL_MAX_P_SIZE];
	struct complex c4[V1_REAL_MAX_P_SIZE];
	struct complex c5[V1_REAL_MAX_P_SIZE];
	struct complex c6[V1_REAL_MAX_P_SIZE];
	struct complex c7[V1_REAL_MAX_P_SIZE];

	struct complex *c0 = p0;

	sum_real(u12, u8, c1, p_size);
	difference_real(u13, u9, c2, p_size);
	difference_real(u12, u8, c3, p_size);
	sum_real(u1, u6, c4, p_size);
	sum_real(u14, u11, c5, p_size);
	sum_real(u13, u9, c6, p_size);
	difference_real(u14, u11, c7, p_size);

	struct complex *c8 = p8;

	sum_real(c, c0, c, p_size);
	sum_real(c+n5, c1, c+n5, p_size);
	sum_real(c+2*n5, c2, c+2*n5, p_size);
	sum_real(c+3*n5, c3, c+3*n5, p_size);
	sum_real(c+4*n5, c4, c+4*n5, p_size);
	sum_real(c+5*n5, c5, c+5*n5, p_size);
	sum_real(c+6*n5, c6, c+6*n5, p_size);
	sum_real(c+7*n5, c7, c+7*n5, p_size);
	sum_real(c+8*n5, c8, c+8*n5, p_size);


	/*for(int i=0; i<2*n-1; i++){
		c[i].re = ((c[i].re%3)+3)%3;
		c[i].im = ((c[i].im%3)+3)%3;
	}*/

}

void lt_real(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define LT_REAL_MAX_P_SIZE 1519 // 2*n-3

	// removed base case

	int n1 = n-1;

	struct complex *a1 = a;
	struct complex *b1 = b;

	struct complex p[LT_REAL_MAX_P_SIZE];
	memset(p, 0, sizeof(p));

	if(n == 761)
		v1_real(a1, b1, p, n1);

	else // n == 9
		ka2_real(a1, b1, p, n1);


	sum_real(c, p, c, 2*n1-1);

	for(int i=0; i<n1 ; ++i){
		c[i+n1].re += a1[i].re * b[n-1].re + b1[i].re * a[n-1].re;
	}

	c[2*n-2].re += a[n-1].re * b[n-1].re;

}




int main(){


	struct complex *a = malloc(761*sizeof(struct complex));
	struct complex *b = malloc(761*sizeof(struct complex));
	struct complex *c = calloc((2*761-1), sizeof(struct complex));

     int i;
////////////////////////reading a and b from the file inp_a3
	FILE *myFile;
    myFile = fopen("inp761test", "r");
    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (0);}

    for (i = 0; i < 761; i++){
        fscanf(myFile, "%d", &a[i].re );
        a[i].im = 0;
    }
    
    for (i = 0; i < 761; i++){
        fscanf(myFile, "%d", &b[i].re );
        b[i].im = 0;
    }
    fclose(myFile);
    

        unsigned long long   LastCycleCount = _rdtsc();
    	clock_t start = clock ();	     
   		for (i=0;i<99999;i++){
        lt_real(a, b, c, 761);
        memset(c, 0, sizeof(struct complex)*(2*761-1));
    	}
    	lt_real(a, b, c, 761);
	    unsigned long long   EndCycleCount = _rdtsc();
        unsigned long long   CyclesElapsed = EndCycleCount - LastCycleCount;
        CyclesElapsed = CyclesElapsed/100000 ;
   
        double Timelapsed=(clock()-start)/(double) CLOCKS_PER_SEC;
        Timelapsed=Timelapsed/100000;
 /////////////////////////////////////////////write on c_a3 file	
	FILE *outFile;
    outFile = fopen("c_lt_hybrid", "w");
    if (outFile == NULL)
	{
    	printf("Cannot Open File\n");
        exit (0);
    }
    
   for (i = 0; i < 2*761-1; i++)
        fprintf(outFile, "%d ", ((c[i].re)%3+3)%3 );
        
	fclose(outFile);
/////////////////////////////////////////////////////////////printing on console
	for( i=0; i<2*761-1; i++)
		printf("%d ", ((c[i].re)%3+3)%3);
    printf("\n\n\n");
	printf("Cycles: %llu\n", CyclesElapsed);
    printf("Multiplication Time:");
    printf( "%f",Timelapsed);
    printf("sec");
    printf("\n \n ");

	free(a);
	free(b);
	free(c);

	return 0;
}
