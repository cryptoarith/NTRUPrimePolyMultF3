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



void sb_comba_complex(struct complex a[], struct complex b[], struct complex c[], int n){
	
	int re;
	int im;

	for(int i=0; i<n; ++i){
		re = im = 0;
		for(int j=0; j<=i; ++j){
			re += a[j].re * b[i-j].re - a[j].im * b[i-j].im;
			im += a[j].re * b[i-j].im + a[j].im * b[i-j].re;
		}
		c[i].re = re;
		c[i].im = im;
	}

	for(int i=n; i<2*n-1; ++i){
		re = im = 0;
		for(int j=i-n+1; j<n; ++j){
			re += a[j].re * b[i-j].re - a[j].im * b[i-j].im;
			im += a[j].re * b[i-j].im + a[j].im * b[i-j].re;
		}
		c[i].re = re;
		c[i].im = im;
	}

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void a3_complex(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define A3_COMPLEX_MAX_AB_SIZE 4 // n/3
	#define A3_COMPLEX_MAX_P_SIZE 7 // 2*n/3-1

	// removed base case

	/* Assumed n = 3^k for some non-negative integer k */
	int n3 = n/3;
	int p_size = 2*n3-1;


	struct complex *a0 = a;
	struct complex *a1 = a+n3;
	struct complex *a2 = a+2*n3;
	struct complex *b0 = b;
	struct complex *b1 = b+n3;
	struct complex *b2 = b+2*n3;

	struct complex ra1[A3_COMPLEX_MAX_AB_SIZE];
	struct complex ra2[A3_COMPLEX_MAX_AB_SIZE];
	struct complex ra3[A3_COMPLEX_MAX_AB_SIZE];
	struct complex ra4[A3_COMPLEX_MAX_AB_SIZE];
	struct complex ra5[A3_COMPLEX_MAX_AB_SIZE];
	struct complex rb1[A3_COMPLEX_MAX_AB_SIZE];
	struct complex rb2[A3_COMPLEX_MAX_AB_SIZE];
	struct complex rb3[A3_COMPLEX_MAX_AB_SIZE];
	struct complex rb4[A3_COMPLEX_MAX_AB_SIZE];
	struct complex rb5[A3_COMPLEX_MAX_AB_SIZE];

	difference(a0, a2, ra1, n3);
	sum(a0, a2, ra2, n3);
	sum(ra2, a1, ra3, n3);
	create_term(ra1, a1, ra4, n3);
	create_term_neg(ra1, a1, ra5, n3);

	difference(b0, b2, rb1, n3);
	sum(b0, b2, rb2, n3);
	sum(rb2, b1, rb3, n3);
	create_term(rb1, b1, rb4, n3);
	create_term_neg(rb1, b1, rb5, n3);


	struct complex p0[A3_COMPLEX_MAX_P_SIZE];
	struct complex p1[A3_COMPLEX_MAX_P_SIZE];
	struct complex p2[A3_COMPLEX_MAX_P_SIZE];
	struct complex p3[A3_COMPLEX_MAX_P_SIZE];
	struct complex p4[A3_COMPLEX_MAX_P_SIZE];
	memset(p0, 0, sizeof(p0));
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));
	memset(p4, 0, sizeof(p4));


	sb_comba_complex(a0, b0, p0, n3);
	sb_comba_complex(ra3, rb3, p1, n3);
	sb_comba_complex(ra4, rb4, p2, n3);
	sb_comba_complex(ra5, rb5, p3, n3);
	sb_comba_complex(a2, b2, p4, n3);

	struct complex u1[A3_COMPLEX_MAX_P_SIZE];
	struct complex u2[A3_COMPLEX_MAX_P_SIZE];
	struct complex u3[A3_COMPLEX_MAX_P_SIZE];
	struct complex u4[A3_COMPLEX_MAX_P_SIZE];
	struct complex u5[A3_COMPLEX_MAX_P_SIZE];

	sum(p0, p4, u1, p_size);
	sum(p2, p3, u2, p_size);
	difference(p2, p3, u3, p_size);
	difference(u2, u1, u4, p_size);
	difference(u4, p1, u5, p_size);


	struct complex c1[A3_COMPLEX_MAX_P_SIZE];
	struct complex c2[A3_COMPLEX_MAX_P_SIZE];
	struct complex c3[A3_COMPLEX_MAX_P_SIZE];

	struct complex *c0 = p0;

	create_term_neg(u5, u3, c1, p_size);
	sum(u1, u2, c2, p_size);
	create_term(u5, u3, c3, p_size);

	struct complex *c4 = p4;

	sum(c, c0, c, p_size);
	sum(c+n3, c1, c+n3, p_size);
	sum(c+2*n3, c2, c+2*n3, p_size);
	sum(c+3*n3, c3, c+3*n3, p_size);
	sum(c+4*n3, c4, c+4*n3, p_size);


	/*for(int i=0; i<2*n-1; ++i){
		c[i].re = ((c[i].re%3)+3)%3;
		c[i].im = ((c[i].im%3)+3)%3;
	}*/


}

void n1_complex(struct complex a[], struct complex b[], struct complex c[], int n){

	#define N1_COMPLEX_MAX_AB_SIZE 48 // n/4
	#define N1_COMPLEX_MAX_P_SIZE 95 // 2*n/4-1

	/* base case */
	if(n == 12){
		sb_comba_complex(a, b, c, n);   ///////////////////////////////////////burasý normalde a3_complex idi deðiþiyor
		return;
	}

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

	n1_complex(a0, b0, p0, n4);
	n1_complex(ra3, rb3, p1, n4);
	n1_complex(ra9, rb9, p2, n4);
	n1_complex(ra11, rb11, p3, n4);	
	n1_complex(ra10, rb10, p4, n4);
	n1_complex(ra12, rb12, p5, n4);
	n1_complex(a3, b3, p6, n4);

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

void ka2_real(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define KA2_REAL_MAX_AB_SIZE 96 // n/2
	#define KA2_REAL_MAX_P_SIZE 191 // n-1

	/* base case */
	if(n == 48){      ///////////////////////////////////////////////////////////////////////////////////////////////////////////bu satýr deðiþiyor
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

	ka2_real(a0, b0, p1, n2);
	ka2_real(sum_a0_a1, sum_b0_b1, p2, n2);
	ka2_real(a1, b1, p3, n2);


	sum_real(c, p1, c, n-1);
	difference_real(c+n2, p3, c+n2, n-1);

	for(int i=n+n2-2; i>=0 ; --i){
		c[i+n2].re -= c[i].re;
	}

	sum_real(c+n2, p2, c+n2, n-1);

}


void n1_real(struct complex a[], struct complex b[], struct complex c[], int n){

	#define N1_REAL_MAX_AB_SIZE 192 // n/4
	#define N1_REAL_MAX_P_SIZE 383 // 2*n/4-1

	/* base case */
	if(n == 1){
		multiplication_real(a, b, c, n);
		return;
	}

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

	struct complex ra0[N1_REAL_MAX_AB_SIZE];
	struct complex ra1[N1_REAL_MAX_AB_SIZE];
	struct complex ra2[N1_REAL_MAX_AB_SIZE];
	struct complex ra3[N1_REAL_MAX_AB_SIZE];
	struct complex ra4[N1_REAL_MAX_AB_SIZE];
	struct complex ra5[N1_REAL_MAX_AB_SIZE];
	struct complex ra6[N1_REAL_MAX_AB_SIZE];
	struct complex ra7[N1_REAL_MAX_AB_SIZE];
	struct complex ra8[N1_REAL_MAX_AB_SIZE];
	struct complex ra9[N1_REAL_MAX_AB_SIZE];
	struct complex ra10[N1_REAL_MAX_AB_SIZE];
	struct complex ra11[N1_REAL_MAX_AB_SIZE];
	struct complex ra12[N1_REAL_MAX_AB_SIZE];

	struct complex rb0[N1_REAL_MAX_AB_SIZE];
	struct complex rb1[N1_REAL_MAX_AB_SIZE];
	struct complex rb2[N1_REAL_MAX_AB_SIZE];
	struct complex rb3[N1_REAL_MAX_AB_SIZE];
	struct complex rb4[N1_REAL_MAX_AB_SIZE];
	struct complex rb5[N1_REAL_MAX_AB_SIZE];
	struct complex rb6[N1_REAL_MAX_AB_SIZE];
	struct complex rb7[N1_REAL_MAX_AB_SIZE];
	struct complex rb8[N1_REAL_MAX_AB_SIZE];
	struct complex rb9[N1_REAL_MAX_AB_SIZE];
	struct complex rb10[N1_REAL_MAX_AB_SIZE];
	struct complex rb11[N1_REAL_MAX_AB_SIZE];
	struct complex rb12[N1_REAL_MAX_AB_SIZE];

	
	difference(a0, a2, ra0, n4);
	difference(a1, a3, ra1, n4);
	sum(a1, a3, ra2, n4);
	sum(a0, ra2, ra3, n4);
	difference(ra1, a2, ra4, n4);
	difference(a0, ra2, ra5, n4);
	sum_neg_neg(a2, ra1, ra6, n4);
	create_term(ra0, ra1, ra7, n4);
	create_term(ra3, ra4, ra8, n4);
	create_term(ra5, ra6, ra9, n4);
	create_term_neg(ra0, ra1, ra10, n4);
	create_term_neg(ra3, ra4, ra11, n4);
	create_term_neg(ra5, ra6, ra12, n4);

	difference(b0, b2, rb0, n4);
	difference(b1, b3, rb1, n4);
	sum(b1, b3, rb2, n4);
	sum(b0, rb2, rb3, n4);
	difference(rb1, b2, rb4, n4);
	difference(b0, rb2, rb5, n4);
	sum_neg_neg(b2, rb1, rb6, n4);
	create_term(rb0, rb1, rb7, n4);
	create_term(rb3, rb4, rb8, n4);
	create_term(rb5, rb6, rb9, n4);
	create_term_neg(rb0, rb1, rb10, n4);
	create_term_neg(rb3, rb4, rb11, n4);
	create_term_neg(rb5, rb6, rb12, n4);


	struct complex p0[N1_REAL_MAX_P_SIZE];
	struct complex p1[N1_REAL_MAX_P_SIZE];
	struct complex p2[N1_REAL_MAX_P_SIZE];
	struct complex p3[N1_REAL_MAX_P_SIZE];
	struct complex p4[N1_REAL_MAX_P_SIZE];
	struct complex p5[N1_REAL_MAX_P_SIZE];
	struct complex p6[N1_REAL_MAX_P_SIZE];
	memset(p0, 0, sizeof(p0));
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));
	memset(p4, 0, sizeof(p4));
	memset(p5, 0, sizeof(p5));
	memset(p6, 0, sizeof(p6));


	n1_complex(ra7, rb7, p0, n4);
	phase_change(p0, p1, p_size);
	n1_complex(ra8, rb8, p2, n4);
	phase_change(p2, p3, p_size);
	n1_complex(ra9, rb9, p4, n4);
	phase_change(p4, p5, p_size);
	ka2_real(a3, b3, p6, n4);


	struct complex u0[N1_REAL_MAX_P_SIZE];
	struct complex u1[N1_REAL_MAX_P_SIZE];
	struct complex u2[N1_REAL_MAX_P_SIZE];
	struct complex u3[N1_REAL_MAX_P_SIZE];
	struct complex u4[N1_REAL_MAX_P_SIZE];
	struct complex u5[N1_REAL_MAX_P_SIZE];


	sum_im_im_real(p2, p4, u0, p_size);
	sum_im_neg_im_real(p4, p2, u1, p_size);
	sum_real(p2, p4, u2, p_size);
	difference_real(p2, p4, u3, p_size);
	difference_real(p6, p0, u4, p_size);
	difference_real(u4, u0, u5, p_size);


	struct complex c0[N1_REAL_MAX_P_SIZE];
	struct complex c1[N1_REAL_MAX_P_SIZE];
	struct complex c2[N1_REAL_MAX_P_SIZE];
	struct complex c3[N1_REAL_MAX_P_SIZE];
	struct complex c4[N1_REAL_MAX_P_SIZE];
	struct complex c5[N1_REAL_MAX_P_SIZE];


	sum_real(u5, u2, c0, p_size);
	sum_re_neg_im_real(u3, p0, c1, p_size);
	sum_real(u0, p6, c2, p_size);
	sum_real(u1, u3, c3, p_size);
	difference_real(u5, u2, c4, p_size);
	sum_re_neg_im_real(u1, p0, c5, p_size);

	struct complex *c6 = p6;

	sum_real(c, c0, c, p_size);
	sum_real(c+n4, c1, c+n4, p_size);
	sum_real(c+2*n4, c2, c+2*n4, p_size);
	sum_real(c+3*n4, c3, c+3*n4, p_size);
	sum_real(c+4*n4, c4, c+4*n4, p_size);
	sum_real(c+5*n4, c5, c+5*n4, p_size);
	sum_real(c+6*n4, c6, c+6*n4, p_size);


	/*for(int i=0; i<2*n-1; i++){
		c[i].re = ((c[i].re%3)+3)%3;
		c[i].im = ((c[i].im%3)+3)%3;
	}*/

}




int main(){


	struct complex *a = malloc(768*sizeof(struct complex));
	struct complex *b = malloc(768*sizeof(struct complex));
	struct complex *c = calloc((2*768-1), sizeof(struct complex));

     int i;
////////////////////////reading a and b from the file inp_a3
	FILE *myFile;
    myFile = fopen("inp768test", "r");
    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (0);}

    for (i = 0; i < 768; i++){
        fscanf(myFile, "%d", &a[i].re );
        a[i].im = 0;
    }
    
    for (i = 0; i < 768; i++){
        fscanf(myFile, "%d", &b[i].re );
        b[i].im = 0;
    }
    fclose(myFile);
    

        unsigned long long   LastCycleCount = _rdtsc();
    	clock_t start = clock ();	     
   	for (i=0;i<99999;i++){
        n1_real(a, b, c, 768);
        memset(c, 0, sizeof(struct complex)*(2*768-1));
    	}
    	n1_real(a, b, c, 768);
	    unsigned long long   EndCycleCount = _rdtsc();
        unsigned long long   CyclesElapsed = EndCycleCount - LastCycleCount;
        CyclesElapsed = CyclesElapsed/100000;
   
        double Timelapsed=(clock()-start)/(double) CLOCKS_PER_SEC;
        Timelapsed=Timelapsed/100000;
 /////////////////////////////////////////////write on c_a3 file	
	FILE *outFile;
    outFile = fopen("c_n1_hybrid2", "w");
    if (outFile == NULL)
	{
    	printf("Cannot Open File\n");
        exit (0);
    }
    
   for (i = 0; i < 2*768-1; i++)
        fprintf(outFile, "%d ", ((c[i].re)%3+3)%3 );
        
	fclose(outFile);
/////////////////////////////////////////////////////////////printing on console
	for( i=0; i<2*768-1; i++)
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
