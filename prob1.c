static double Cd = 0.25 ;
static double A = 0.7 ;
static double air_density = 1.293 ;
static double g = 9.81 ;
static double  m = 65 ; //mass
static double h = 0.01 ; //step size
static double H = 1 ; // step size
static double beta = 0.5*Cd*A*air_density ;

// ma = mg - 1/2 (air_density * A * Cd * v*v)
// a = dv/dt , v = dx/dt

double fx(double x, double v) {
	return v ; 
}

double fv(double x, double v) {
	return g - 0.5*((air_density * A * Cd * v * v) / m) ;
}

double fax(double x, double v) {
	return fx(x+(2./3.)*h,v+(2./3.)*h*fx(x,v)) ;
}

double fav(double x, double v) {
	return fv(x+(2./3.)*h,v+(2./3.)*h*fv(x,v)) ;
}

void prob1_1()
{
	double v, t, x, N, n ;
	N = 1./h ;
	for (double end = 1 ; end <= 10 ; end++ ) {
	x = 0 ; //reset x
	v = 0 ;	//reset v
		for (n = 0 ; n < 100*end ; n++ ) { //reset and count n
		x += h*(0.25*fx(x,v) + 0.75*fax(x,v)) ;
		v += h*(0.25*fv(x,v) + 0.75*fav(x,v)) ;
		}

	printf("x : %f,\t v : %f,\t n : %f\n",x,v,n/100) ; //time t = n*h
	}	

}

void prob1_2() {
	double x = 0. , v = 0. ;
	double f1, f2, f3, f4 ;
	double F1, F2, F3, F4 ;
	int N = 10 ;
	double t_[N+1], v_[N+1] ;

	t_[0] = 0 ;
	v_[0] = 0 ;

	for (int j = 0 ; j < 10 ; j++) {
		for ( int i = 0 ; i < N ; i++) {
			
			f1 = fx(x,v) ;
			F1 = fv(x,v) ;
			f2 = fx(x+0.5*f1,v+0.5*F1) ;
			F2 = fv(x+0.5*f1,v+0.5*F1) ;
			f3 = fx(x+0.5*f2,v+0.5*F2) ;
			F3 = fv(x+0.5*f2,v+0.5*F2) ;
			f4 = fx(x+f3,v+F3) ;
			F4 = fv(x+f3,v+F3) ;
		
			x += H/6.*(f1 + f2*2 + f3*2 + f4) ;
			v += H/6.*(F1 + F2*2 + F3*2 + F4) ;

			}
		t_[j+1] = (j+1) ;
		v_[j+1] = v ;
		
		printf("%f\t%f\t%d\n",x,v,(j+1)*10) ;
		}

	TGraph *p = new TGraph(N+1,t_,v_) ;
	p -> Draw("APC");

	TF1 *P = new TF1("P", "sqrt((9.81*65)/0.1131375)*TMath::TanH(sqrt(65*0.1131375*9.81)*x) ",0,100) ;
	P -> Draw("same");

}
