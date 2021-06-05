static double hbar = 1 ;
static double m = 0.5 ;

double fx(double x, double psi, double D_psi) {
	return D_psi ;
}

double dfx(double x, double psi, double D_psi) {
	double E ;
	E = (TMath::Pi()*TMath::Pi())/4. ;
	return ((2.*m*E)/(hbar*hbar))*psi ;
}


void prob6()
{
	double h = 0.1 ;
	int N = 10 ;
	double psi, x, D_psi, E ;
	double fx1, fx2, fx3, fx4, dfx1, dfx2, dfx3, dfx4 ;
	double PSI[N+1], T[N+1], M_T[N+1] ;
	psi = 1 ;
	D_psi = 0 ;

	PSI[0] = 1 ;
	T[0] = 0 ;
	M_T[0] = 0 ;


	for (int n = 0 ; n<N ; n++ ) {

		x = (n+1)*h ;
		

		fx1 = fx(x,psi,D_psi) ;
		dfx1 = dfx(x,psi,D_psi) ;

		fx2 = fx(x+0.5*h,psi+0.5*h*fx1,D_psi+0.5*h*dfx1) ;
		dfx2 = dfx(x+0.5*h,psi+0.5*h*fx1,D_psi+0.5*h*dfx1);

		fx3 = fx(x+0.5*h,psi+0.5*h*fx2,D_psi+0.5*h*dfx2) ;
		dfx3 = dfx(x+0.5*h,psi+0.5*h*fx2,D_psi+0.5*h*dfx2) ;

		fx4 = fx(x+h,psi+h*fx3,D_psi+h*dfx3) ;
		dfx4 = dfx(x+h,psi+h*fx3,D_psi+h*dfx3) ;
		
		psi += (h/6)*(fx1+2*fx2+2*fx3+fx4) ;
		D_psi += -(h/6)*(dfx1+2*dfx2+2*dfx3+dfx4) ;
		
		PSI[n+1] = psi ;
		T[n+1] = x ;
		M_T[n+1] = -x ;
		printf("%f\t%f\t%f\t%d\n",psi,D_psi,x,n) ;
	}

	TGraph *p = new TGraph(N+1,T,PSI) ;
	p -> Draw("ACP");	

	p -> GetXaxis() -> SetLimits(-1.1,1.1) ;
	TGraph *g = new TGraph(N+1,M_T,PSI) ;
	g -> Draw("same");

	TF1 *k = new TF1("k","cos(TMath::Pi()/2*x)",-1,1) ;
	k -> Draw("same");

}
