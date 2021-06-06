static double la = 0.693/5018 ; // ramda = log2/half_life_time
static double lb = 0.693/138376 ;

void prob5()
{
	double Na = 100000, Nb = 0, Nc = 0 ;
	double h = 100 ;
	double NA[50000],NB[50000],NC[50000], T[50000] ;
	int t ;


	for (int n = 0 ; n<(5*10e+5/100) ; n++ ) {
		t = n*h ;
		Na += -h*la*Na ;
		Nb += -h*(-la*Na + lb*Nb) ;
		Nc += h*(lb*Nb) ;
		NA[n] = (Na) ;
		NB[n] = (Nb) ;
		NC[n] = (Nc) ;
		T[n] = t ;
		printf("Na : %f,\t Nb : %f,\t Nc : %f,\t Na+Nb+Nc : %f,\t time t : %d\n",Na,Nb,Nc,Na+Nb+Nc,t) ;
	}

	TGraph *p = new TGraph(5000,T,NA) ;
	p -> Draw("AP") ;
	TGraph *g = new TGraph(5000,T,NB) ;
	g -> Draw("same") ;
	TGraph *k = new TGraph(5000,T,NC) ;
	k -> Draw("same") ;

	g -> SetLineColor(4);
	k -> SetLineColor(7);

}

