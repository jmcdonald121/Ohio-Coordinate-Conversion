void south2north(float *x,float *y)
	{
	double pi=3.141592654;
	double a=20925832.16;        /*radius of the earth r=6378206.4 m*/
	double e2=0.006768658;
	double e;
	double phi0n=0.692313936;
	double phi1n=0.705694794;
	double phi2n=0.727802298;
	double lamda0n=1.439896633;
	double phi0s=0.663225115;
	double phi1s=0.676024196;
	double phi2s=0.698713477;
	double lamda0s=1.439896633;
	double n,f,rho,rho0,theta;
	double phi;
	double lamda;
	double phiold;
	double m,m1,m2;
	double t,t0,t1,t2,tnew;
	/*****************************/
	/* input the x-y coordinates */
	/*****************************/
	*x=2000000-*x;
	/**************************************/
	/* compute the latitude and longitude */
	/**************************************/
	e=sqrt(e2);
	m=cos(phi)/sqrt(1.0-e2*pow(sin(phi),2.0));
	m1=cos(phi1s)/sqrt(1.0-e2*pow(sin(phi1s),2.0));
	m2=cos(phi2s)/sqrt(1.0-e2*pow(sin(phi2s),2.0));
	t=tan(pi/4.0-phi/2.0)/pow((1.0-e*sin(phi))/(1.0+e*sin(phi)),e/2.0);
	t0=tan(pi/4.0-phi0s/2.0)/pow((1.0-e*sin(phi0s))/(1.0+e*sin(phi0s)),e/2.0);
	t1=tan(pi/4.0-phi1s/2.0)/pow((1.0-e*sin(phi1s))/(1.0+e*sin(phi1s)),e/2.0);
	t2=tan(pi/4.0-phi2s/2.0)/pow((1.0-e*sin(phi2s))/(1.0+e*sin(phi2s)),e/2.0);
	n=(log(m1)-log(m2))/(log(t1)-log(t2));
	f=m1/(n*pow(t1,n));
	rho0=a*f*pow(t0,n);
	rho=sqrt(pow((double)*x,2.0)+pow((rho0-(double)*y),2.0));
	theta=atan2((double)*x,(rho0-(double)*y));
	tnew=pow(rho/(a*f),1.0/n);
	phi=pi/2.0-2.0*atan(tnew);
	phiold=0.0;
	while(fabs(phi-phiold)>=0.00000000000001)
		{
		phiold=phi;
		phi=pi/2.0-2.0*atan(tnew*pow((1.0-e*sin(phi))/(1.0+e*sin(phi)),e/2.0));
		}
	lamda=theta/n+lamda0s;
	/**************************************************/
	/* compute the south zone state plane coordinates */
	/**************************************************/
	e=sqrt(e2);
	m=cos(phi)/sqrt(1.0-e2*pow(sin(phi),2.0));
	m1=cos(phi1n)/sqrt(1.0-e2*pow(sin(phi1n),2.0));
	m2=cos(phi2n)/sqrt(1.0-e2*pow(sin(phi2n),2.0));
	t=tan(pi/4.0-phi/2.0)/pow((1.0-e*sin(phi))/(1.0+e*sin(phi)),e/2.0);
	t0=tan(pi/4.0-phi0n/2.0)/pow((1.0-e*sin(phi0n))/(1.0+e*sin(phi0n)),e/2.0);
	t1=tan(pi/4.0-phi1n/2.0)/pow((1.0-e*sin(phi1n))/(1.0+e*sin(phi1n)),e/2.0);
	t2=tan(pi/4.0-phi2n/2.0)/pow((1.0-e*sin(phi2n))/(1.0+e*sin(phi2n)),e/2.0);
	n=(log(m1)-log(m2))/(log(t1)-log(t2));
	f=m1/(n*pow(t1,n));
	rho0=a*f*pow(t0,n);
	rho=a*f*pow(t,n);
	theta=n*(lamda-lamda0n);
	*x=rho*sin(theta);
	*x=2000000-*x;
	*y=rho0-rho*cos(theta);
	}
