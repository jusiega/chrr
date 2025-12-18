#ifndef CH_CHI_POTENSZJALS2_H
#define CH_CHI_POTENSZJALS2_H


//0
double pustaprzestrzeń(double x, const double *pp)
{
    return 0.0;
}

//a(x-x0)^2
double oscylator(double x, const double *pp)//[a, x0]
{
    return pp[0]*(x-pp[1])*(x-pp[1]);
}

//ax
double skośny(double x, const double *pp)//[a]
{
    return pp[0]*x;
}

//polynomial
//ax^6+bx^5+cx^4+dx^3+ex^2+fx^1+g
double wielomianowy(double x, const double *pp)//[g,f,e,f,c,b,a]
{
    //jaki indeks parametru taka potęga iksa
    return pp[6]*x*x*x*x*x*x+pp[5]*x*x*x*x*x+pp[4]*x*x*x*x+pp[3]*x*x*x+pp[2]*x*x+pp[1]*x+pp[0];
}


double woodssaxon(double x, const double *pp)
{
    //double U0=pp[00];
    //double x0=pp[1];
    //double a=pp[2];
    return -pp[0]/(1.+exp((fabs(x)-pp[1])/pp[2]));
}

//woods-saxon step potential
double stepwoodssaxon(double x, const double *pp)
{
    //double U0=pparameters[0];
    //double x0=pparameters[1];
    //double a=pparameters[2];
    //double ref=pparameters[3];
    return -pp[0]/(1.+exp((x-pp[1])/pp[2]))+pp[3];
}

//power-exponential
double wykładniczopotęgowy(double x, const double *pp)
{
    //aż cztery parametry
    //U0=pp[0]
    //x0=pp[1]
    //a=pp[2]
    //n=pp[3]
    return pp[0]*exp(-pow((x-pp[1])/pp[2],pp[3]));
}

#endif
