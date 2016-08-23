#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>
using namespace std;
int main ()
{ char command[20],*filename;
  char *station,*roast,*estreg,*estkrig,*moranreg,*output;
  ifstream obsfile;
  int index;

  void dmatern(double,double *,double,double *),d2matern(double,double *,double,double *),transpose(double *,double *,int,int),prod(double *,double *,double *,int,int,int),sweep(double *,int,int),inverse(double *,double *,int),ls(double *,int),gls(double *,double *,double *,int,int,double *),betacondtheta(double *,double *,double *,double *,double,int,int,double*),moranI(double *,double *,int,double*),neighbour(double *,int,double *),colrow(char*,int*,int*);
  int loglikelihood(double *,double *,double *,double *,double,int,int,int,double *),thetaest(double *,double *,double *,double *,double,int,int),mle(double *,double *,double *,double,int,int,double *,double *),preddata(char *,char *,char *,char *,char *,char *),estimation(char *),cv(char *,char *,char *,char *,char *);
  double spheric(double *,double *),mgamma(double),pnorm(double),qnorm(double),besselJ(double,double),besselY(double,double),besselI(double,double),besselK(double,double),matern(double,double *,double),trace(double *,int),detlog(double *,int),regression(double *,double *,double *,int,int,double *,double*),krig(double *,double *,double *,double *,double *,double *,double *,double,int,int,double *);
 
  //estimation
 station="obsnew02072011.txt";   //observations
 index=estimation(station);
 cout <<"Estimation index is:"<<index<<endl;
 
 //prediction
 station="obs02072011.txt";
 roast="roast.txt";
 estreg="estimate.reg";
 estkrig="estimate.krig";
 moranreg="moranI.reg";
 output="prediction.txt";
 index=preddata(station,roast,estreg,estkrig,moranreg,output);
 cout <<"Pred index is "<<index<<endl;

 //cross validation
 /*output="cv.txt";
 index=cv(station,estreg,estkrig,moranreg,output);
 cout <<"CV index is "<<index<<endl;
 */
 return 0;
}


/* spheric distance, input longitude and latitude of two points and output their spheric distance */
double spheric(double *point1,double *point2)
{double p1[3],p2[3],S,R=6378.1,PI=3.1415926535897932;
 p1[0]=cos(point1[1]*PI/180)*cos(point1[0]*PI/180);
 p1[1]=cos(point1[1]*PI/180)*sin(point1[0]*PI/180);
 p1[2]=sin(point1[1]*PI/180);
 p2[0]=cos(point2[1]*PI/180)*cos(point2[0]*PI/180);
 p2[1]=cos(point2[1]*PI/180)*sin(point2[0]*PI/180);
 p2[2]=sin(point2[1]*PI/180);
 S=2*R*asin(0.5*sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2])));
 return S;
}

/* logrithm of gamma function */
double gamma(double alpha)
{int k;
 double S=0.0,a=alpha,S1=0.0,f,x=30.0;
 while (a>1){
 S=S+log(a-1);
 a=a-1;}

 if (a<1) {
 f=a*log(x)-x-log(a);
 S1=S1+exp(f);
 a=a+1;
 while((exp(f)>pow(10,-12.0))||(a<x)){
 f=f+log(x)-log(a);
 S1=S1+exp(f);
 a=a+1;
 }
 S=S+log(S1);
 }
 return S;
}

/*
FInd the p-value of standard gamma distribution
*/

double pgamma(double x,double alpha)
{int k,k0;
 double S=0.0,f1,f; 
 if (x<=0) return 0.0;
 f=alpha*log(x)-x-log(alpha)-gamma(alpha);
 S=exp(f);
 if (x>alpha)
 {k=k0=(int)floor(x-alpha+1);
  f=(alpha+k)*log(x)-x-gamma(alpha+1+k);
  S=S+exp(f);
  while(fabs(exp(f)/S)>=1e-12)
      { k++;
        f=f+log(x)-log(alpha+k);
        S=S+exp(f);
      }
   k=k0-1;
   if (k>=1)
   {f=(alpha+k)*log(x)-x-gamma(alpha+1+k);
    S=S+exp(f);
   }

   while((fabs(exp(f)/S)>=1e-12)&&(k>1))
     {k--;
      f=(alpha+k)*log(x)-x-gamma(alpha+1+k);
      S=S+exp(f); 
     }
 }
 if (x<=alpha)
   {k=0;f=0.0;
   while(fabs(exp(f)/S)>=1e-12)
     {k++;
      f=(alpha+k)*log(x)-x-gamma(alpha+1+k);
      S=S+exp(f);
     }
   }
 return S;
}

/*
Find the q-value of standard gamma distribution
*/
double qgamma(double p,double alpha)
{double x1,x2,x,f,f1;
 x1=0.0;x2=alpha;
 f=pgamma(x2,alpha);
 while (f<p) 
   {x1=x2;x2=2*x1;f=pgamma(x2,alpha);
   }
 x=(x1+x2)/2;
 while((x2-x1)/x2>0.01)
   {f=pgamma(x,alpha);
   if (f<p) x1=x; else x2=x; x=(x1+x2)/2;
   }
 x=(x1+x2)/2;f=pgamma(x,alpha)-p;
 while(fabs(f/x)>=1e-12)
   {f1=exp(-x+(alpha-1)*log(x)-gamma(alpha));
   f=pgamma(x,alpha)-p;
   x=x-f/f1;
   }
 return x;
}

/*
P-value of standard normal distribution
*/

double intgral(int k,double b)
{int i;
 double s=0.0,f,a=1/(sqrt(b+1)+sqrt(b));
 if (b==0.0) return (1/(k+0.5));
 for(i=0;i<=k;i++)
   {f=(k+i+1)*log(a)+(k-i)*log(2.0)+0.5*(k-i)*log(b)+gamma(k+1.0)-gamma(i+1.0)-gamma(k-i+1.0);
    s=s+2*exp(f)/(k+i+1.0);
   }
 return s;
}

double pnorm(double x)
{int k,k0;
 double s,f,f0,a,a1,b,PI=3.1415926535897932;
 if (x<=-33) {s=0.0;return s;}
 if (x>=33) {s=1.0;return s;}
 if (x==0) {s=0.5;return s;}
 k=k0=(int)floor((x*x-1)/2)+1;
 a=-x*x/2-0.5*log(2*PI);
 f0=a+2*k*log(fabs(x))-gamma(2*k+2.0)+k*log(2.0)+gamma(k+1.0);
 if (x==0) return 0.5;
 if(x>-4.0)  
 {f=x*exp(f0);s=0.5+f;
 while(fabs(f/s)>1e-10)
   {k++;f0=a+2*k*log(fabs(x))-gamma(2*k+2.0)+k*log(2.0)+gamma(k+1.0);
   f=x*exp(f0);s=s+f;
   }
 k=k0-1;if(k<0) return s;
 f0=a+2*k*log(fabs(x))-gamma(2*k+2.0)+k*log(2.0)+gamma(k+1.0);
 f=x*exp(f0);s=s+f;

 while((fabs(f/s)>1e-10)&&(k>0))
   {k--;f0=a+2*k*log(fabs(x))-gamma(2*k+2.0)+k*log(2.0)+gamma(k+1.0);
   f=x*exp(f0);s=s+f;
   }
 }
 else
   {a=x*x/2;k=k0=0;s=0.0;b=a;a1=exp(-a)/(2*sqrt(PI)); 
    f=f0=intgral(k0,b)*a1;
    while(fabs(f0/f)>1e-12)
      {k0++;f0=a1*intgral(k0,b)/exp(gamma(k0+1.0));
       if ((k0/2.0)==floor(k0/2.0)) f=f+f0;
       else f=f-f0;}
    s=s+f;
    while(fabs(f/s)>1e-12)
      {k++;b=a+k;k0=0;f=f0=intgral(k0,b)*a1*exp(-1.0*k);
        while(fabs(f0/f)>1e-12)
         {k0++;f0=a1*exp(-1.0*k)*intgral(k0,b)/exp(gamma(k0+1.0));
          if ((k0/2.0)==floor(k0/2.0)) f=f+f0;
          else f=f-f0;
	 }
       s=s+f;      
      }
   }
 return s;
}

/*
Input: alpha,
Output: standard normal q-value. PHI^{-1)(alpha).
*/

double qnorm(double alpha)
{int i;
 double x0,y0,k,x1,x,PI=3.1415926535897932;
 if (alpha==0.5) return 0.0;
 x1=x0=0.0;
 if(alpha>0.5)
   {x1=0.0;
    while((pnorm(x1)<alpha))
       x1=x1+1;  x0=x1-1;
      }
 else 
   {x0=0.0;
    while((pnorm(x0)>alpha))
      x0=x0-1;  x1=x0+1;
       } 
 x=(x1+x0)/2;
 while((x1-x0)>0.01)
   {if (pnorm(x)<alpha) x0=x; else x1=x;x=(x1+x0)/2;}
 x1=x0=x; 
 k=exp(-x0*x0/2)/sqrt(2*PI);
 y0=pnorm(x0)-alpha;
 x0=x0-y0/k;
 while(fabs((x1-x0)/x0)>1e-12)
   {x1=x0;x0=x0-y0/k;
    y0=pnorm(x0)-alpha;
    k=exp(-x0*x0/2)/sqrt(2*PI);
    }
 return x0;
}

/*
the density of p-distribution
*/

double dt(double x,int n)
{double a,b,S,PI=3.1415926535897932;
 a=gamma((n+1)/2.0)-gamma(n/2.0)-log(n*PI)/2;
 b=-(n+1)/2.0*log(1+x*x/n);
 S=exp(a+b);
 return S;
}

/*
Input: x, value;r, df.
Output: p-value of t.
for(i=0;i<(r-4);i++)*/
  /* while((fabs(f)>0.000001)&&(i<(r-3)))*/

double pt(double x,int n)
{int i,k;
 double S=0.0,t,theta,f=1.0,z,z2,y,t_x,a,d,p_z,PI=3.1415926535897932;
 if (x==0) {S=0.5; return S;}
 if (x<0) t=x; else t=-x; 
 theta=atan(t/sqrt(n*1.0)); 
 if (n==1) 
    {S=0.5+theta/PI;if (x<0) return S; else return 1-S;}

 if ((n<=20)&&(t>=-4))
   {if ((n%2)==1)
     {S=0.5+theta/PI; f=sin(theta)*cos(theta)/PI;
      S=S+f;
      for(k=1;k<=(n-3)/2;k++)
        {f=f*(2*k)/(2*k+1)*cos(theta)*cos(theta);
         S=S+f;
	}  
     }
   else 
     { f=sin(theta)/2; S=n/(2*sqrt(n+t*t)*(sqrt(n+t*t)-t));
       for(k=1;k<=(n-2)/2;k++)
         {f=f*(2*k-1)*cos(theta)*cos(theta)/(2*k);
          S=S+f;
	 }
     }
  if (x<0) return S; else return 1-S;
   }
 
 if ((n>20)&&(t>=-4))
   {a=n-0.5; d=48*a*a; z=sqrt(a*log(1+t*t/n)); z2=z*z;
    t_x=z;
    p_z=(z2+3)*z/d;
    t_x=t_x+p_z;
    p_z=(((4*z2+33)*z2+240)*z2+855)*z/(10*d*d);
    t_x=t_x-p_z;
    p_z=(((((64*z2+788)*z2+9801)*z2+89775)*z2+543375)*z2+1788885)*z/(210*d*d*d);
    t_x=t_x+p_z;
    S=pnorm(-t_x);
    if (x<0) return S; else return 1-S;
   }

   y=1/sqrt(1+t*t/n); 
    a=gamma((n+1.0)/2)-log(n*PI)/2-gamma(n/2.0)-(n+1.0)/2*log(1+t*t/n)+log(sqrt(n/(y*y)));
    S=1.0/n; f=y*y/2;k=1;S=S+f/(n+2*k);
    while((f/((n+2*k)*S))>1e-9)
      {k++; 
       f=f*y*y*(2*k-1)/(2*k);  S=S+f/(n+2*k);
       }
    S=exp(a)*S;
   if (x<0) return S; else return 1-S;
   
}

/*
Input: level alpha, degree of freedom, r;
Output: qnantile value
*/

double qt(double alpha0,int n)
{double k0,k,x0,x1,x2,y0,alpha,PI=3.1415926535897932;
 int i;
 if (alpha0==0.5) return 0.0;
 x0=0.0;
 if (alpha0<0.5) alpha=alpha0; else alpha=1-alpha0;
 x1=0.0;x2=-1.0; y0=pt(x2,n);
 while(y0>alpha)
   {x1=x2;x2=x1-1;y0=pt(x2,n);
    } 
 
 x0=(x1+x2)/2;
 y0=log(pt(x0,n))-log(alpha);
 k0=exp(gamma((n+1.0)/2)-gamma(n/2.0)-0.5*log(n*PI));
 k=k0/pt(x0,n);
 while(fabs(y0/x0)>1e-8)
  {x0=x0-y0/k;
   y0=log(pt(x0,n))-log(alpha); 
   k=k0*pow((1+x0*x0/n),-(n+1.0)/2)/pt(x0,n);
  }
 if (alpha0<0.5) return x0; else return -x0;
}

/* 
Compute the p-value of 
beta distribution 
*/

double pbeta(double t,double r,double s)
{int i,k;
 double x,a,b,coeff,  S,a1,A0,A1,A2,B0,B1,B2,V;
 if (t==0) return 0.0; 
 if (t==1) return 1.0;
 if (t<=(r-1)/(r+s-2)) {x=t;a=r;b=s;} else {x=1-t;b=r;a=s;}
 coeff=exp(gamma(a+b)-gamma(a+1)-gamma(b)+a*log(x)+b*log(1-x));
 S=0.0;
 A0=B1=1.0;A1=B0=0.0;i=0;
 A2=A1+A0;B2=B1+B0;S=A2/B2;V=0.0;
 A0=A1;B0=B1;A1=A2;B1=B2;

 while(fabs((V-S)/S)>1e-12)
   {i++;k=i/2;V=S;
   if ((i%2)==0) a1=k*(b-k)*x/((a+2*k-1)*(a+2*k));
     else a1=-(a+k)*(a+b+k)*x/((a+2*k)*(a+2*k+1));
   A2=1+a1*A0/A1;B2=B1/A1+a1*B0/A1;S=A2/B2;
   A0=1.0;B0=B1/A1;B1=B2;A1=A2;
   }
 if (t<=(r-1)/(r+s-2)) S=coeff*S; else S=1-coeff*S;
 return S;
}

/* 
Compute the q-value of 
beta distribution 
*/

double qbeta(double alpha,double r,double s)
{double f,y,x1,x2,x;
 x1=0.0;x2=1.0;x=(x1+x2)/2;y=pbeta(x,r,s);
 while(((x2-x1)/x)>1e-12)
 {if ((y-alpha)>0) x2=x; else x1=x;
 x=(x1+x2)/2;y=pbeta(x,r,s);
  }
 return x;
}

/*
Compute the P-value of F-distribution
where m and n could be positive real number
*/

double pf(double x,double m,double n)
{double S,t,r,s;
 r=m/2;s=n/2;
 t=m*x/n;
 t=t/(1+t);
 S=pbeta(t,r,s);
 return S;
}

/*
Compute the Q-value of F-distribution
where m and n could be positive real number
*/
double qf(double alpha,double m,double n)
{double x,t,r,s;
 r=m/2;s=n/2;
 t=qbeta(alpha,r,s);
 x=n*t/(m*(1-t));
 return x;
}

/* Bessel J function */
double besselJ(double x,double alpha)
{int m,k;
  double S,aa,ff,absalpha,alphasq,xsq,PI=3.1415926535897932;
  if (x>20) {
    aa=sqrt(2/(PI*x))*cos(x-PI*(2*alpha+1)/4);
    xsq=x*x;alphasq=alpha*alpha;
    ff=1-(9-40*alphasq+16*alphasq*alphasq)/(128*xsq)+(11025-51664*alphasq+31584*alphasq*alphasq-5376*alphasq*alphasq*alphasq+256*alphasq*alphasq*alphasq*alphasq)/(98304*xsq*xsq);
    S=ff*aa;
    aa=sqrt(2/(PI*x))*((1-4*alphasq)/(8*x))*sin(x-PI*(2*alpha+1)/4);
    ff=1-(225-136*alphasq+16*alphasq*alphasq)/(384*xsq)+(893025-656784*alphasq+137824*alphasq*alphasq-10496*alphasq*alphasq*alphasq+256*alphasq*alphasq*alphasq*alphasq)/(491520*xsq*xsq);
    S+=ff*aa;
    return S;
  }

  if (alpha==floor(alpha))
    {absalpha=fabs(alpha);
      aa=exp(-gamma(absalpha+1)+absalpha*log(x/2));
      S=aa;m=1;
      aa=aa*(-1)*(x/2)*(x/2)/(m*(m+absalpha));
      S+=aa;
      while((fabs(aa/S)>1e-12)||(m<=(x/2)))
	{m++;aa=aa*(-1)*(x/2)*(x/2)/(m*(m+absalpha));
         S+=aa;
	}
      if ((alpha<0)&&(alpha/2!=floor(alpha/2))) S=-S;  
    return S;
    }

  if (alpha>-1){
   aa=exp(-gamma(alpha+1)+alpha*log(x/2));
  }
   else{
     ff=1;
     for(k=0;k<floor(-alpha);k++){
       ff=ff*(k+alpha+1);
     }
     aa=exp(-gamma(1+alpha+floor(-alpha))+alpha*log(x/2))*ff;
   }
  S=aa;
  m=1;
  aa=aa*(-1)*(x/2)*(x/2)/(m*(m+alpha));
  S+=aa;

  while((fabs(aa/S)>1e-12)||(m<=(x/2))){
    m++;
    aa=aa*(-1)*(x/2)*(x/2)/(m*(m+alpha));
    S+=aa;
  }
  return S;
}

/*Bessel Y function */
double besselY(double x,double alpha)
{int k;
  double absalpha,S,term1,term2,aa1,aa2,aa3,aa,alphasq,ff,xsq,PI=3.1415926535897932,Euler=0.57721566490153286;
  if (x>20){
    alphasq=alpha*alpha;
    xsq=x*x;
    aa=sqrt(2/(PI*x))*sin(x-PI*alpha/2-PI/4);
    ff=1-(16*alphasq*alphasq-40*alphasq+9)/(128*xsq)+(16*alphasq*(2*alphasq*(8*(alphasq-21)+987)-3229)+11025)/(98304*xsq*xsq);
    S=aa*ff;
    aa=sqrt(2/(PI*x))*((4*alphasq-1)/(8*x))*cos(x-PI*alpha/2-PI/4);
    ff=1+(16*alphasq*alphasq-40*alphasq+9)/(128*xsq)+(16*alphasq*(2*alphasq*(8*(alphasq-21)+987)-3229)+11025)/(98304*xsq*xsq);
    S+=aa*ff;
    return S;
  }
  if (alpha!=floor(alpha)){
    S=(besselJ(x,alpha)*cos(alpha*PI)-besselJ(x,-alpha))/sin(alpha*PI);
    return S;
  }

  absalpha=fabs(alpha);
  S=(2/PI)*log(x/2)*besselJ(x,alpha);
  aa=-exp(-absalpha*log(x/2)+gamma(absalpha))/PI;
  term1=aa;k=0;
  for(k=1;k<absalpha;k++){
    aa=aa*(x/2)*(x/2)/(k*(absalpha-k));
    term1+=aa;  
  }
  S+=term1;
  aa1=aa2=-Euler;
  for(k=1;k<=absalpha;k++){
    aa2+=1.0/k;
   }
  aa3=-exp(absalpha*log(x/2)-gamma(absalpha+1))/PI;
  aa=(aa1+aa2)*aa3;
  term2=aa;
  k=0;
  while((fabs(aa/term2)>1e-12)||(k<=(x/2))){
    k++;
    aa1+=1.0/k;
    aa2+=1.0/(k+absalpha);
    aa3*=(-1)*(x/2)*(x/2)/(k*(k+absalpha));
    aa=(aa1+aa2)*aa3;
    term2+=aa;
  }
  S+=term2;
  return S;
}

/*Bessel I function */
double besselI(double x,double alpha)
{int m,k;
  double S,aa,ff,absalpha,alphasq,xsq,PI=3.1415926535897932;
  if (x>20){
    alphasq=alpha*alpha;
    xsq=x*x;
    aa=exp(x)/sqrt(2*PI*x);
    ff=1+(1-4*alphasq)/(8*x)+(9-40*alphasq+16*alphasq*alphasq)/(128*xsq);
    S=aa*ff;
  }

  if (alpha==floor(alpha))
    {absalpha=fabs(alpha);
      aa=exp(-gamma(absalpha+1)+absalpha*log(x/2));
      S=aa;m=1;
      aa=aa*(x/2)*(x/2)/(m*(m+absalpha));
      S+=aa;
      while((fabs(aa/S)>1e-12)||(m<=(x/2)))
	{m++;aa=aa*(x/2)*(x/2)/(m*(m+absalpha));
         S+=aa;
	}
    return S;
    }

  if (alpha>-1){
   aa=exp(-gamma(alpha+1)+alpha*log(x/2));
  }
   else{
     ff=1;
     for(k=0;k<floor(-alpha);k++){
       ff=ff*(k+alpha+1);
     }
     aa=exp(-gamma(1+alpha+floor(-alpha))+alpha*log(x/2))*ff;
   }
  S=aa;
  m=1;
  aa=aa*(x/2)*(x/2)/(m*(m+alpha));
  S+=aa;

  while((fabs(aa/S)>1e-12)||(m<=(x/2))){
    m++;
    aa=aa*(x/2)*(x/2)/(m*(m+alpha));
    S+=aa;
  }
  return S;
}

/* Bessel K function, alpha can not be 0*/
double besselK(double x,double alpha)
{int k;
  double absalpha,S,term1,term2,aa1,aa2,aa3,aa,alphasq,ff,PI=3.1415926535897932,Euler=0.57721566490153286;
  if (x>10){
    alphasq=alpha*alpha;
    aa=sqrt(PI/(2*x))*exp(-x);
    ff=1+(4*alphasq-1)/(8*x)+(16*alphasq*alphasq-40*alphasq+9)/(128*x*x);
    S=aa*ff;
    return S;  
  }  
  if (alpha!=floor(alpha)){
    S=(PI/2)*(besselI(x,-alpha)-besselI(x,alpha))/sin(alpha*PI);
    return S;
  }

   absalpha=fabs(alpha);
   if((absalpha/2)==floor(absalpha/2))
    S=-besselI(x,absalpha)*log(x/2);
   else 
    S=besselI(x,absalpha)*log(x/2); 
   aa=exp(-absalpha*log(x/2)+gamma(absalpha))/2;
   term1=aa;k=0;
 
   for(k=1;k<absalpha;k++){
     aa=aa*(-1)*(x/2)*(x/2)/(k*(absalpha-k));
     term1+=aa;
   }
   S+=term1;
 
   aa1=aa2=-Euler;
   for(k=1;k<=absalpha;k++){
     aa2+=1.0/k;
   }

   if ((absalpha/2)==floor(absalpha/2))
     aa3=exp(absalpha*log(x/2)-gamma(absalpha+1))/2;
   else 
     aa3=(-1)*exp(absalpha*log(x/2)-gamma(absalpha+1))/2;
   aa=(aa1+aa2)*aa3;
   term2=aa;
   k=0;
   while((fabs(aa/term2)>1e-12)||(k<=(x/2))){
     k++;
     aa1+=1.0/k;
     aa2+=1.0/(k+absalpha);
     aa3*=(x/2)*(x/2)/(k*(k+absalpha));
     aa=(aa1+aa2)*aa3;
     term2+=aa;
   }
   S+=term2;
   return S;
}

/* Matern Correlation function, x is a real number, theta=(theta[0],theta[1])
 is the parameters in correlation function */
double matern(double x,double *theta,double alpha)
{double S,f1,f2;
 S=theta[0];
 if (x==0) return S;
 f1=theta[0]*exp(alpha*log(x/theta[1])-(alpha-1)*log(2.0)-gamma(alpha));
 if (alpha==0) 
   f2=besselK(x/theta[1],2)-besselK(x/theta[1],1)*2/(x/theta[1]);
 else f2=besselK(x/theta[1],alpha);
 S=f1*f2;
 return S;
}

/* Partial derivative of Matern correlation function 
 output vector of S, where S[0]=derivative of theta[0],
 and S[1]=derivative of theta[1]
*/
void dmatern(double x,double *theta,double alpha,double *S)
{double f1,f2;
 if (x==0) {S[0]=1;S[1]=0;return;}
 f1=exp(alpha*log(x/theta[1])-(alpha-1)*log(2.0)-gamma(alpha));
 if (alpha==0)
   f2=besselK(x/theta[1],2)-besselK(x/theta[1],1)*2/(x/theta[1]);
 else f2=besselK(x/theta[1],alpha);
 S[0]=f1*f2;
 f1=theta[0]*exp((alpha+1)*log(x)-(alpha+2)*log(theta[1])-(alpha-1)*log(2.0)-gamma(alpha));
 if (alpha==1)
   f2=besselK(x/theta[1],2)-besselK(x/theta[1],1)*2/(x/theta[1]);
 else f2=besselK(x/theta[1],alpha-1);
 S[1]=f1*f2;
 return;
}

/* Second order partial derivative of Matern correlation function 
 output vector is S with dimension 3, where S[0]=derivative of 
 (theta[0],theta[0]), S[1]=derivative of (theta[1],theta[1]), 
 S[2]=derivative of (theta[0],theta[1])
*/
void d2matern(double x,double *theta,double alpha,double *S)
{double f01,f02,f1,f2;
 if (x==0) {S[0]=S[1]=S[2]=0;return;}
 S[0]=0;
 f01=theta[0]*exp((alpha+2)*log(x)-(alpha+4)*log(theta[1])-(alpha-1)*log(2.0)-gamma(alpha));
 f02=theta[0]*(2*alpha+1)*exp((alpha+1)*log(x)-(alpha+3)*log(theta[1])-(alpha-1)*log(2.0)-gamma(alpha));
 if (alpha==0) f1=besselK(x/theta[1],2)-besselK(x/theta[1],1)*2/(x/theta[1]);
 else f1=besselK(x/theta[1],alpha);
 if (alpha==1) f2=besselK(x/theta[1],2)-besselK(x/theta[1],1)*2/(x/theta[1]);
 else f2=besselK(x/theta[1],alpha-1);
 S[1]=f01*f1-f02*f2;
 f1=exp((alpha+1)*log(x)-(alpha+2)*log(theta[1])-(alpha-1)*log(2.0)-gamma(alpha));
 if (alpha==1) f2=besselK(x/theta[1],2)-besselK(x/theta[1],1)*2/(x/theta[1]);
 else f2=besselK(x/theta[1],alpha-1);
 S[2]=f1*f2;
 return;
}

/*
sort, order, and rank of vector in increasing order, input vector A, output vector B, length of vector nn
*/
void mysort(double *A,double *B,int nn)
{int i,j;
 double f;
 for(i=0;i<nn;i++)
   B[i]=A[i];
 for(i=0;i<nn-1;i++)
   {for(j=i+1;j<nn;j++)
     if (B[i]>B[j]) {f=B[i];B[i]=B[j];B[j]=f;}
   }
 return;
}

void myorder(double *A,int *B,int nn)
{int i,j,k;
 double f,*C;
 C=new double[nn];
 for(i=0;i<nn;i++)
   {B[i]=i+1;C[i]=A[i];}
 for(i=0;i<nn-1;i++)
   {for(j=i+1;j<nn;j++)
     if (C[i]>C[j]) 
       {f=C[i];C[i]=C[j];C[j]=f;
       k=B[i];B[i]=B[j];B[j]=k;}
   }
 delete []C;
 return;
}

void myrank(double *A,double *B,int nn)
{int i,j,less,equal;
 for(i=0;i<nn;i++)
   {less=0;equal=0;
   for(j=0;j<nn;j++) 
     {less+=(A[j]<A[i]);equal+=(A[j]==A[i]);}
   B[i]=less+(equal+1)/2.0;
   }
 return;
}

/*
generate the nn*nn weight matrix of Moran's I
*/
void neighbour(double *location,int nn,double *weight)
{int i,j,kmin,*distorder;
 double *hh,*dist,p1[2],p2[2],S;
 hh=new double[nn*nn];
 dist=new double[nn];
 distorder=new int[nn];
 for(i=0;i<nn;i++)
   {p1[0]=location[2*i];p1[1]=location[2*i+1];
     for(j=0;j<nn;j++)
       {p2[0]=location[2*j];p2[1]=location[2*j+1];
	 hh[i*nn+j]=spheric(p1,p2);
       }
   }
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     weight[i*nn+j]=0;

 kmin=(int)floor(sqrt(nn*1.0));
 if (kmin>10) kmin=10;
 
 for(i=0;i<nn;i++)
   {for(j=0;j<nn;j++)
     dist[j]=hh[i*nn+j];
   myorder(dist,distorder,nn);
   for(j=0;j<=kmin;j++)
     {weight[i*nn+distorder[j]-1]=1.0/kmin;
      weight[i*nn+i]=0;
     }
   }

 delete []hh;
 delete []dist;
 delete []distorder;
 return;
}

/*
Moran's I statistic: output 5 dimensional vector: I-value, expected-value, variance, z-value, p-value
*/
void moranI(double *yobs,double *ww,int nn,double *output)
{int i,j;
 double S0,S1,S2,b2,b4,wi,ybar,term1,term2,term3;
 S0=S1=S2=b2=b4=ybar=0;
 for(i=0;i<nn;i++) 
   {ww[i*nn+i]=0;ybar+=yobs[i]/nn;}
 for(i=0;i<nn;i++)
   {wi=0;
    b2+=(yobs[i]-ybar)*(yobs[i]-ybar)/nn;
    b4+=(yobs[i]-ybar)*(yobs[i]-ybar)*(yobs[i]-ybar)*(yobs[i]-ybar)/nn;
    for(j=0;j<nn;j++)
      {wi+=ww[i*nn+j]+ww[j*nn+i];
       S1+=(ww[i*nn+j]+ww[j*nn+i])*(ww[i*nn+j]+ww[j*nn+i])/2; 
      }
	S0+=wi/2;
    S2+=wi*wi;
   }

 for(i=0;i<4;i++)
   output[i]=0;

 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     output[0]+=ww[i*nn+j]*(yobs[i]-ybar)*(yobs[j]-ybar)/(S0*b2);

 output[1]=-1.0/(nn-1);
 term1=S1*(nn*b2*b2-b4)/(S0*S0*b2*b2*(nn-1));
 term2=(S2-2*S1)*(2*b4-nn*b2*b2)/(S0*S0*b2*b2*(nn-1)*(nn-2));
 term3=(S0*S0-S2+S1)*(3*nn*b2*b2-6*b4)/(S0*S0*b2*b2*(nn-1)*(nn-2)*(nn-3));
 output[2]=term1+term2+term3-output[1]*output[1];
 output[3]=(output[0]-output[1])/sqrt(output[2]);
 output[4]=2*(1-pnorm(fabs(output[3])));
 return;
}

/* transpose of a matrix(first input matrix not change,second output matrix change)
 dimension of input matrix is nn*mm
 dimension of output matrix is mm*nn
*/
void transpose(double *A,double *B,int nn,int mm)
{int i,j;
 double aa;
 for(i=0;i<nn;i++)
   for(j=0;j<mm;j++)
     B[j*nn+i]=A[i*mm+j];
 return;
}

/* product of matrix(first, second input not change,third out put matrix change) dimension of input matrix A is nn*kk
 dimension of input matrix B is kk*mm
 dimension of output matrix C is nn*mm
*/
void prod(double *A,double *B,double *C,int nn,int kk,int mm)
{int i,j,k;
 for(i=0;i<nn;i++)
  for(j=0;j<mm;j++)
  {C[i*mm+j]=0;
     for(k=0;k<kk;k++)
      C[i*mm+j]+=A[i*kk+k]*B[k*mm+j];
   }
 return;
}

/* sweep of a square matrix(input matrix will change) 
 dimension of input(output) matrix A is nn*nn;
 sweep only to be done for (row,column)=(ii,ii);
*/
void sweep(double *A,int ii,int nn)
{int i,j;
 A[ii*nn+ii]=-1/A[ii*nn+ii];
 for(i=0;i<nn;i++)
   if (i!=ii) {A[ii*nn+i]=-A[ii*nn+i]*A[ii*nn+ii];
               A[i*nn+ii]=-A[i*nn+ii]*A[ii*nn+ii];}
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     if ((i!=ii)&&(j!=ii)) A[i*nn+j]=A[i*nn+j]+A[i*nn+ii]*A[ii*nn+j]/A[ii*nn+ii];
 return;
}

/* inverse of matrix(first input matrix not change, second output 
 matrix is the inverse)
  dimension of input matrix A is nn*nn
  dimension of output matrix B is nn*nn
  */ 
void inverse(double *A,double *B,int nn)
{int i,j,k;
 i=0;
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     B[i*nn+j]=-A[i*nn+j];
 for(i=0;i<nn;i++)
   sweep(B,i,nn);
 return;
}

/* least square method(input matrix willl change)
   dimension of matrix of A is nn*nn */
void ls(double *A,int nn)
{int i,j;
  for(i=0;i<nn-1;i++)
    sweep(A,i+1,nn);
  for(i=0;i<nn;i++)
    for(j=0;j<nn;j++)
      if ((i!=0)&&(j!=0)) A[i*nn+j]=-A[i*nn+j];
  return;
}

/*regression analysis (input response and design matrix of independent 
  variables) last two columns are output. Location input ll will be used in calculating moran's I statistic.
Output: The first one is 2*pp+1 dimensional vector. It counts  pp estimates, pp corresponding std, and squared root of the MSE (sigma). The second one is the model residual. It is a nn-dimensional vector of model residuals.
*/
double regression(double *yobs,double *xobs,double *ll,int nn,int pp,double *output,double *residual)
{int i,j,k;
 ofstream outfile;
 double *xy,*txy,*out,*yhat,*rsort,*SS1,*SS2,*SSout,*ww,*estimate;
 double tvalue,pvalue,ybar,SSR,SSE,SST,dfR,dfE,dfT,fvalue,loglike;

 xy=new double[nn*(pp+1)];
 txy=new double[nn*(pp+1)];
 for(i=0;i<nn;i++)
   {xy[i*(pp+1)]=yobs[i];
    for(j=0;j<pp;j++)
      xy[i*(pp+1)+j+1]=xobs[i*pp+j];
   }

 transpose(xy,txy,nn,pp+1);

 out=new double[(pp+1)*(pp+1)];
 prod(txy,xy,out,pp+1,nn,pp+1);

 SS1=new double[pp];
 SS2=new double[pp];
 SSout=new double[pp*pp];
 
 for(k=2;k<pp+1;k++)
   {for(i=0;i<pp;i++)
     for(j=0;j<pp;j++)
       {if((i<k)&&(j<k)) SSout[i*pp+j]=out[i*(pp+1)+j];
       else if ((i>=k)&&(j<k)) SSout[i*pp+j]=out[(i+1)*(pp+1)+j];
       else if ((i<k)&&(j>=k)) SSout[i*pp+j]=out[i*(pp+1)+j+1];
       else if ((i>=k)&&(j>=k)) SSout[i*pp+j]=out[(i+1)*(pp+1)+j+1];
       }
   ls(SSout,pp);
   SS2[k-1]=SSout[0];
   }
 
 for(i=0;i<pp;i++)
   {sweep(out,i+1,pp+1);
   SS1[i]=out[0];
   }
 for(i=0;i<pp-1;i++)
   SS1[i]=SS1[i]-SS1[i+1];
    
 for(i=0;i<pp+1;i++)
   for(j=0;j<pp+1;j++)
    if ((i!=0)&&(j!=0)) out[i*(pp+1)+j]=-out[i*(pp+1)+j];
 
 output[2*pp]=sqrt(out[0]/(nn-pp));
 for(i=0;i<pp;i++)
   {output[i]=out[i+1];    
   output[pp+i]=sqrt(output[2*pp]*out[(i+1)*(pp+1)+(i+1)]);
   }

 yhat=new double[nn];
 rsort=new double[nn];
 for(i=0;i<nn;i++)
   {yhat[i]=0;residual[i]=0;
   for(j=0;j<pp;j++)
     yhat[i]+=xobs[i*pp+j]*output[j];
   residual[i]=yobs[i]-yhat[i];
   }
 mysort(residual,rsort,nn);

 ybar=0;
 SST=SSE=SSR=0;
 dfT=nn-1;dfE=nn-pp;dfR=pp-1;
 for(i=0;i<nn;i++)
   ybar+=yobs[i]/nn;

 for(i=0;i<nn;i++)
   {SST+=(yobs[i]-ybar)*(yobs[i]-ybar);
    SSE+=(yobs[i]-yhat[i])*(yobs[i]-yhat[i]);
    SSR+=(yhat[i]-ybar)*(yhat[i]-ybar);
   } 

 outfile.open("coeff.reg",ios::out);
 outfile <<"Estimate"<<'\t'<<"Std"<<'\t'<<"t-value"<<'\t'<<"p-value"<<endl;
 for(i=0;i<pp;i++)
   {tvalue=output[i]/output[pp+i];
    pvalue=2*(1-pt(fabs(tvalue),nn-pp));
    outfile <<output[i]<<'\t'<<output[pp+i]<<'\t'<<tvalue<<'\t'<<pvalue<<endl;
   } 
 outfile.close();

 outfile.open("mse.reg",ios::out);
 outfile <<output[2*pp]*output[2*pp]<<endl;
 outfile.close();

 loglike=-nn*log(2*3.1415926)/2-nn*log(output[2*pp])-(nn-pp)/2;
 outfile.open("aicbic.reg",ios::out);
 outfile <<"AIC"<<'\t'<<"BIC"<<endl;
 outfile <<-2*loglike+2*pp<<'\t'<<-2*loglike+pp*log(nn*1.0);
 outfile.close();

 outfile.open("R2.reg",ios::out);
 outfile <<SSR/SST<<endl;
 outfile.close();
  
 outfile.open("residualplot.reg",ios::out);
 outfile <<"predict"<<'\t'<<"residual"<<endl;
 for(i=0;i<nn;i++)
   outfile <<yhat[i]<<'\t'<<residual[i]<<endl;
 outfile.close();
 
 outfile.open("qqplot.reg",ios::out);
 outfile <<"normal-quantile"<<'\t'<<"residual"<<endl;
 for(i=0;i<nn;i++)
   outfile <<qnorm((i+1.0)/(nn+1.0))<<'\t'<<rsort[i]<<endl;
 outfile.close();

 outfile.open("anova.reg",ios::out);
 outfile <<"\tdf\tSS\tMS\tFvalue\tPvalue\n";
 for(i=0;i<pp-1;i++)
   {fvalue=SS1[i]/(SSE/dfE);
    pvalue=1-pf(fvalue,1,dfE);
    outfile <<"NO."<<i+1<<"\t"<<1<<'\t'<<SS1[i]<<'\t'<<SS1[i]<<'\t'<<fvalue<<'\t'<<pvalue<<endl;
   }
 outfile <<"Error\t"<<dfE<<'\t'<<SSE<<'\t'<<SSE/dfE<<endl;
 outfile <<"Total\t"<<dfT<<'\t'<<SST<<endl;
 outfile.close();

 outfile.open("significance.reg",ios::out);
 outfile <<"\tdf\tSS\tMS\tFvalue\tPvalue\n";
 for(i=0;i<pp-1;i++)
   {fvalue=(SS2[i+1]-SS2[0])/(SSE/dfE);
    pvalue=1-pf(fvalue,1,dfE);
    outfile <<"NO."<<i+1<<"\t"<<1<<'\t'<<SS2[i+1]-SS2[0]<<'\t'<<SS2[i+1]-SS2[0]<<'\t'<<fvalue<<'\t'<<pvalue<<endl;
   }
 outfile <<"Error\t"<<dfE<<'\t'<<SSE<<'\t'<<SSE/dfE<<endl;
 outfile <<"Total\t"<<dfT<<'\t'<<SST<<endl;
 outfile.close();

 ww=new double[nn*nn];
 neighbour(ll,nn,ww);
 estimate=new double[5];
 moranI(residual,ww,nn,estimate);
 outfile.open("moranI.reg",ios::out);
 outfile <<"Ivalue\t"<<"Mean\t"<<"Std\t"<<"Z_value\t"<<"P-value\n";
 for(i=0;i<5;i++)
   outfile <<estimate[i]<<'\t'; outfile <<endl;
 outfile.close();
 pvalue=estimate[4];

 outfile.open("estimate.reg",ios::out);
 outfile <<nn<<'\t'<<pp<<endl;
 for(i=0;i<pp+1;i++)
   {for(j=0;j<pp+1;j++)
     outfile <<out[i*(pp+1)+j]<<'\t';
   outfile <<endl;
   }
 outfile.close();

 delete []ww;
 delete []estimate;
 delete []xy;
 delete []txy;
 delete []out;
 delete []yhat;
 delete []rsort;
 delete []SS1;
 delete []SS2;
 delete []SSout;
 return pvalue;
}



/* trace of matrix(input matrix not change)
   dimension of matrix of A is nn*nn */
double trace(double *A,int nn)
{int i,j;
 double S(0);
 for(i=0;i<nn;i++)
   S+=A[i*nn+i];
 return S;
}

/* log of the determinant of matrix with sign(input matrix not change)
  dimension of matrix of A is nn*nn   */
double detlog(double *A,int nn)
{int i,j;
  double S(0),*AA;
  AA=new double[nn*nn];
  for(i=0;i<nn;i++)
    for(j=0;j<nn;j++)
      AA[i*nn+j]=A[i*nn+j];
 for(i=0;i<nn;i++)
  {S+=log(AA[i*nn+i]);
  sweep(AA,i,nn);
  }
 delete []AA;
 return S;
}

/* profile loglikelihood function matern correlation function, and first order, second order partial loglikelihood in krig method
 matrix used:  
 Profile likelihood:  
              Mtemp1=R.inverse%*%X, 
              Mtemp2=t(X), 
              Mtemp3=t(X)%*%R.inverse%*%X, 
              Mtemp4=solve(Mtemp3),
              Mtemp5=Mtemp1%*%Mtemr4, 
              Mterm6=Mterm5*%*t(Mtemp1), 
              Mterm7=t(Mtemp1),
  Partial derivative of Profile likelihood: */
int loglikelihood(double *yy,double *xx,double *hh,double *theta,double alpha,int nn,int pp,int index,double *f)
{int i,j,k;
 double S,YMY(0),dYMY[2],d2YMY[3],V1[2],V2[3],PI=3.1415926535897932;
 double *RR,*RRinverse,*MM,*Mtemp1,*Mtemp2,*Mtemp3,*Mtemp4,*Mtemp5,*Mtemp6,*Mtemp7,*Mtemp8,*Mtemp9,*Mtemp10,*Mtemp11,*Mtemp12;
 f[0]=f[1]=f[2]=f[3]=f[4]=f[5]=0.0;
 S=-(nn/2.0)*(1+log(2*PI/nn));
 RR=new double[nn*nn];
  
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     if (i==j) RR[i*nn+j]=1.0; else RR[i*nn+j]=matern(hh[i*nn+j],theta,alpha);
 S+=-detlog(RR,nn)/2;

 RRinverse=new double[nn*nn];
 inverse(RR,RRinverse,nn);

 Mtemp2=new double[nn*pp];
 transpose(xx,Mtemp2,nn,pp);
 
 Mtemp1=new double[nn*pp];
 prod(RRinverse,xx,Mtemp1,nn,nn,pp);

 Mtemp7=new double[nn*pp];
 transpose(Mtemp1,Mtemp7,nn,pp);

 Mtemp3=new double[pp*pp];
 prod(Mtemp2,Mtemp1,Mtemp3,pp,nn,pp);
 delete []Mtemp2;
 
 Mtemp4=new double[pp*pp];
 inverse(Mtemp3,Mtemp4,pp);
 delete []Mtemp3;

 Mtemp5=new double[nn*pp];
 prod(Mtemp1,Mtemp4,Mtemp5,nn,pp,pp);
 delete []Mtemp1;
 delete []Mtemp4;

 Mtemp6=new double[nn*nn];
 prod(Mtemp5,Mtemp7,Mtemp6,nn,pp,nn);

 MM=new double[nn*nn];
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     MM[i*nn+j]=RRinverse[i*nn+j]-Mtemp6[i*nn+j];

 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     YMY+=yy[i]*MM[i*nn+j]*yy[j];
 S+=-log(YMY)*nn/2;
 
 f[0]=S;

 delete []Mtemp5;
 delete []Mtemp6;
 delete []Mtemp7;
 if (index==0) 
  {k=0;
  delete []RR; delete []RRinverse; delete []MM; return k;
 } 
 k=1;

  Mtemp1=new double[nn*nn]; //dRR1
  Mtemp2=new double[nn*nn]; //dRR2 
  Mtemp6=new double[nn*nn];  //RRinverse%*%dRR1
  Mtemp7=new double[nn*nn];  //RRinverse%*%dRR2 
  for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     {dmatern(hh[i*nn+j],theta,alpha,V1);
      if (i==j) Mtemp1[i*nn+j]=0; else Mtemp1[i*nn+j]=V1[0];
      Mtemp2[i*nn+j]=V1[1];
     }
  prod(RRinverse,Mtemp1,Mtemp6,nn,nn,nn); 
  prod(RRinverse,Mtemp2,Mtemp7,nn,nn,nn);

  Mtemp8=new double[nn*nn];  //MM%*%dRR1%*%MM
  Mtemp9=new double[nn*nn];  //MM%*%dRR2%*%MM
  Mtemp12=new double[nn*nn]; //Intermediate Matrix
  prod(MM,Mtemp1,Mtemp12,nn,nn,nn);
  prod(Mtemp12,MM,Mtemp8,nn,nn,nn);
  prod(MM,Mtemp2,Mtemp12,nn,nn,nn);
  prod(Mtemp12,MM,Mtemp9,nn,nn,nn);

  dYMY[0]=dYMY[1]=d2YMY[0]=d2YMY[1]=d2YMY[2]=0;
  for(i=0;i<nn;i++)
    for(j=0;j<nn;j++)
      {dYMY[0]+=yy[i]*Mtemp8[i*nn+j]*yy[j];
       dYMY[1]+=yy[i]*Mtemp9[i*nn+j]*yy[j];
      }

  f[1]=-trace(Mtemp6,nn)/2+dYMY[0]*nn/(2*YMY);
  f[2]=-trace(Mtemp7,nn)/2+dYMY[1]*nn/(2*YMY);
  f[3]=dYMY[0]*dYMY[0]*nn/(2*YMY*YMY);
  f[4]=dYMY[1]*dYMY[1]*nn/(2*YMY*YMY);
  f[5]=dYMY[1]*dYMY[0]*nn/(2*YMY*YMY);

  Mtemp10=new double[nn*nn]; //MM%*%dRR1%*%MM%*%dRR1%*%MM; or MM%*%dRR2%*%MM%*%dRR1%*%MM; or MM%*%dRR1%*%MM%*%dRR2%*%MM;
  Mtemp11=new double[nn*nn]; //RRinverse%*%dRR1%*%RRinverse%*%dRR1; or RRinverse%*%dRR2%*%RRinverse%*%dRR2; or  RRinverse%*%dRR1%*%RRinverse%*%dRR2;
  prod(Mtemp6,Mtemp6,Mtemp11,nn,nn,nn);
  f[3]+=trace(Mtemp11,nn)/2;
  prod(Mtemp7,Mtemp7,Mtemp11,nn,nn,nn);
  f[4]+=trace(Mtemp11,nn)/2; 
  prod(Mtemp6,Mtemp7,Mtemp11,nn,nn,nn);
  f[5]+=trace(Mtemp11,nn)/2;

  prod(Mtemp8,Mtemp1,Mtemp12,nn,nn,nn);  //(dRR1,dRR1)
  prod(Mtemp12,MM,Mtemp10,nn,nn,nn);
  for(i=0;i<nn;i++)
    for(j=0;j<nn;j++)
     d2YMY[0]+=yy[i]*Mtemp10[i*nn+j]*yy[j];
  f[3]+=-d2YMY[0]*nn/YMY;

  prod(Mtemp9,Mtemp2,Mtemp12,nn,nn,nn);  //(dRR2,dRR2)
  prod(Mtemp12,MM,Mtemp10,nn,nn,nn);
  for(i=0;i<nn;i++)
    for(j=0;j<nn;j++)
     d2YMY[1]+=yy[i]*Mtemp10[i*nn+j]*yy[j];
  f[4]+=-d2YMY[1]*nn/YMY;

  prod(Mtemp8,Mtemp2,Mtemp12,nn,nn,nn);  //(dRR1,dRR2)
  prod(Mtemp12,MM,Mtemp10,nn,nn,nn);
  for(i=0;i<nn;i++)
    for(j=0;j<nn;j++)
     d2YMY[2]+=yy[i]*Mtemp10[i*nn+j]*yy[j];
  f[5]+=-d2YMY[2]*nn/YMY;

  delete []Mtemp1;
  delete []Mtemp2;
  delete []Mtemp6;
  delete []Mtemp7;
  delete []Mtemp8;
  delete []Mtemp9;
  delete []Mtemp10;

 Mtemp3=new double[nn*nn]; //d2RR11
 Mtemp4=new double[nn*nn]; //d2RR22
 Mtemp5=new double[nn*nn]; //d2RR12
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     {d2matern(hh[i*nn+j],theta,alpha,V2);
      Mtemp3[i*nn+i]=V2[0];
      Mtemp4[i*nn+j]=V2[1];
      Mtemp5[i*nn+j]=V2[2];
     }

 prod(RRinverse,Mtemp3,Mtemp11,nn,nn,nn); //Mtemp11=RRinverse%*%d2R11 
 f[3]+=-trace(Mtemp11,nn)/2;
 prod(RRinverse,Mtemp4,Mtemp11,nn,nn,nn); //Mtemp11=RRinverse%*%d2R22 
 f[4]+=-trace(Mtemp11,nn)/2;
 prod(RRinverse,Mtemp5,Mtemp11,nn,nn,nn); //Mtemp11=RRinverse%*%d2R12 
 f[5]+=-trace(Mtemp11,nn)/2;

 d2YMY[0]=d2YMY[1]=d2YMY[2]=0;
 prod(MM,Mtemp3,Mtemp12,nn,nn,nn); 
 prod(Mtemp12,MM,Mtemp11,nn,nn,nn); //Mtemp11=MM%*%d2R11%*%MM
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     d2YMY[0]+=yy[i]*Mtemp11[i*nn+j]*yy[j];
 f[3]+=d2YMY[0]*nn/(2*YMY);
 prod(MM,Mtemp4,Mtemp12,nn,nn,nn); 
 prod(Mtemp12,MM,Mtemp11,nn,nn,nn); //Mtemp11=MM%*%d2R22%*%MM
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     d2YMY[1]+=yy[i]*Mtemp11[i*nn+j]*yy[j];
 f[4]+=d2YMY[1]*nn/(2*YMY);
 prod(MM,Mtemp5,Mtemp12,nn,nn,nn); 
 prod(Mtemp12,MM,Mtemp11,nn,nn,nn); //Mtemp11=MM%*%d2R12%*%MM
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     d2YMY[2]+=yy[i]*Mtemp11[i*nn+j]*yy[j];
 f[5]+=d2YMY[2]*nn/(2*YMY);

  delete []RR;
  delete []RRinverse;
  delete []MM;
  delete []Mtemp3;
  delete []Mtemp4;
  delete []Mtemp5;
  delete []Mtemp11;
  delete []Mtemp12;

 return k;
}
/* Profile Estimate of theta for matern correlation function
 Input: yy=nn-dimensional response, xx=nn*pp dimensional independent variables
        hh=nn*nn dimensional distance maternx, 
        theta=initital guess (will be the output of the estimate of theta)
        alpha=smooth parameter
        nn=row, pp=column
Output: theta (estimate)
Return value: 0=convergence, 1=not convergence.
*/
int thetaest(double *yy,double *xx,double *hh,double *theta,double alpha,int nn,int pp)
{int i,k,index;
 double ff[6],diff;
 k=0;index=1; i=0; 
 k=loglikelihood(yy,xx,hh,theta,alpha,nn,pp,index,ff);
 diff=fabs(ff[1])+fabs(ff[2]);

 theta[0]=theta[0]-(ff[4]*ff[1]-ff[5]*ff[2])/(ff[3]*ff[4]-ff[5]*ff[5]);
 theta[1]=theta[1]-(ff[3]*ff[2]-ff[5]*ff[1])/(ff[3]*ff[4]-ff[5]*ff[5]);

 while ((i<10)&&(diff>1e-7)){
   i++;
   k=loglikelihood(yy,xx,hh,theta,alpha,nn,pp,index,ff);
   diff=fabs(ff[1])+fabs(ff[2]);
   theta[0]=theta[0]-(ff[4]*ff[1]-ff[5]*ff[2])/(ff[3]*ff[4]-ff[5]*ff[5]);
   theta[1]=theta[1]-(ff[3]*ff[2]-ff[5]*ff[1])/(ff[3]*ff[4]-ff[5]*ff[5]);
 }

 k=(diff>1e-7);
 return k;
}

/* Generalized Least Square Method for Multiple Regression Model
 where dimension of input vector yy is nn
       dimension of input matrix xx is nn*pp
       dimension of input Sigma is nn*pp
       dimension of output matrix is (pp+1)*(pp+1)
 The output matrix is composed of 
  SS       hatbeta'
  hatbeta  cov(hatbeta)
  where hatbeta is the column vector of the estimate of coefficients
*/
void gls(double *yy,double *xx,double *Sigma,int nn,int pp,double *output)
{int i,j;
 double *xy,*Mtemp1,*Mtemp2,*Mtemp3; 

 xy=new double[nn*(pp+1)];
 for(i=0;i<nn;i++)
   for(j=0;j<=pp;j++)
     if (j<pp) xy[i*(pp+1)+j]=xx[i*pp+j];else xy[i*(pp+1)+j]=yy[i];

 Mtemp1=new double[nn*nn]; // Mtemp1=solve(Sigma)
 inverse(Sigma,Mtemp1,nn);

 Mtemp2=new double[nn*(pp+1)]; //Mtemp=t(xy)
 transpose(xy,Mtemp2,nn,pp+1);

 Mtemp3=new double[nn*(pp+1)]; //Mtemp3=t(xy)%*%solve(Sigma)
 prod(Mtemp2,Mtemp1,Mtemp3,(pp+1),nn,nn);

 prod(Mtemp3,xy,output,(pp+1),nn,(pp+1));
 for(i=0;i<pp;i++)
   sweep(output,i,(pp+1));

 for(i=0;i<pp;i++)
   for(j=0;j<pp;j++)
     output[i*(pp+1)+j]=-output[i*(pp+1)+j];

 delete []xy;
 delete []Mtemp1;
 delete []Mtemp2;
 delete []Mtemp3;
 return;
}
/* Estimate of beta and sigmasq given theta for matern correlation function 
  input vecotr yy is nn-dimensional response, 
  input matrix xx is nn*pp-matrix of independent variable
  input matrix location is nn*2 matrix of location index
  output the estimate of sigma, beta, and standard error of beta
  dimension of output 2*pp+1, first pp are estimates of beta
                              next pp are standard error of the estimates
                              last one, MSE.
 */
void betacondtheta(double *yy,double *xx,double *location,double *theta,double alpha,int nn,int pp,double *output)
{int i,j;
 double p1[2],p2[2];
 double *CORR,*est;

 CORR=new double[nn*nn];
 for(i=0;i<nn;i++)
   {p1[0]=location[2*i];p1[1]=location[2*i+1];
     for(j=0;j<nn;j++)
       {p2[0]=location[2*j];p2[1]=location[2*j+1];
	 CORR[i*nn+j]=spheric(p1,p2);
	 if (i==j) CORR[i*nn+j]=1;
         else CORR[i*nn+j]=matern(CORR[i*nn+j],theta,alpha);
       }
   }

 est=new double[(pp+1)*(pp+1)];
 gls(yy,xx,CORR,nn,pp,est);

 output[2*pp]=sqrt(est[(pp+1)*(pp+1)-1]/nn);
 for(i=0;i<pp;i++)
   {output[i]=est[i*(pp+1)+pp];
    output[pp+i]=sqrt(est[i*(pp+1)+i])*output[2*pp];
   }
 
  delete []CORR;
  delete []est;
  return;
}

/* MLE of theta, beta and sigmasq for matern correlation function 
 we assume the MLE of theta[0] in [0.01,0.99]
    and theta[1] in [10,5000]
 if any is at the boundary, then output 1 for not converence
 if both witn then out 0 for convergence 
 output the estimate of theta, beta, standard error of beta, and sigma
  dimension of output 2*pp+3.
 residual is a nn-dimensional vector of residual. It contains the residual of the model when theta is given. I may also want to use Moran's I to detect dependence.   
*/
int mle(double *yy,double *xx,double *location,double alpha,int nn,int pp,double *output,double *residual)
{int i,j,k,l,index;
 ofstream outfile;
 double *xy,*txy,*CORR,*CORRinverse,*out,*out1,*SS1,*SS2,*yhat,*rsort,*SSout,*Liketemp,*hh,*indexmax;
 double tvalue,pvalue,fvalue,diff,ff[6],theta[2],p1[2],p2[2],bound[4],step[2],SSE,dfE,SST,dfT;

 hh=new double[nn*nn];
 for(i=0;i<nn;i++)
   {p1[0]=location[2*i];p1[1]=location[2*i+1];
     for(j=0;j<nn;j++)
       {p2[0]=location[2*j];p2[1]=location[2*j+1];
	 hh[i*nn+j]=spheric(p1,p2);
       }
   }

 Liketemp=new double[75];
 indexmax=new double[3];
 k=0;index=0;

 bound[0]=0.0;bound[1]=0.99;bound[2]=1.0;bound[3]=5000.0;
 step[0]=(bound[1]-bound[0])/4;
 step[1]=(bound[3]-bound[2])/4;

 for(l=0;l<12;l++){
 for(i=0;i<5;i++)
   for(j=0;j<5;j++)
     {theta[0]=Liketemp[3*(5*i+j)]=bound[0]+i*step[0];
      theta[1]=Liketemp[3*(5*i+j)+1]=bound[2]+j*step[1];
      k=loglikelihood(yy,xx,hh,theta,alpha,nn,pp,index,ff);
      Liketemp[3*(5*i+j)+2]=ff[0];
	 }

 indexmax[0]=Liketemp[0];indexmax[1]=Liketemp[1];indexmax[2]=Liketemp[2];
 for(i=0;i<25;i++)
   if (indexmax[2]<=Liketemp[3*i+2]) 
     {indexmax[0]=Liketemp[3*i];
      indexmax[1]=Liketemp[3*i+1];
      indexmax[2]=Liketemp[3*i+2];
     }
    
 if ((indexmax[0]==bound[0])||(indexmax[0]==bound[1]))
   {bound[0]=indexmax[0]-2*step[0];bound[1]=indexmax[0]+2*step[0];}
 else 
   {bound[0]=indexmax[0]-step[0];bound[1]=indexmax[0]+step[0];}

 if ((indexmax[1]==bound[2])||(indexmax[1]==bound[3]))
   {bound[2]=indexmax[1]-2*step[1];bound[3]=indexmax[1]+2*step[1];}
 else 
   {bound[2]=indexmax[1]-step[1];bound[3]=indexmax[1]+step[1];}

 if (bound[0]<0) bound[0]=0;
 if (bound[1]>0.99) bound[1]=0.99;
 if (bound[2]<10) bound[2]=1;
 if (bound[3]>5000) bound[3]=5000;
	      
 step[0]=(bound[1]-bound[0])/4;
 step[1]=(bound[3]-bound[2])/4;
 }

 k=0;index=1; i=0; 
 k=loglikelihood(yy,xx,hh,theta,alpha,nn,pp,index,ff);
 diff=fabs(ff[1])+fabs(ff[2]);
 theta[0]=theta[0]-(ff[4]*ff[1]-ff[5]*ff[2])/(ff[3]*ff[4]-ff[5]*ff[5]);
 theta[1]=theta[1]-(ff[3]*ff[2]-ff[5]*ff[1])/(ff[3]*ff[4]-ff[5]*ff[5]);

 while ((i<10)&&(diff>1e-7)){
   i++;
   k=loglikelihood(yy,xx,hh,theta,alpha,nn,pp,index,ff);
   diff=fabs(ff[1])+fabs(ff[2]);
   theta[0]=theta[0]-(ff[4]*ff[1]-ff[5]*ff[2])/(ff[3]*ff[4]-ff[5]*ff[5]);
   theta[1]=theta[1]-(ff[3]*ff[2]-ff[5]*ff[1])/(ff[3]*ff[4]-ff[5]*ff[5]);
 }
 k=(diff>1e-7);

 delete []Liketemp;
 Liketemp=new double[2*pp+1];

 xy=new double[nn*(pp+1)];
 txy=new double[nn*(pp+1)];
 CORR=new double[nn*nn];
 CORRinverse=new double[nn*nn];
 out1=new double[nn*(pp+1)];
 out=new double[(pp+1)*(pp+1)];
 for(i=0;i<nn;i++)
   {xy[i*(pp+1)]=yy[i];
    for(j=0;j<pp;j++)
      xy[i*(pp+1)+j+1]=xx[i*pp+j];
   }
 transpose(xy,txy,nn,pp+1);
 for(i=0;i<nn;i++)
   {p1[0]=location[2*i];p1[1]=location[2*i+1];
     for(j=0;j<nn;j++)
       {p2[0]=location[2*j];p2[1]=location[2*j+1];
	 CORR[i*nn+j]=spheric(p1,p2);
	 if (i==j) CORR[i*nn+j]=1;
         else CORR[i*nn+j]=matern(CORR[i*nn+j],theta,alpha);
       }
   }
 inverse(CORR,CORRinverse,nn);
 prod(txy,CORRinverse,out1,pp+1,nn,nn);
 prod(out1,xy,out,pp+1,nn,pp+1);
 
 SS1=new double[pp];
 SS2=new double[pp];
 SSout=new double[pp*pp];
 for(k=2;k<pp+1;k++)
   {for(i=0;i<pp;i++)
     for(j=0;j<pp;j++)
       {if((i<k)&&(j<k)) SSout[i*pp+j]=out[i*(pp+1)+j];
       else if ((i>=k)&&(j<k)) SSout[i*pp+j]=out[(i+1)*(pp+1)+j];
       else if ((i<k)&&(j>=k)) SSout[i*pp+j]=out[i*(pp+1)+j+1];
       else if ((i>=k)&&(j>=k)) SSout[i*pp+j]=out[(i+1)*(pp+1)+j+1];
       }
   ls(SSout,pp);
   SS2[k-1]=SSout[0];
   }
 for(i=0;i<pp;i++)
   {sweep(out,i+1,pp+1);
   SS1[i]=out[0];
   }
 for(i=0;i<pp-1;i++)
   {SS1[i]=SS1[i]-SS1[i+1];
   }
 SST=SS1[0];dfT=nn-1;SSE=SS1[pp-1];dfE=nn-pp;

 for(i=0;i<pp+1;i++)
   for(j=0;j<pp+1;j++)
    if ((i!=0)&&(j!=0)) out[i*(pp+1)+j]=-out[i*(pp+1)+j];

 betacondtheta(yy,xx,location,theta,alpha,nn,pp,Liketemp);
 output[0]=theta[0];output[1]=theta[1];output[2*(pp+1)]=sqrt(out[0]/(nn-pp));
 for(i=0;i<pp;i++)
   {
     output[2+i]=out[i+1];
     output[2+pp+i]=output[2*(pp+1)]*sqrt(out[(i+1)*(pp+1)+(i+1)]);
   }

 yhat=new double[nn];
 rsort=new double[nn];
 for(i=0;i<nn;i++)
   {yhat[i]=0;residual[i]=0;
   for(j=0;j<pp;j++)
     yhat[i]+=xx[i*pp+j]*output[2+j];
   residual[i]=yy[i]-yhat[i];
   }
 mysort(residual,rsort,nn);

 outfile.open("coeff.krig",ios::out);
 outfile <<"Estimate"<<'\t'<<"Std"<<'\t'<<"t-value"<<'\t'<<"p-value"<<endl;
 outfile <<"Matern Correlation Function Used\n";
 outfile <<"Smooth Paramete is 1"<<endl;
 outfile <<"Estimate of Correlation Parameter is:"<<output[0]<<endl;
 outfile <<"Estimate of Scale Parameter is: "<<output[1]<<endl;

 for(i=0;i<pp;i++)
   {tvalue=output[2+i]/output[2+pp+i];
    pvalue=2*(1-pt(fabs(tvalue),nn-pp));
    outfile <<output[2+i]<<'\t'<<output[2+pp+i]<<'\t'<<tvalue<<'\t'<<pvalue<<endl;
   } 
 outfile.close();

 outfile.open("mse.krig",ios::out);
 outfile <<output[2*pp+2]*output[2*pp+2]<<endl;
 outfile.close();

 outfile.open("aicbic.krig",ios::out);
 outfile <<"AIC"<<'\t'<<"BIC"<<endl;
 outfile <<-2*ff[0]+2*pp<<'\t'<<-2*ff[0]+pp*log(nn*1.0);
 outfile.close();

 outfile.open("residualplot.krig",ios::out);
 outfile <<"predict"<<'\t'<<"residual"<<endl;
 for(i=0;i<nn;i++)
   outfile <<yhat[i]<<'\t'<<residual[i]<<endl;
 outfile.close();

 outfile.open("qqplot.krig",ios::out);
 outfile <<"normal-quantile"<<'\t'<<"residual"<<endl;
 for(i=0;i<nn;i++)
   outfile <<qnorm((i+1.0)/(nn+1.0))<<'\t'<<rsort[i]<<endl;
 outfile.close();

 outfile.open("anova.krig",ios::out);
 outfile <<"\tdf\tSS\tMS\tFvalue\tPvalue\n";
 for(i=0;i<pp-1;i++)
   {fvalue=SS1[i]/(SSE/dfE);
    pvalue=1-pf(fvalue,1,dfE);
    outfile <<"NO."<<i+1<<"\t"<<1<<'\t'<<SS1[i]<<'\t'<<SS1[i]<<'\t'<<fvalue<<'\t'<<pvalue<<endl;
   }
 outfile <<"Error\t"<<dfE<<'\t'<<SSE<<'\t'<<SSE/dfE<<endl;
 outfile <<"Total\t"<<dfT<<'\t'<<SST<<endl;
 outfile.close();

 outfile.open("estimate.krig",ios::out);
 outfile <<nn<<'\t'<<pp<<endl;
 outfile <<theta[0]<<'\t'<<theta[1]<<'\t'<<alpha<<endl;
 for(i=0;i<pp+1;i++)
   {for(j=0;j<pp+1;j++)
     outfile <<out[i*(pp+1)+j]<<'\t';
   outfile <<endl;
   }
 outfile.close();

 delete []xy;
 delete []txy;
 delete []CORR; 
 delete []CORRinverse;
 delete []SS1;
 delete []SS2;
 delete []SSout;
 delete []hh;
 delete []Liketemp;
 delete []indexmax;
 delete []yhat;
 delete []rsort;
 delete []out1;
 delete []out;
 return k;
}

/* krig prediction given theta and beta 
 input variable: yy, nn-dimensional response vector of observations
                 xx, nn*pp dimensional independent variable of observations
           location, nn*2 dimensional (longitude,latitude) for observations
                xx0, independent vector of the unobserved location
          location0, (longitude,latitude) for unobserved location
               beta, estimate of coefficients of independent variables,
                     their std, and MSE 
              theta, 2-dimensional vecotr for matern correlation function
              alpha, smooth parameter for matern correlation fnction
                nn,  row-dimension of xx 
                pp, column-dimension of xx
             output, 2-dimensional vector, the first is estimate, 
                     the second is the standard error of the estimate
                 
*/
double krig(double *yy,double *xx,double *location,double *xx0,double *location0,double *beta,double *theta,double alpha,int nn,int pp,double *output)
{int i,j;
  double S,p1[2],p2[2],*c0,*CORR,*CORRinverse,*error;
 
 CORR=new double[nn*nn];
 c0=new double[nn];
 for(i=0;i<nn;i++)
   {p1[0]=location[2*i];p1[1]=location[2*i+1];
    c0[i]=spheric(location0,p1);
    c0[i]=matern(c0[i],theta,alpha);
    for(j=0;j<nn;j++)
       {p2[0]=location[2*j];p2[1]=location[2*j+1];
	 CORR[i*nn+j]=spheric(p1,p2);
	 if (i==j) CORR[i*nn+j]=1;
         else CORR[i*nn+j]=matern(CORR[i*nn+j],theta,alpha);
       }
   }

 CORRinverse=new double[nn*nn];
 inverse(CORR,CORRinverse,nn);

 S=0; output[0]=0;
 for(j=0;j<pp;j++)
   {S+=xx0[j]*beta[j];
    output[0]+=xx0[j]*beta[j];
   }

 error=new double[nn];
 for(i=0;i<nn;i++)
   error[i]=yy[i];

 for(i=0;i<nn;i++)
   for(j=0;j<pp;j++)
     error[i]+=-xx[i*pp+j]*beta[j];
 
 output[1]=1;
 for(i=0;i<nn;i++)
   for(j=0;j<nn;j++)
     {S+=c0[i]*CORRinverse[i*nn+j]*error[j];
      output[0]+=c0[i]*CORRinverse[i*nn+j]*error[j];
      output[1]+=-(c0[i]*CORRinverse[i*nn+j]*c0[j]);
     }

 output[1]=beta[2*pp]*sqrt(output[1]);
 delete []c0;
 delete []CORR;
 delete []CORRinverse;
 delete []error;
 return S;
}

/*
output the column and row of a data file
*/
void colrow(char *filename,int *cols,int *rows)
{
		char str[4096], *ptr = NULL;
		*cols =1; *rows =1;
		ifstream f (filename);
		f.getline(str, 4096);
		ptr = str;
		while (ptr && (ptr = strstr(ptr, "\t"))) {
				(*cols)++; ptr++;
		}
		while (f.getline(str, 4096)) {
				(*rows)++;
		}
		f.close();
}

/*
regression and kriging estimation of the models
*/
int estimation(char *filename)
{int cols,rows,i,j,k,nn,pp,index;
 char line[256];
 ifstream obsfile;
 double *id,*time,*location,*yobs,*xobs,*beta,*estimate,*residual,*hh,*theta;
 double pvalue,alpha=1.0;

 colrow(filename,&cols,&rows);

 nn=rows-1;pp=cols-5; 
 obsfile.open(filename);
 index=obsfile.fail();
 if (index==1) 
   {cout <<"Observation file is not available\n"; return index;}
  
 for(i=0;i<pp+5;i++)
   {obsfile >>line;
   }

  id=new double[nn];
  time=new double[nn];
  location=new double[2*nn];
  yobs=new double[nn];
  xobs=new double[nn*pp];
  beta=new double[2*pp+1];
  residual=new double[nn];

  for(i=0;i<nn;i++)
    {obsfile >>id[i];
     obsfile >>time[i];
     obsfile >>location[2*i]>>location[2*i+1];
     obsfile >>yobs[i];
     for(j=0;j<pp;j++)
       obsfile >>xobs[i*pp+j];
    }

 estimate=new double[2*pp+1];
 pvalue=regression(yobs,xobs,location,nn,pp,estimate,residual);
 delete []estimate;

 if (pvalue<0.05)
  {estimate=new double[2*pp+3];
   theta=new double[2];
   for(i=0;i<2*pp+3;i++)
     estimate[i]=0;
   k=mle(yobs,xobs,location,alpha,nn,pp,estimate,residual);
   theta[0]=estimate[0];
   theta[1]=estimate[1];

  delete []estimate;
  delete []theta;
  }

 delete []id;
 delete []time;
 delete []location;
 delete []yobs;
 delete []xobs;
 delete []beta;
 delete []residual;
 obsfile.close();
 return index;
}


/*
regression and kriging prediction given station data, roast data and estimated values
name of station data: obs.txt
name of roast data: roast.txt
name of estimated values: estimate.reg (regression method)
                          estimate.krig (kriging method)
output file: prediction.txt (the last two columns are predicted values)
*/

int preddata(char *station,char *roast,char *estreg,char *estkrig,char *moranreg,char *output)
{ifstream stationfile,roastfile,estregfile,estkrigfile,moranfile;
 ofstream outputfile;
 char names[20];
 int rows,cols,i,j,nn,pp,nnkrig,ppkrig,nnstation,ppstation,index;
 double *id,*time,*yy,*xx,*xxinv,*xrxinv,*xrxrxinv,*xxtemp1,*xxtemp2,*xxtemp3,*xxtemp4,*location,id0,time0,*yy0,*xx0,*location0,*theta,alpha,*resultreg,*resultkrig,pvaluereg;
 double *betareg,*betakrig,*covreg,*covkrig,msereg,msekrig,*c0,*CORR,*CORRinverse,*error,*p1,*p2;

 moranfile.open(moranreg);
 index=moranfile.fail();
 if (index==1) 
   {cout <<"Moran I file of regression is not available\n"; return index;}

 for(i=0;i<5;i++)
   {moranfile >>names;
   }
  for(i=0;i<5;i++)
   {moranfile >>pvaluereg;
   }
 moranfile.close();

 estregfile.open(estreg);
 index=estregfile.fail();
 if (index==1) 
   {cout <<"Input regression estimation file is not available\n"; return index;}
 estregfile >>nn>>pp;
 resultreg=new double[(pp+1)*(pp+1)];
 betareg=new double[pp];
 covreg=new double[pp*pp];
 for(i=0;i<pp+1;i++)
   for(j=0;j<pp+1;j++)
     estregfile >>resultreg[i*(pp+1)+j];
 msereg=resultreg[0]/(nn-pp);
 for(i=0;i<pp;i++)
   betareg[i]=resultreg[i+1];
 for(i=0;i<pp;i++)
   for(j=0;j<pp;j++)
     covreg[i*pp+j]=resultreg[(i+1)*(pp+1)+j+1];
 estregfile.close();

 colrow(station,&cols,&rows);
 nnstation=rows-1;ppstation=cols-5;
 index=ppstation-pp;

 if (index!=0)
  {cout <<"Input station file does not match\n"; return index;}

 stationfile.open(station);
 index=stationfile.fail();
 if (index==1) 
   {cout <<"Input station dataset is not available\n"; return index;}

 id=new double[nnstation];
 time=new double[nnstation];
 location=new double[2*nnstation];
 yy=new double[nnstation];
 xx=new double[nnstation*ppstation]; 

 for(i=0;i<ppstation+5;i++)
   stationfile >>names;
 for(i=0;i<nnstation;i++)
   {stationfile >>id[i]>>time[i]>>location[2*i]>>location[2*i+1]>>yy[i];
   for(j=0;j<pp;j++)
     stationfile >>xx[i*pp+j];
   }
 stationfile.close();

 xxinv=new double[ppstation*ppstation];
 xrxinv=new double[ppstation*ppstation];
 xrxrxinv=new double[nnstation*nnstation];

 xxtemp1=new double[ppstation*ppstation];
 xxtemp2=new double[ppstation*nnstation];
 xxtemp3=new double[nnstation*ppstation];
 transpose(xx,xxtemp2,nnstation,ppstation);
 prod(xxtemp2,xx,xxtemp1,ppstation,nnstation,ppstation);
 inverse(xxtemp1,xxinv,ppstation);

 delete [] xxtemp1;
 delete [] xxtemp2;
 delete [] xxtemp3;

 if (pvaluereg<0.05)
  {theta=new double[2];
   estkrigfile.open(estkrig);
   index=estkrigfile.fail();
   if (index==1) 
    {cout <<"Input kriging estimation file is not available\n"; return index;}
  estkrigfile >>nnkrig>>ppkrig;
  estkrigfile >>theta[0]>>theta[1]>>alpha;
  resultkrig=new double[(ppkrig+1)*(ppkrig+1)];
  betakrig=new double[ppkrig];
  covkrig=new double[ppkrig*ppkrig];
  for(i=0;i<ppkrig+1;i++)
    for(j=0;j<ppkrig+1;j++)
      estkrigfile >>resultkrig[i*(ppkrig+1)+j];
  msekrig=resultkrig[0]/(nnkrig-ppkrig);
  for(i=0;i<ppkrig;i++)
    betakrig[i]=resultkrig[i+1];
  for(i=0;i<ppkrig;i++)
    for(j=0;j<ppkrig;j++)
      covkrig[i*ppkrig+j]=resultkrig[(i+1)*(ppkrig+1)+j+1];
  estkrigfile.close();

  p1=new double[2];p2=new double[2];
  CORR=new double[nnstation*nnstation];
  CORRinverse=new double[nnstation*nnstation];
  error=new double[nnstation];
  for(i=0;i<nnstation;i++)
    {p1[0]=location[2*i];p1[1]=location[2*i+1];
     for(j=0;j<nnstation;j++)
        {p2[0]=location[2*j];p2[1]=location[2*j+1];
 	 CORR[i*nnstation+j]=spheric(p1,p2);
	  if (i==j) CORR[i*nnstation+j]=1;
          else CORR[i*nnstation+j]=matern(CORR[i*nnstation+j],theta,alpha);
        }
    }
   inverse(CORR,CORRinverse,nnstation);

   xxtemp1=new double[ppstation*ppstation];
   xxtemp2=new double[ppstation*nnstation];
   xxtemp3=new double[ppstation*nnstation];
   xxtemp4=new double[nnstation*ppstation];
   
   transpose(xx,xxtemp2,nnstation,ppstation);
   prod(xxtemp2,CORRinverse,xxtemp3,ppstation,nnstation,nnstation);
   prod(xxtemp3,xx,xxtemp1,ppstation,nnstation,ppstation);
   transpose(xxtemp3,xxtemp4,ppstation,nnstation);
   inverse(xxtemp1,xrxinv,ppstation);      

   prod(xrxinv,xxtemp3,xxtemp2,ppstation,ppstation,nnstation);
   prod(xxtemp4,xxtemp2,xrxrxinv,nnstation,ppstation,nnstation);

   delete [] xxtemp1;
   delete [] xxtemp2;
   delete [] xxtemp3;
   delete [] xxtemp4;
  }

 roastfile.open(roast);
 index=roastfile.fail();
 if (index==1) 
   {cout <<"Input roast dataset is not available\n"; return index;}
 outputfile.open("prediction.txt",ios::out);
 location0=new double[2];
 xx0=new double[pp];
 for(i=0;i<pp+4;i++)
   {roastfile >>names;
   outputfile <<names<<'\t';
   }
 if (pvaluereg<0.05) outputfile <<"predreg\tstdreg\tpredkrig\tstdkrig\n"; 
 else outputfile <<"predreg\tstdreg\n"; 

 yy0=new double[4];
 c0=new double[nnstation];

 //loop for roast data below
 while(!roastfile.eof())
   {roastfile >>id0>>time0>>location0[0]>>location0[1];
    for(i=0;i<pp;i++)
      roastfile >>xx0[i];

   outputfile <<id0<<'\t'<<time0<<'\t'<<location0[0]<<'\t'<<location0[1]<<'\t';
     for(i=0;i<pp;i++)
       outputfile <<xx0[i]<<'\t';

    yy0[0]=yy0[1]=yy0[2]=yy0[3]=0.0;
    for(i=0;i<pp;i++)
      yy0[0]+=betareg[i]*xx0[i];

    for(i=0;i<pp;i++)
      for(j=0;j<pp;j++)
	yy0[2]+=xx0[i]*xxinv[i*pp+j]*xx0[j];
	yy0[2]=(1+yy0[2])*msereg;

	if (pvaluereg>=0.05)
	 {
	   if (!roastfile.eof()) outputfile <<yy0[0]<<'\t'<<sqrt(yy0[2])<<endl;
      else outputfile <<yy0[0]<<'\t'<<sqrt(yy0[2]);
     }
    else{
    for(i=0;i<pp;i++)
      yy0[1]+=betakrig[i]*xx0[i]; 
      
	 for(i=0;i<nnstation;i++)
       {p1[0]=location[2*i];p1[1]=location[2*i+1];
       c0[i]=spheric(location0,p1);
       c0[i]=matern(c0[i],theta,alpha);
	   }
     for(i=0;i<nnstation;i++)
       error[i]=yy[i];       
     for(i=0;i<nnstation;i++)
       for(j=0;j<pp;j++)
	 error[i]+=-xx[i*pp+j]*betakrig[j];
     for(i=0;i<nnstation;i++)
       for(j=0;j<nnstation;j++)
	 yy0[1]+=c0[i]*CORRinverse[i*nnstation+j]*error[j];

	 //old kriging prediction standard error below
	 yy0[3]=1.0;
     for(i=0;i<nnstation;i++)
       for(j=0;j<nnstation;j++)
         yy0[3]+=-c0[i]*CORRinverse[i*nnstation+j]*c0[j];
   
     yy0[3]=yy0[3]*msekrig;

     if (!roastfile.eof()) outputfile <<yy0[0]<<'\t'<<sqrt(yy0[2])<<'\t'<<yy0[1]<<'\t'<<sqrt(yy0[3])<<endl;
     else outputfile <<yy0[0]<<'\t'<<sqrt(yy0[2])<<'\t'<<yy0[1]<<'\t'<<sqrt(yy0[3]);
	}
   }

 index=0;
 roastfile.close();
 outputfile.close();

 delete []id;
 delete []time;
 delete []yy;
 delete []xx;
 delete []xxinv;
 delete []xrxinv;
 delete []xrxrxinv;
 delete []location;
 delete []yy0;
 delete []xx0;
 delete []location0;
 delete []theta;
 delete []resultreg;
 delete []resultkrig;
 delete []betareg;
 delete []betakrig;
 delete []covreg;
 delete []covkrig;
 delete []c0;
 delete []CORR;
 delete []CORRinverse;
 delete []error;
 delete []p1;
 delete []p2;
 return index;
}

/*
cross validation for kriging method given data and the regression method 
*/
int cv(char *station,char *estreg,char *estkrig,char *moranreg,char *output)
{ifstream stationfile,estregfile,estkrigfile,moranfile;
 ofstream outputfile;
 char names[20];
 int i,j,k,nn,pp,nnkrig,ppkrig,index,rowindex;
 double *id,*time,*yy,*xx,*location,*theta,alpha,y0,y0hat,*yobs,*xobs,*llobs,*x0,*ll0,*beta,*outtemp,*cvresult,*xy,*txy,*out,*betareg,*betakrig,*yhat;

 double id0,time0,*yy0,*xx0,*location0,*resultreg,*resultkrig,pvaluereg;
 double *covreg,*covkrig,msereg,msekrig,*c0,*CORR,*CORRinverse,*error,*p1,*p2;

 moranfile.open(moranreg);
 index=moranfile.fail();
 if (index==1) 
   {cout <<"Moran I file of regression is not available\n"; return index;}

 for(i=0;i<5;i++)
   {moranfile >>names;
   }
  for(i=0;i<5;i++)
   {moranfile >>pvaluereg;
   }
 moranfile.close();

 estregfile.open(estreg);
 index=estregfile.fail();
 if (index==1) 
   {cout <<"Input regression estimation file is not available\n"; return index;}
 estregfile >>nn>>pp;
 resultreg=new double[(pp+1)*(pp+1)];
 betareg=new double[pp];
 covreg=new double[pp*pp];
 for(i=0;i<pp+1;i++)
   for(j=0;j<pp+1;j++)
     estregfile >>resultreg[i*(pp+1)+j];
 msereg=resultreg[0]/(nn-pp);
 for(i=0;i<pp;i++)
   betareg[i]=resultreg[i+1];
 for(i=0;i<pp;i++)
   for(j=0;j<pp;j++)
     covreg[i*pp+j]=resultreg[(i+1)*(pp+1)+j+1];
 estregfile.close();
 
 stationfile.open(station);
 index=stationfile.fail();
 if (index==1) 
   {cout <<"Input station dataset is not available\n"; return index;}
 id=new double[nn];
 time=new double[nn];
 location=new double[2*nn];
 yy=new double[nn];
 xx=new double[nn*pp]; 

 for(i=0;i<pp+5;i++)
   stationfile >>names;
 for(i=0;i<nn;i++)
   {stationfile >>id[i]>>time[i]>>location[2*i]>>location[2*i+1]>>yy[i];
   for(j=0;j<pp;j++)
     stationfile >>xx[i*pp+j];
   }
 stationfile.close();

 delete []betareg;

 yobs=new double[nn-1];
 xobs=new double[(nn-1)*pp];
 llobs=new double[(nn-1)*2];
 x0=new double[pp];
 ll0=new double[2];
 betakrig=new double[2*pp+1];
 betareg=new double[2*pp+1];
 outtemp=new double[2];
 cvresult=new double[4]; 
 theta=new double[2];
 xy=new double[(nn-1)*(pp+1)];
 txy=new double[(nn-1)*(pp+1)];
 out=new double[(pp+1)*(pp+1)];
 cvresult[0]=cvresult[1]=cvresult[2]=cvresult[3]=0.0;

 if (pvaluereg<0.05){
  estkrigfile.open(estkrig);
  index=estkrigfile.fail();
  if (index==1)
    {cout <<"Input kriging estimation file is not available\n"; return index;}
  estkrigfile >>nnkrig>>ppkrig;
  estkrigfile >>theta[0]>>theta[1]>>alpha;
  estkrigfile.close();
 }

 for(k=0;k<nn;k++){
   y0=yy[k];
   ll0[0]=location[2*k];ll0[1]=location[2*k+1];
   for(j=0;j<pp;j++)
     x0[j]=xx[k*pp+j];

   for(i=0;i<nn-1;i++)
     {if (i<k) rowindex=i; else rowindex=i+1;
     yobs[i]=yy[rowindex];
     llobs[2*i]=location[2*rowindex];llobs[2*i+1]=location[2*rowindex+1];
     for(j=0;j<pp;j++)
       xobs[i*pp+j]=xx[rowindex*pp+j];
     }

  for(i=0;i<nn-1;i++)
    {xy[i*(pp+1)]=yobs[i];
    for(j=0;j<pp;j++)
      xy[i*(pp+1)+j+1]=xobs[i*pp+j];
    }

  transpose(xy,txy,nn-1,pp+1);
  prod(txy,xy,out,pp+1,nn-1,pp+1);
  ls(out,pp+1);

  for(i=0;i<pp;i++)
    {betareg[i]=out[i+1];
     betareg[pp+i]=sqrt(output[2*pp]*out[(i+1)*(pp+1)+(i+1)]);
    }
  betareg[2*pp]=out[0]/(nn-pp);

  y0hat=0;
  for(i=0;i<pp;i++)
    y0hat+=betareg[i]*x0[i];
  cvresult[0]+=fabs(y0hat-y0)/nn;
  cvresult[1]+=(y0hat-y0)*(y0hat-y0)/nn;

   if(pvaluereg<0.05){
     betacondtheta(yobs,xobs,llobs,theta,alpha,nn-1,pp,betakrig);
     y0hat=krig(yobs,xobs,llobs,x0,ll0,betakrig,theta,alpha,nn-1,pp,outtemp);
     cvresult[2]+=fabs(y0hat-y0)/nn;
     cvresult[3]+=(y0hat-y0)*(y0hat-y0)/nn;
   }
 }

 outputfile.open("cv.txt",ios::out);
 outputfile <<"CVregL1\tCVregL2\tCVkrigL2\tCVkrigL2\n";
 outputfile <<cvresult[0]<<'\t'<<cvresult[1]<<'\t'<<cvresult[2]<<'\t'<<cvresult[3]<<endl;
 outputfile.close();
 
 delete []resultreg;
 delete []covreg;
 delete []id;
 delete []time;
 delete []location;
 delete []yy;
 delete []xx;
 delete []cvresult;
 delete []yobs;
 delete []xobs;
 delete []llobs;
 delete []x0;
 delete []ll0;
 delete []betakrig;
 delete []betareg;
 delete []outtemp;
 delete []theta;
 delete []xy;
 delete []txy;
 delete []out;
 return index;
}
