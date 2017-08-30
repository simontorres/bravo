#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>

using namespace std;

int nlin;





float pdm(double time[], double signal[], float f0, float fn, float df, int Nbin, int Ncover, float D, int rr,ofstream &dumFile){

 int k,j,l,m,f,ind,ntot;
 float frq;
 float phase,cuam,promag,theta,freqfin,ttot=9999;
 float quad[Nbin*Ncover];
 float prom[Nbin*Ncover];
 int num[Nbin*Ncover];
 float binsize = 1.0/Nbin ;
 float covershift = 1.0/(Nbin*Ncover) ;


 for(k = 0; k < rr; k++)
    {
//      frq=pow(10,f0+(fn-f0)/float(rr)*k);
      frq=1/(f0+(fn-f0)/float(rr)*k);
//      frq=(f0+(fn-f0)/float(rr)*k);
      cuam=0;
      promag=0;
      for(f=0;f<Nbin*Ncover;f++){quad[f]=0;prom[f]=0;num[f]=0;}
      theta=0;
      ntot=0;
      for(j = 0; j < nlin; j++)
	{
	cuam+=pow(signal[j],2);
	promag+=signal[j];
	phase=fmod((time[j] - time[0]) * frq + D/2.*pow(time[j],2), 1.0);
	for(l = 0; l < Ncover; l++)
	  {
	    ind=int(fmod((phase-covershift*l)/binsize+float(Nbin),float(Nbin)))+l*Nbin;
	    quad[ind]+=pow(signal[j],2);
	    prom[ind]+=signal[j];
	    num[ind]++;
	  }
	}
     for(m = 0; m < Nbin*Ncover; m++)
	{
	 if(num[m] > 0)
	 {
	 theta+=(quad[m]/num[m]-pow(prom[m]/num[m],2))*(num[m]-1);
	 ntot++;
	 }
	}
     theta=theta/(Ncover*nlin-ntot);
     theta=theta/((cuam/nlin-pow(promag/nlin,2)));
     dumFile << setprecision(10) << theta <<" "<< frq << "\n";
     if(theta<ttot && theta>0){ttot=theta;freqfin=frq;}
     }




 return 1/freqfin;

}





int main (int argc, char* argv[]) {
	FILE * datos;
        long taman;
        int i,n,m,nit;
        double b,u;
        datos = fopen (argv[1],"r");
        if (datos==NULL) {fputs ("File error",stderr); exit (1);}
        nlin= atoi(argv[2]);
        double jd[nlin],mag[nlin],mager[nlin];
	ofstream dumFile (string((argv[1])+string(".pdm")).c_str());

        for(n = 0; n < nlin; n++)
            {
              fscanf(datos,"%lE", &jd[n] );
              fscanf(datos,"%lE",&mag[n]);
	      fscanf(datos,"%lE",&mager[n]);
	      fscanf(datos,"%lE",&u);
	      fscanf(datos,"%lE",&u);
            }
//	cout << pdm(jd,mag,-1,0,0.0001,10,5,0,4000,dumFile) << "\n";
	cout << setprecision(10) << pdm(jd,mag,atof(argv[3]),atof(argv[4]),0.0001,atof(argv[5]),atof(argv[6]),0,atof(argv[7]),dumFile) << "\n";
	}
	
