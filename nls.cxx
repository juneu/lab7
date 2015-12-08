#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void function(const double eta, double* y, double* f);

int main(){

  const double eta = 0.5;
  const int xmin = 0;
  const int xmax = 100;
  const double dx = 0.0001;
  double x = xmin;
  const int Nit = (xmax-xmin)/dx + 1;
  
  const int N = 2;
  double y[N];
  double f[N];
  double k1[N];
  double k2[N];
  double k3[N];
  double ytemp[N];
  
  //initial values
    y[0] = 1e-5;
    y[1] = sqrt(eta) * y[0];
  
  ofstream out ("lab7_output_dx_0.0001.dat");
  
  out << x << "\t" << y[0] << endl;
  
  for(int i=0; i<Nit; i++){

    x += dx;
    
    function(eta,y,f);
    k1[0] = f[0];
    k1[1] = f[1];
  
    ytemp[0] = y[0] + dx/2*k1[0];
    ytemp[1] = y[1] + dx/2*k1[1];
  
    function(eta,ytemp,f);
    k2[0] = f[0];
    k2[1] = f[1];
    
    ytemp[0] = y[0] + dx*(-k1[0] + 2*k2[0]);
    ytemp[1] = y[1] + dx*(-k1[1] + 2*k2[1]);
    
    function(eta,ytemp,f);
    k3[0] = f[0];
    k3[1] = f[1];
    
    y[0] += dx/6 * (k1[0] + 4*k2[0] + k3[0]);
    y[1] += dx/6 * (k1[1] + 4*k2[1] + k3[1]);
    
    out << x << "\t" << y[0] << endl;
    
  }
  
  out.close();
  
  return 0;
}

void function(const double eta, double* y, double* f){

  f[0] = y[1];
  f[1] = y[0]*(eta - y[0]*y[0]);
  
}
