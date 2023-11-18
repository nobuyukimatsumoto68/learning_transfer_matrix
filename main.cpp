#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>
#include <cmath>
#include <Eigen/Dense>

const int L=100;
const double beta = atanh(sqrt(2.0)-1.0);
const double beta_tilde = -0.5*log( tanh(beta) );

using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);

#include "header.h"

double epsilon_correlator(const int p, const TransferMatrix& T, const int t){
  ProductState vac;
  vac.set2vacuum( T );

  GeneralState v0;
  v0.epsilon_tilde(p, vac);
  const Complex vev = v0.dot(vac)*(1.0/L);

  GeneralState vt = v0;
  vt.evolve(t, T);
  const Complex eet = v0.dot(vt) * (1.0/L/L);

  const Complex res = eet-vev*vev;
  assert( abs(res.imag()) < 1.0e-12 );
  return res.real();
}


int main(int argc, char* argv[])
{
  // prep
  std::cout << std::scientific << std::setprecision(15);
  TransferMatrix T;

  // main part
  int p = 0;
  if(argc==2){ p = atoi(argv[1]); }

  const int tmax = 24;
  std::cout << "# L = " << L << ", p = " << p << std::endl;
  for(int t=0; t<tmax; t++) std::cout << t << ' ' << epsilon_correlator( p, T, t) << std::endl;

  return 0;
}
