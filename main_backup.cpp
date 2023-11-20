#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>
#include <cmath>
#include <Eigen/Dense>

#include "fftw3.h"

const int L=100;
const double beta = atanh(sqrt(2.0)-1.0);
const double beta_tilde = -0.5*log( tanh(beta) );

using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);

#include "header.h"

double epsilon_tilde_correlator(const int p, const TransferMatrix& T, const int t){
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
  for(int t=0; t<tmax; t++) std::cout << t << ' ' << epsilon_tilde_correlator( p, T, t ) << std::endl;

    {
    const Idx N = 100;

    fftw_complex *in, *out;
    fftw_plan p;
    // ...
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    //
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    //
    // initialize in
    //
    for(Idx i=0; i<N; ++i){
      in[i][0] = std::cos(2.0*M_PI*i/N);
      in[i][1] = std::sin(2.0*M_PI*i/N);
    }
    //
    fftw_execute(p); /* repeat as needed */
    for(Idx i=0; i<N; ++i){
      if(std::abs(out[i][0]<1.0e-14)) out[i][0] = 0.0;
      if(std::abs(out[i][1]<1.0e-14)) out[i][1] = 0.0;
      std::cout << "out[" << i << "] = " << out[i][0] << " + i* " << out[i][1] << std::endl;
    }
    // ...
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

  return 0;
}
