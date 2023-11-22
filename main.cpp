#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>
#include <cmath>
#include <Eigen/Dense>

#include "fftw3.h"

#ifdef _OPENMP
#include "omp.h"
#else
int omp_get_thread_num() { return 0; }
int omp_get_num_threads() { return 1; }
int omp_get_max_threads() { return 1; }
#endif

const int L=128;
const int tmax = 128;
std::string desc = "L"+std::to_string(L)+"_T"+std::to_string(tmax);
const std::string pcorr = "./data/fourier_"+desc+".dat";
const std::string xcorr = "./data/spatial_"+desc+".dat";

const double beta = atanh(sqrt(2.0)-1.0);
const double beta_tilde = -0.5*log( tanh(beta) );

const int nthread = 8;

using Complex = std::complex<double>;
const Complex I = Complex(0.0, 1.0);
using Idx = long unsigned int;

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

void FFT( const uint N, const Complex* in_, Complex* out_ ){

  fftw_complex *in, *out;
  fftw_plan p;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  for(Idx i=0; i<N; ++i){
    in[i][0] = in_[i].real();
    in[i][1] = in_[i].imag();
  }

  fftw_execute(p); /* repeat as needed */

  for(Idx i=0; i<N; ++i){
    out_[i] = out[i][0] + I*out[i][1];
  }

  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}

void FFT2real( const uint N, const Complex* in_, double* out_ ){

  fftw_complex *in, *out;
  fftw_plan p;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  for(Idx i=0; i<N; ++i){
    in[i][0] = in_[i].real();
    in[i][1] = in_[i].imag();
  }

  fftw_execute(p); /* repeat as needed */

  for(Idx i=0; i<N; ++i){
    assert( abs(out[i][1])<1.0e-14 );
    out_[i] = out[i][0];
  }

  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}


int main(int argc, char* argv[])
{
  // prep
  std::cout << std::scientific << std::setprecision(15);
  TransferMatrix T;
#ifdef _OPENMP
  omp_set_num_threads(nthread);
#endif

  Complex corr_tilde[tmax][L];
  double corr[tmax][L];

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int t=0; t<tmax; t++) {
    for(int p=0; p<L; p++) corr_tilde[t][p] = epsilon_tilde_correlator( p, T, t );
    FFT2real(L, corr_tilde[t], corr[t]);
  }

  write2file2d<tmax, L>( corr_tilde, pcorr, "t\tk\tcorr (Re,Im)" );
  write2file2d<tmax, L>( corr, xcorr, "t\tx\tcorr" );

  return 0;
}
