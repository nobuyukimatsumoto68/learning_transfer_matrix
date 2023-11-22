#pragma once

#include <iomanip>
#include <complex>
#include <sstream>
#include <fstream>
#include <vector>
#include <cassert>
#include <initializer_list>
#include <string>
#include <cmath>
#include <Eigen/Dense>

using CM2 = Eigen::Matrix2cd;
using RV2 = Eigen::Vector2d;
using CV2 = Eigen::Vector2cd;

inline int mod(const int a, const int b, const int c=0){ return ((a+c)%b+b)%b-c; }

struct TransferMatrix {
  CM2 id, sigma1, sigma2, sigma3;

  std::vector<CM2> evec_TE; // | 0 >, |k,-k> sector
  std::vector<RV2> eval_TE;
  //
  std::vector<CM2> evec_TO; // | k >, |-k > sector
  std::vector<RV2> eval_TO;
  //
  CM2 evec_T0;
  RV2 eval_T0;
  //
  CM2 evec_TL;
  RV2 eval_TL;

  void initialize_sigma(){
    id << 1.0, 0.0, 0.0, 1.0;
    sigma1 << 0.0, 1.0, 1.0, 0.0;
    sigma2 << 0.0, -I, I, 0.0;
    sigma3 << 1.0, 0.0, 0.0, -1.0;
  }

  TransferMatrix()
    : evec_TE(L)
    , evec_TO(L)
    , eval_TE(L)
    , eval_TO(L)
  {
    assert(L%2==0);
    initialize_sigma();
    compute_eigensystems();
  }

  // | 0 >, |k,-k> sector
  CM2 thetaE(const int k){
    CM2 res = cosh(2.0*beta)*id + sinh(2.0*beta)*( cos(M_PI*k/L)*sigma3 + sin(M_PI*k/L)*sigma2 );
    res *= exp(-2.0*beta*cos(M_PI*k/L) );
    return res;
  }

  CM2 thetaE_half(const int k){
    CM2 res = cosh(beta)*id + sinh(beta)*( cos(M_PI*k/L)*sigma3 + sin(M_PI*k/L)*sigma2 );
    res *= exp(-beta*cos(M_PI*k/L) );
    return res;
  }

  CM2 theta_tildeE(const int k){
    CM2 res = cosh(2.0*beta_tilde)*id - sinh(2.0*beta_tilde)*sigma3;
    res /= 2.0*sinh(2.0*beta_tilde);
    return res;
  }

  CM2 E(const int k){
    assert(k%L!=0);
    const CM2 theta_half = thetaE_half(k);
    const CM2 theta_tilde = theta_tildeE(k);
    return theta_half*theta_tilde*theta_half;
  }

  // | k >, |-k > sector
  CM2 thetaO(const int k){
    CM2 res = id;
    res *= exp(-2.0*beta*cos(M_PI*k/L) );
    return res;
  }

  CM2 thetaO_half(const int k){
    CM2 res = id;
    res *= exp(-beta*cos(M_PI*k/L) );
    return res;
  }

  CM2 theta_tildeO(const int k){
    CM2 res = id;
    res /= 2.0*sinh(2.0*beta_tilde);
    return res;
  }

  CM2 O(const int k){
    assert(k%L!=0);
    const CM2 theta_half = thetaO_half(k);
    const CM2 theta_tilde = theta_tildeO(k);
    return theta_half*theta_tilde*theta_half;
  }

  CM2 Zero() {
    CM2 res;
    res << 1.0, 0.0, 0.0, exp(-2.0*beta + 2.0*beta_tilde);
    return res;
  }

  CM2 Ell() {
    CM2 res;
    res << 1.0, 0.0, 0.0, exp(2.0*beta + 2.0*beta_tilde);
    res /= 2.0*exp(2.0*beta_tilde)*sinh(2.0*beta_tilde);
    return res;
  }

  // -------------


  void compute_eigensystems() {

    // TE; 2x2
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int k=1; k<L; k++){
      CM2 Tk = E(k);
      Eigen::SelfAdjointEigenSolver< CM2 > solver(Tk);
      evec_TE[k] = solver.eigenvectors().rowwise().reverse();
      eval_TE[k] = solver.eigenvalues().reverse();
    }

    // TO; 2x2
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int k=1; k<L; k++){
      CM2 Tk = O(k);
      Eigen::SelfAdjointEigenSolver< CM2 > solver(Tk);
      evec_TO[k] = solver.eigenvectors().rowwise().reverse();
      eval_TO[k] = solver.eigenvalues().reverse();
    }

    // Zero
    {
      CM2 Tk = Zero();
      Eigen::SelfAdjointEigenSolver< CM2 > solver(Tk);
      evec_T0 = solver.eigenvectors().rowwise().reverse();
      eval_T0 = solver.eigenvalues().reverse();
    }

    // Ell
    {
      CM2 Tk = Ell();
      Eigen::SelfAdjointEigenSolver< CM2 > solver(Tk);
      evec_TL = solver.eigenvectors().rowwise().reverse();
      eval_TL = solver.eigenvalues().reverse();
    }

  }


  // level=0: lower energy
  CV2 evec( const int k, const bool is_02sector, const int level ) const {
    CV2 res;
    if(k==0) res = evec_T0.col(level);
    else if(k==L) res = evec_TL.col(level);
    else if(is_02sector==1) res = evec_TE[k].col(level);
    else res = evec_TO[k].col(level);
    return res;
  }

  // level=0: lower energy
  double eval( const int k, const int is_02sector, const int level ) const {
    double res;
    if(k==0) res = eval_T0(level);
    else if(k==L) res = eval_TL(level);
    else if(is_02sector==1) res = eval_TE[k](level);
    else res = eval_TO[k](level);
    return res;
  }

};



struct KState {
  // 0, 2, +, -
  std::vector<Complex> coeffs;
  int parity; // 0,2=>+1, +,-=>-1, otherwise 0

  KState()
    : coeffs(4, 0.0)
    , parity(0.0)
  {}

  void set( const std::initializer_list<Complex> args, const int parity_ ) {
    auto itr = args.begin();
    for( int i=0; i<4; i++ ) coeffs[i] = *(itr+i);
    parity = parity_;
  }

  std::string print() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);
    for( Complex elem : coeffs ) ss << elem << " ";
    return ss.str();
  }

  void a_tilde(const int sign=1){
    parity *= -1;
    if(sign==1) coeffs = std::vector<Complex>{coeffs[2], 0.0, 0.0, coeffs[1]};
    else if(sign==-1) coeffs = std::vector<Complex>{coeffs[3], 0.0, -coeffs[1], 0.0};
    else assert(false);
  }

  void a_tilde_dagger(const int sign=1){
    parity *= -1;
    if(sign==1) coeffs = std::vector<Complex>{0.0, coeffs[3], coeffs[0], 0.0};
    else if(sign==-1) coeffs = std::vector<Complex>{0.0, -coeffs[2], 0.0, coeffs[0]};
    else assert(false);
  }

  void mult( const Complex scalar ){ for( int i=0; i<4; i++ ) coeffs[i] *= scalar; }

  Complex dot( const KState& other ) const {
    Complex res = 0.0;
    for(int i=0; i<4; i++) res += conj(coeffs[i])*other.coeffs[i];
    return res;
  }

  Complex squared_norm() const { return dot( *this ); }

  void evolve( const int t, const int k, const TransferMatrix& T ) {
    const double lambda0 = T.eval_TE[k][0];

    if( parity==1 ){
      const CM2 eigvecs = T.evec_TE[k];
      const RV2 eigvals = T.eval_TE[k];

      CV2 v;
      v << coeffs[0], coeffs[1];
      const RV2 lambda_t = (eigvals/lambda0).array().pow(t);
      v = eigvecs.adjoint() * v;
      v = lambda_t.asDiagonal() * v;
      v = eigvecs * v;
      coeffs[0] = v(0);
      coeffs[1] = v(1);
    }
    else if(parity==-1){
      const CM2 eigvecs = T.evec_TO[k];
      const RV2 eigvals = T.eval_TO[k];

      CV2 v;
      v << coeffs[2], coeffs[3];
      const RV2 lambda_t = (eigvals/lambda0).array().pow(t);
      v = eigvecs.adjoint() * v;
      v = lambda_t.asDiagonal() * v;
      v = eigvecs * v;
      coeffs[2] = v(0);
      coeffs[3] = v(1);
    }
    else assert(false);
  }

};


struct ProductState {
  std::vector<KState> k_states;
  bool parity;
  // parity=true: k=1,3,5,7,...
  // parity=false: k=[0L],2,4,6,...

  ProductState( const bool parity_=1 )
    : k_states(L/2)
    , parity(parity_)
  {}

  int kidx( const int k ) const {
    assert(k>=0);

    int res;
    if(parity==1) {
      assert(k%2==1);
      res = (k-1)/2;
    }
    else {
      assert(k%2==0);
      res = (k%L)/2-1;
    }

    return res;
  }

  KState& operator()(const int k) { return k_states[kidx(k)]; }
  KState operator()(const int k) const { return k_states[kidx(k)]; }

  std::string print() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);
    for( const KState& elem : k_states ) ss << elem.print() << std::endl;
    return ss.str();
  }

  void set2vacuum( const TransferMatrix& T ) {
    parity=true;
    for(int k=1; k<L; k+=2){
      const CV2 even = T.evec(k, true, 0);
      (*this)(k).set({even(0), even(1), 0.0, 0.0}, +1);
    }
  }

  void set2empty( const bool is_02sector=1 ) {
    parity=true;
    for(int k=1; k<L; k+=2)(*this)(k).set({1.0, 0.0, 0.0, 0.0}, +1);
  }

  void a_tilde( const int k_ ){
    const int k=mod(k_,2*L,L);
    const int absk = abs(k);
    const int signk = k/absk;
    (*this)(absk).a_tilde(signk);

    int sign=1;
    for(int ik=0; ik<kidx(absk); ik++) sign*=k_states[ik].parity;
    (*this)(absk).mult(sign);
  }

  void a_tilde_dagger( const int k_ ){
    const int k=mod(k_,2*L,L);
    const int absk = abs(k);
    const int signk = k/absk;
    (*this)(absk).a_tilde_dagger(signk);

    int sign=1;
    for(int ik=0; ik<kidx(absk); ik++) sign*=k_states[ik].parity;
    (*this)(absk).mult(sign);
  }

  Complex dot( const ProductState& other ) const {
    Complex res = 1.0;
    for(int ik=0; ik<k_states.size(); ik++) res *= k_states[ik].dot( other.k_states[ik] );
    return res;
  }

  void evolve( const int t, const TransferMatrix& T ) {
    if(parity) for(int k=1; k<=L-1; k+=2) (*this)(k).evolve(t, k, T);
    else assert(false); // to be implemented
  }

  Complex squared_norm() const { return dot( *this ); }

};


struct GeneralState {
  std::vector<ProductState> states;
  std::vector<Complex> coeffs;

  GeneralState()
    : states()
    , coeffs()
  {}

  GeneralState( const ProductState& state )
    : states({state})
    , coeffs(1,1.0)
  {}

  std::string print() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);
    for( int i=0; i<states.size(); i++ ) {
      ss << coeffs[i] << std::endl;
      ss << states[i].print() << std::endl;
    }
    return ss.str();
  }

  void a_tilde( const int k, const GeneralState& other ){
    states.clear();
    coeffs.clear();

    for(int i=0; i<other.states.size(); i++){
      ProductState tmp = other.states[i];
      tmp.a_tilde(k);

      states.push_back(tmp);
      coeffs.push_back(other.coeffs[i]);
    }
  }

  void a_tilde_dagger( const int k, const GeneralState& other ){
    states.clear();
    coeffs.clear();

    for(int i=0; i<other.states.size(); i++){
      ProductState tmp = other.states[i];
      tmp.a_tilde_dagger(k);

      states.push_back(tmp);
      coeffs.push_back(other.coeffs[i]);
    }
  }

  void epsilon_tilde( const int p, const ProductState& state ){
    states.clear();
    coeffs.clear();

    for(int k=-L+1; k<=L-1; k+=2){
      {
        ProductState tmp = state;
        tmp.a_tilde(2.0*p-k);
        tmp.a_tilde(k);
        const double mom_phys = M_PI*(2.0*p-k)/L;
        states.push_back( tmp );
        coeffs.push_back( exp(-I*mom_phys) );
      }
      {
        ProductState tmp = state;
        tmp.a_tilde_dagger(-2.0*p+k);
        tmp.a_tilde(k);
        const double mom_phys = M_PI*(2.0*p-k)/L;
        states.push_back( tmp );
        coeffs.push_back( exp(-I*mom_phys) );
      }
      {
        ProductState tmp = state;
        tmp.a_tilde(2.0*p+k);
        tmp.a_tilde_dagger(k);
        const double mom_phys = M_PI*(2.0*p+k)/L;
        states.push_back( tmp );
        coeffs.push_back( -exp(-I*mom_phys) );
      }
      {
        ProductState tmp = state;
        tmp.a_tilde_dagger(-2.0*p-k);
        tmp.a_tilde_dagger(k);
        const double mom_phys = M_PI*(2.0*p+k)/L;
        states.push_back( tmp );
        coeffs.push_back( -exp(-I*mom_phys) );
      }
    }
  }



  // void sigma( const int x, const ProductState& state ){
  //   states.clear();
  //   coeffs.clear();

  //   for(int k=-L+1; k<=L-1; k+=2){
  //     {
  //       ProductState tmp = state;
  //       tmp.a_tilde(2.0*p-k);
  //       tmp.a_tilde(k);
  //       const double mom_phys = M_PI*(2.0*p-k)/L;
  //       states.push_back( tmp );
  //       coeffs.push_back( exp(-I*mom_phys) );
  //     }
  //     {
  //       ProductState tmp = state;
  //       tmp.a_tilde_dagger(-2.0*p+k);
  //       tmp.a_tilde(k);
  //       const double mom_phys = M_PI*(2.0*p-k)/L;
  //       states.push_back( tmp );
  //       coeffs.push_back( exp(-I*mom_phys) );
  //     }
  //     {
  //       ProductState tmp = state;
  //       tmp.a_tilde(2.0*p+k);
  //       tmp.a_tilde_dagger(k);
  //       const double mom_phys = M_PI*(2.0*p+k)/L;
  //       states.push_back( tmp );
  //       coeffs.push_back( -exp(-I*mom_phys) );
  //     }
  //     {
  //       ProductState tmp = state;
  //       tmp.a_tilde_dagger(-2.0*p-k);
  //       tmp.a_tilde_dagger(k);
  //       const double mom_phys = M_PI*(2.0*p+k)/L;
  //       states.push_back( tmp );
  //       coeffs.push_back( -exp(-I*mom_phys) );
  //     }
  //   }
  // }



  // void a_tilde( const int k ){
  //   const GeneralState other = *this;
  //   a_tilde( k, other );
  // }

  // void a_tilde_dagger( const int k ){
  //   const GeneralState other = *this;
  //   a_tilde_dagger( k, other );
  // }

  // void epsilon_tilde( const int p ){
  //   const ProductState other = *this;
  //   epsilon_tilde( p, other );
  // }

  void evolve( const int t, const TransferMatrix& T ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<states.size(); i++) states[i].evolve(t, T);
  }

  Complex dot( const ProductState& state ) const {
    Complex res = 0.0;
    for(int i=0; i<states.size(); i++) res += conj(coeffs[i]) * states[i].dot(state);
    return res;
  }

  Complex dot( const GeneralState& other ) const {
    Complex res = 0.0;
    for(int i=0; i<states.size(); i++) {
      for(int j=0; j<other.states.size(); j++) {
        res += conj(coeffs[i])*other.coeffs[j] * states[i].dot(other.states[j]);
      }}
    return res;
  }

  Complex squared_norm() const { return dot( *this ); }

};



template <int X, int Y>
void write2file2d( const double obj[X][Y], const std::string& filename, const std::string& memo="" ){
  std::ofstream ofs(filename);
  ofs << std::scientific << std::setprecision(15);
  ofs << "# " << memo << std::endl;
  for(int x=0; x<X; x++){
    for(int y=0; y<Y; y++) {
      ofs << x << '\t' << y << '\t' << obj[x][y] << std::endl;
    }
  }
}

template <int X, int Y>
void write2file2d( const Complex obj[X][Y], const std::string& filename, const std::string& memo="" ){
  std::ofstream ofs(filename);
  ofs << std::scientific << std::setprecision(15);
  ofs << "# " << memo << std::endl;
  for(int x=0; x<X; x++){
    for(int y=0; y<Y; y++) {
      ofs << x << '\t' << y << '\t' << obj[x][y].real() << '\t' << obj[x][y].imag() << std::endl;
    }
  }
}


