
  // GeneralState v0;
  // v0.epsilon_tilde(0, vac);
  // Complex vev = v0.dot(vac)*(1.0/L);
  // std::cout << vev << std::endl;
  // // GeneralState v(vac);
  // // // const int k=1;
  // // const int t=1;
  // // ProductState state1 = vac;
  // // std::cout << "vac: " << state1.print() << std::endl;
  // // state1.a_tilde_dagger(1);
  // // state1.a_tilde_dagger(-1);
  // // std::cout << "t0: " << state1.print() << std::endl;
  // // state1.evolve( t, T);
  // // std::cout << "t1: " << state1.print() << std::endl;
  // // std::cout << "norm state =" << state1.squared_norm() << std::endl;
  // std::cout << "norm eps =" << v0.squared_norm() << std::endl;
  // // std::cout << v.print() << std::endl;
  // const int t=2;
  // GeneralState vt = v0;
  // vt.evolve(t, T);
  // // std::cout << "v0 = " << v0.print() << std::endl;
  // // std::cout << "vt = " << vt.print() << std::endl;
  // std::cout << v0.dot(v0) << std::endl;
  // std::cout << vt.dot(vt) << std::endl;
  // const Complex eet = v0.dot(vt) * (1.0/L/L);
  // std::cout << eet << std::endl;
  // std::cout << eet-vev*vev << std::endl;
  // // std::cout << v0.dot(vt)*(1.0/L/L) << std::endl;
  // // std::cout << v0.dot(vt)-vev*vev << std::endl;
  // // std::cout << v.print() << std::endl;















// ProdState vec=vac;
  // vec.a_tilde( 1 );
  // std::cout << vec.squared_norm() << std::endl;

  // GeneralState tmp( vac );
  // GeneralState test;
  // test.a_tilde( -3, tmp );
  // std::cout << test.squared_norm() << std::endl;
  // std::cout << test.print() << std::endl;

  // GeneralState test2;
  // test2.a_tilde( 3, test );
  // std::cout << test2.squared_norm() << std::endl;
  // std::cout << test2.print() << std::endl;










// struct SigmaState {
//   std::vector<Complex> coeffs;

//   SigmaState()
//     : coeffs(2)
//   {}

//   void set( const std::initializer_list<Complex> args ) {
//     auto itr = args.begin();
//     for( int i=0; i<2; i++ ) coeffs[i] = *(itr+i);
//   }

//   void sigma1(){ coeffs[1] *= -1.0; }

//   void sigma3(){
//     const Complex tmp = coeffs[0];
//     coeffs[0] = coeffs[1];
//     coeffs[1] = tmp;
//   }

//   void sigma2(){
//     sigma3();
//     sigma1();
//     for( int i=0; i<2; i++ ) coeffs[i] *= -I;
//   }

//   void sigmaP(){
//     const Complex tmp0 = coeffs[0];
//     const Complex tmp1 = coeffs[1];
//     coeffs[0] = 0.5*( tmp0-tmp1 );
//     coeffs[1] = 0.5*( tmp0-tmp1 );
//   }

//   void sigmaM(){
//     const Complex tmp0 = coeffs[0];
//     const Complex tmp1 = coeffs[1];
//     coeffs[0] = 0.5*( tmp0+tmp1 );
//     coeffs[1] = -0.5*( tmp0+tmp1 );
//   }

//   std::string print() const {
//     std::stringstream ss;
//     ss << std::scientific << std::setprecision(15);
//     for( Complex elem : coeffs ) ss << elem << " ";
//     return ss.str();
//   }

//   Complex dot( const SigmaState& other ) const {
//     Complex res = 0.0;
//     for(int i=0; i<2; i++) res += conj(coeffs[i])*other.coeffs[i];
//     return res;
//   }


// };

// struct ProdSigmaState {
//   std::vector<SigmaState> sigma_states;

//   ProdSigmaState()
//     : sigma_states(L)
//   {
//     set2empty();
//   }

//   void set2empty() {
//     for(int x=0; x<L; x++) sigma_states[x].set({1.0/sqrt(2), -1.0/sqrt(2)});
//   }

//   void a(const int x_){
//     const int x=(x_+L)%L;
//     for(int y=0; y<x; y++) sigma_states[y].sigma3();
//     sigma_states[x].sigmaM();
//   }

//   void a_dagger(const int x_){
//     const int x=(x_+L)%L;
//     for(int y=0; y<x; y++) sigma_states[y].sigma3();
//     sigma_states[x].sigmaP();
//   }

//   std::string print() const {
//     std::stringstream ss;
//     ss << std::scientific << std::setprecision(15);
//     for( SigmaState elem : sigma_states ) ss << elem.print() << std::endl;
//     return ss.str();
//   }

//   Complex dot( const ProdSigmaState& other ) const {
//     Complex res = 1.0;
//     for(int x=0; x<L; x++) {
//       res *= sigma_states[x].dot( other.sigma_states[x] );
//     }
//     return res;
//   }

//   Complex squared_norm() const { return dot( *this ); }

// };


// struct GenSigmaState {
//   std::vector<ProdSigmaState> states;
//   std::vector<Complex> coeffs;

//   GenSigmaState()
//     : states()
//     , coeffs()
//   {}

//   GenSigmaState( const ProdSigmaState& state )
//     : states({state})
//     , coeffs({0.0})
//   {}

//   void set_w_a_tilde( const int k_, const ProdSigmaState& state ){
//     const int k=(k_+2*L)%L;
//     const double k_phys = M_PI*k/L;

//     states.clear();
//     coeffs.clear();
//     for(int x=0; x<L; x++){
//       ProdSigmaState tmp = state;
//       tmp.a(x);

//       states.push_back(tmp);
//       coeffs.push_back( exp( I*(k_phys*x) )/sqrt(L) );
//     }
//   }

//   void set_w_a_tilde_dagger( const int k_, const ProdSigmaState& state ){
//     const int k=(k_+2*L)%L;
//     const double k_phys = M_PI*k/L;

//     states.clear();
//     coeffs.clear();
//     for(int x=0; x<L; x++){
//       ProdSigmaState tmp = state;
//       tmp.a_dagger(x);

//       states.push_back(tmp);
//       coeffs.push_back( exp( -I*(k_phys*x) )/sqrt(L) );
//     }
//   }

//   void a_tilde( const int k_, const GenSigmaState& other ){
//     const int k=(k_+2*L)%L;
//     const double k_phys = M_PI*k/L;

//     states.clear();
//     coeffs.clear();

//     for(int i=0; i<other.states.size(); i++){
//       for(int x=0; x<L; x++){
//         ProdSigmaState tmp = other.states[i];
//         tmp.a(x);

//         states.push_back(tmp);
//         coeffs.push_back( other.coeffs[i]*exp( I*(k_phys*x) )/sqrt(L) );
//       }
//     }
//   }

//   void a_tilde_dagger( const int k_, const GenSigmaState& other ){
//     const int k=(k_+2*L)%L;
//     const double k_phys = M_PI*k/L;

//     states.clear();
//     coeffs.clear();

//     for(int i=0; i<other.states.size(); i++){
//       for(int x=0; x<L; x++){
//         ProdSigmaState tmp = other.states[i];
//         tmp.a_dagger(x);

//         states.push_back(tmp);
//         coeffs.push_back( other.coeffs[i]*exp( -I*(k_phys*x) )/sqrt(L) );
//       }
//     }
//   }

//   Complex dot( const GenSigmaState& other ) const {
//     Complex res = 0.0;
//     for(int i=0; i<states.size(); i++) {
//       for(int j=0; j<other.states.size(); j++) {
//         res += conj(coeffs[i])*other.coeffs[j] * states[i].dot(other.states[j]);
//       }}
//     return res;
//   }

//   Complex squared_norm() const { return dot( *this ); }

// };







// std::cout << T.id << std::endl
//           << T.sigma1 << std::endl
//           << T.sigma2 << std::endl
//           << T.sigma3 << std::endl;
// std::cout << "beta = " << beta << std::endl
//           << "beta_tilde = " << beta_tilde << std::endl;

// const int k=3;

// std::cout << T.thetaE(k) << std::endl;
// std::cout << T.thetaE_half(k) << std::endl;
// std::cout << T.theta_tildeE(k) << std::endl;
// std::cout << T.E(k) << std::endl;
// std::cout << T.O(k) << std::endl;

// for(int k=1; k<L; k++){
//   std::cout << "TE, k=" << k << ": " << T.eval_TE[k] << std::endl;
//   M2 mat = T.E(k);

//   V2 eval = T.eval_TE[k];
//   M2 evec = T.evec_TE[k];
//   std::cout << mat * evec << std::endl;
//   std::cout << evec * eval.asDiagonal() << std::endl;

//   CV2 v = T.evec(k, true, 0);
//   double lambda = T.eval(k, true, 0);
//   std::cout << mat * v << std::endl;
//   std::cout << lambda * v << std::endl;

//   //}
// }

// ProdState state( true );
// state.set2vacuum( T );
// std::cout << state.print() << std::endl;

// state(k).set({1.0 + 0.1*I, 0.1, 0.2+0.8*I, 0.8});
// std::cout << state(k).print() << std::endl;;

// state(k).a_dagger();
// std::cout << state(k).print() << std::endl;;

// ProdSigmaState empty;
// // std::cout << sigma.print() << std::endl;
// std::cout << "norm = " << empty.squared_norm() << std::endl;

// ProdSigmaState a_empty = empty;
// a_empty.a(1);
// // std::cout << sigma.print() << std::endl;
// std::cout << "norm = " << a_empty.squared_norm() << std::endl;

// ProdSigmaState a_dagger_empty = empty;
// a_dagger_empty.a_dagger(1);
// // std::cout << sigma.print() << std::endl;
// std::cout << "norm = " << a_dagger_empty.squared_norm() << std::endl;

// GenSigmaState a_tilde_empty;
// a_tilde_empty.set_w_a_tilde( 1, empty );
// std::cout << "norm = " << a_tilde_empty.squared_norm() << std::endl;

// GenSigmaState a_tilde_dagger_empty;
// a_tilde_dagger_empty.set_w_a_tilde_dagger( 1, empty );
// std::cout << "norm = " << a_tilde_dagger_empty.squared_norm() << std::endl;
