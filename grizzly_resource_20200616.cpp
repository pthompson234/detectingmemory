#include "TMB.hpp"
#include "distributions_R_new.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_simmed ); //number of simulated steps per actual step
  DATA_VECTOR( r_cases ); // observed step lengths
  DATA_VECTOR( phi_cases ); // observed turning angles
  DATA_ARRAY( res_cases ); //environmental covariate values that were observed
  DATA_ARRAY( res_controls ); //environmental covariate values that were simulated
  DATA_VECTOR( delta ); // initial condition for HMM
  DATA_SCALAR( sl_stationary ); // we fix it for now
  
  // Parameters
  PARAMETER( steplength );
  PARAMETER( angularconcentration );
  PARAMETER_VECTOR( beta_res );
  PARAMETER( logit_lambda ); //stationary; so the true probability is always in [0,1]
  PARAMETER( logit_gamma ); //non-stationary; so the true probability is always in [0,1]
  
  int n_cases = r_cases.size();
  int k = res_cases.cols();
  Type lambda = 1 / (1 + exp(-logit_lambda));
  Type gamma = 1 / (1 + exp(-logit_gamma));
  matrix<Type> A(2, 2);
  A << lambda, 1-gamma,
       1-lambda, gamma;
  
  // Objective function
  Type val = 0.0;
  vector<Type> ncll_ns(n_cases); // negative conditional log-likelihood for non-stationary state
  vector<Type> ncll_tilde(n_cases); // for calculation of PHI
  vector<Type> ncll_s(n_cases); // negative conditional log-likelihood for stationary state
  
  for( int i = 0; i < n_cases; i++) {
    ncll_ns(i) = -dexp(r_cases(i), 1.0 / steplength, true);
    ncll_ns(i) -= dvonmises(phi_cases(i), Type(0), angularconcentration, true);
    ncll_s(i) = -log(Type(2.0)) - dnorm(r_cases(i), Type(0), sl_stationary, true);
    ncll_tilde(i) = ncll_ns(i);
    
    Type cs = 0;
    for (int p = 0; p < k; p++) {
      cs += beta_res(p) * res_cases(i, p);
    }
    
    vector<Type> arr(n_simmed);
    for (int j = 0; j < n_simmed; j++ ) {
      arr(j) = 0;
      for (int p = 0; p < k; p++) {
        arr(j) += beta_res(p) * res_controls(i * n_simmed + j, p);
      }
    }
    
    Type sum = logspace_add(arr(0), arr(1));
    for (int j = 2; j < n_simmed; j++) {
      sum = logspace_add(sum, arr(j));
    }
    
    ncll_ns(i) -= (cs - (sum - log(n_simmed)));
  }
  
  matrix<Type> PHI_tilde(1, 2); 
  // this will represent the relative probability of being in any state, which varies in time
  
  PHI_tilde(0,0) = delta(0);
  PHI_tilde(0,1) = delta(1); 
  // start out at our initial probability, which is just a guess
  
  //do the first point here b/c we don't need to multiply by A
  val -= log(PHI_tilde(0,0) * exp(-ncll_s(0)) + PHI_tilde(0,1) * exp(-ncll_ns(0)));
  
  PHI_tilde(0,0) = PHI_tilde(0,0) * exp(-ncll_s(0));
  PHI_tilde(0,1) = PHI_tilde(0,1) * exp(-ncll_tilde(0));
  PHI_tilde *= 1.0 / PHI_tilde.sum(); 
  // reset PHI so it represents probabilities for the next iteration
  
  for (int i = 1; i < n_cases; i++) {
    PHI_tilde = PHI_tilde * A; //multiply by the transition matrix to get new state probabilities
    
    val -= log(PHI_tilde(0,0) * exp(-ncll_s(i)) + PHI_tilde(0,1) * exp(-ncll_ns(i)));
    
    PHI_tilde(0,0) = PHI_tilde(0,0) * exp(-ncll_s(i));
    PHI_tilde(0,1) = PHI_tilde(0,1) * exp(-ncll_tilde(i));
    PHI_tilde *= 1.0 / PHI_tilde.sum();
  }
  
  // Reporting
  REPORT(val);
  REPORT( steplength );
  REPORT( angularconcentration );
  REPORT( beta_res );
  REPORT( lambda );
  REPORT( gamma );
  
  ADREPORT( steplength );
  ADREPORT( angularconcentration );
  ADREPORT( beta_res );
  ADREPORT( lambda );
  ADREPORT( gamma );
  
  return val;
}
