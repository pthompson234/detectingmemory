#include "TMB.hpp"
#include "distributions_R_new.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_simmed ); //number of simulated steps per actual step
  DATA_VECTOR( r_cases ); // observed step lengths. For data being fit only!
  DATA_VECTOR( phi_cases ); // observed turning angles. For data being fit only!
  DATA_VECTOR( x_cases ); // x-coordinate of observed steps. For all data!
  DATA_VECTOR( y_cases ); // y-coordinate of observed steps. For all data!
  DATA_VECTOR( x_controls ); // x-coordinate of simulated steps. For data being fit only!
  DATA_VECTOR( y_controls ); // y-coordinate of simulated steps. For data being fit only!
  DATA_FACTOR( strata_controls ); // for controls, the time index of the simulated points. For data being fit only!
  DATA_VECTOR( delta ); // initial condition for HMM
  DATA_SCALAR( sl_stationary ); // we fix it for now
  
  // Parameters
  PARAMETER( steplength );
  PARAMETER( angularconcentration );
  PARAMETER( dist_coef ); // length n_peaks - weight for each component of the lagged distance
  PARAMETER( mean_tau ); // mean of memory-based movement cycle; length n_peaks
  PARAMETER( sd_tau ); // fix for now; variation in memory-based movement cycle
  PARAMETER( logit_lambda ); //stationary; so the true probability is always in [0,1]
  PARAMETER( logit_gamma ); //non-stationary; so the true probability is always in [0,1]
  PARAMETER( log_alpha ); // we fix it for now
  
  int n_cases = r_cases.size();
  int n_controls = x_controls.size();
  
  Type exp_alpha = exp(log_alpha);
  Type lambda = 1 / (1 + exp(-logit_lambda));
  Type gamma = 1 / (1 + exp(-logit_gamma));
  matrix<Type> A(2, 2);
  A << lambda, 1-gamma,
       1-lambda, gamma;
  
  // Objective function
  Type val = 0.0;
  
  //get "distance from point at t - tau" values for every case and control
  vector<Type> dist_cases(n_cases);
  vector<Type> dist_controls(n_controls);
  
  int last_val = strata_controls(n_controls - 1) + 1;
  vector<Type> dnorms(last_val);
  for (int i = 0; i < last_val; i++) {
    dnorms(i) = dnorm(Type(i), mean_tau, sd_tau);
  }
  
  for (int i = 0; i < n_cases; i++) {
    int t = strata_controls(i * n_simmed);
    
    vector<Type> dnorms_sub = dnorms.segment(1, t);
    Type sumNorm = dnorms_sub.sum();
    vector<Type> dist_vars_cs(t);
    for (int tau = 1; tau <= t; tau++) {
      // tau ranges from 1 to the "maximum" value for tau, i.e., the current value of t for this point in space
      // tau = 0 makes no sense!
      dist_vars_cs(tau - 1) = sqrt(pow(x_cases(t) - x_cases(t - tau), 2) + pow(y_cases(t) - y_cases(t - tau), 2));
    }
    dist_cases(i) = (dnorms_sub * exp(-exp_alpha * dist_vars_cs)).sum() * (1.0 / sumNorm);
    
    for (int j = 0; j < n_simmed; j++) {
      vector<Type> dist_vars_ct(t);
      for (int tau = 1; tau <= t; tau++) {
        dist_vars_ct(tau - 1) = sqrt(pow(x_controls(i * n_simmed + j) - x_cases(t - tau), 2) + pow(y_controls(i * n_simmed + j) - y_cases(t - tau), 2));
      }
      dist_controls(i * n_simmed + j) = (dnorms_sub * exp(-exp_alpha * dist_vars_ct)).sum() * (1.0 / sumNorm);
    }
  }
  
  vector<Type> ncll_tilde = -dexp(r_cases, 1.0 / steplength, true) - dvonmises(phi_cases, Type(0.0), angularconcentration, true);
  // negative conditional log-likelihood for non-stationary state
  vector<Type> ncll_s = -log(Type(2.0)) - dnorm(r_cases, Type(0.0), sl_stationary, true); // negative conditional log-likelihood for stationary state
  
  vector<Type> cs_vector = dist_coef * dist_cases;
  vector<Type> sums(n_cases);
  
  for( int i = 0; i < n_cases; i++) {
    vector<Type> arr = dist_coef * dist_controls.segment(i * n_simmed, n_simmed);
    
    Type sm = logspace_add(arr(0), arr(1));
    for (int j = 2; j < n_simmed; j++) {
      sm = logspace_add(sm, arr(j));
    }
    
    sums(i) = sm;
  }
  
  vector<Type> ncll_ns = ncll_tilde - (cs_vector - (sums - log(n_simmed)));
  
  matrix<Type> PHI_tilde(1, 2); 
  // this will represent the relative probability of being in any state, which varies in time
  
  PHI_tilde(0,0) = delta(0);
  PHI_tilde(0,1) = delta(1); 
  // start out at our initial probability, which is just a guess
  
  //do the first point here b/c we don't need to multiply by A
  val -= logspace_add(log(PHI_tilde(0,0)) - ncll_s(0), log(PHI_tilde(0,1)) - ncll_ns(0));
  
  PHI_tilde(0,0) = PHI_tilde(0,0) * exp(-ncll_s(0));
  PHI_tilde(0,1) = PHI_tilde(0,1) * exp(-ncll_tilde(0));
  PHI_tilde *= 1.0 / PHI_tilde.sum(); 
  // reset PHI so it represents probabilities for the next iteration
  
  for (int i = 1; i < n_cases; i++) {
    PHI_tilde = PHI_tilde * A; //multiply by the transition matrix to get new state probabilities
    
    val -= logspace_add(log(PHI_tilde(0,0)) - ncll_s(i), log(PHI_tilde(0,1)) - ncll_ns(i));
    
    PHI_tilde(0,0) = PHI_tilde(0,0) * exp(-ncll_s(i));
    PHI_tilde(0,1) = PHI_tilde(0,1) * exp(-ncll_tilde(i));
    PHI_tilde *= 1.0 / PHI_tilde.sum();
  }
  
  // Reporting
  REPORT( val );
  REPORT( steplength );
  REPORT( angularconcentration );
  REPORT( dist_coef );
  REPORT( mean_tau );
  REPORT( sd_tau );
  REPORT( lambda );
  REPORT( gamma );
  REPORT( log_alpha );
  
  REPORT(ncll_ns);
  REPORT(ncll_s);
  REPORT(ncll_tilde);
  REPORT(dist_cases);
  REPORT(dist_controls);
  REPORT(cs_vector);
  REPORT(sums);
  
  ADREPORT( steplength );
  ADREPORT( angularconcentration );
  ADREPORT( dist_coef );
  ADREPORT( mean_tau );
  ADREPORT( sd_tau );
  ADREPORT( lambda );
  ADREPORT( gamma );
  ADREPORT( log_alpha );
  
  return val;
}
