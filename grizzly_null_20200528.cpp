#include "TMB.hpp"
#include "distributions_R_new.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  //"case" here means the observed points (i.e., case-control design)
  DATA_VECTOR( r_cases ); // observed step lengths
  DATA_VECTOR( phi_cases ); // observed turning angles
  DATA_VECTOR( delta ); // initial condition for HMM
  DATA_SCALAR( sl_stationary ); // step length parameter for stationary distribution

  // Parameters
  PARAMETER( steplength ); //exponential distribution rate parameter representing mean step length
  PARAMETER( angularconcentration ); //"kappa" parameter for von Mises distribution
  PARAMETER( logit_lambda ); //stationary state; so the true probability is always in [0,1]
  PARAMETER( logit_gamma ); //non-stationary state; so the true probability is always in [0,1]
  
  int n_cases = r_cases.size();
  Type lambda = 1 / (1 + exp(-logit_lambda));
  Type gamma = 1 / (1 + exp(-logit_gamma));
  matrix<Type> A(2, 2);
  A << lambda, 1-gamma,
       1-lambda, gamma;
  //convert our parameters to [0,1] values for the actual matrix
  
  // Objective function - "val" is our negative log likelihood for the parameters given the data
  Type val = 0.0;
  vector<Type> ncll_ns(n_cases); // negative conditional log-likelihood for non-stationary state
  vector<Type> ncll_s(n_cases); // negative conditional log-likelihood for stationary state
  
  //first we fill in the conditional log-likelihoods for each data point
  for( int i = 0; i < n_cases; i++) {
    //non-stationary state - the animal follows an exponential step length distribution and a
    ////von Mises turning angle distribution
    ncll_ns(i) = -dexp(r_cases(i), 1.0 / steplength, true);
    ncll_ns(i) -= dvonmises(phi_cases(i), Type(0), angularconcentration, true);
    //stationary state - the animal follows a half-Gaussian step length distribution and a uniform
    //von Mises turning angle distribution
    ncll_s(i) = -log(Type(2.0)) - dnorm(r_cases(i), Type(0), sl_stationary, true);
    
    //the weighting functions all turn to 0 here in the log so we don't have anything to do here
  }
  
  matrix<Type> PHI(1, 2); 
  // this will represent the relative probability of being in any state, which varies in time
  //PHI is iteratively updated based on the observed conditional likelihoods

  PHI(0,0) = delta(0);
  PHI(0,1) = delta(1); // start out at our initial probability
  
  //multiply state probabilities by continuous likelihoods (given the state)
  //exp(- (negative log likelihood)) = the likelihood
  PHI(0,0) = PHI(0,0) * exp(-ncll_s(0));
  PHI(0,1) = PHI(0,1) * exp(-ncll_ns(0));
  
  //do the first point outside the loop b/c we don't need to multiply by A
  val -= log(PHI.sum());
  PHI *= 1.0 / PHI.sum(); 
  // normalize PHI so it represents probabilities for the next iteration
  
  for (int i = 1; i < n_cases; i++) {
    PHI = PHI * A; //multiply by the transition matrix to get new state probabilities
    
    PHI(0,0) = PHI(0,0) * exp(-ncll_s(i));
    PHI(0,1) = PHI(0,1) * exp(-ncll_ns(i));
    
    val -= log(PHI.sum());
    
    PHI *= 1.0 / PHI.sum(); //normalize PHI for the next step
  }
  
  // Reporting
  REPORT(val);
  REPORT( steplength );
  REPORT( angularconcentration );
  REPORT( logit_lambda );
  REPORT( logit_gamma );
  REPORT( r_cases ); // observed step lengths
  REPORT( phi_cases ); // observed turning angles
  REPORT( delta ); // initial condition for HMM
  REPORT( sl_stationary ); // we fix it for now
  ADREPORT( steplength );
  ADREPORT( angularconcentration );
  ADREPORT( logit_lambda );
  ADREPORT( logit_gamma );
  
  return val;
}
