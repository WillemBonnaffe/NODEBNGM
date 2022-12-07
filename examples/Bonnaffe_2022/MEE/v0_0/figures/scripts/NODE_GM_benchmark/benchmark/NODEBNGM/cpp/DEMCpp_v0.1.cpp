//* MCMCMH ALGORITHM *//

// goal: code an MCMC algorithm implemented with MH

// updates:
// @ v0.1 - 17/11/2019 - implemented automatic cacluation of gamma


/* ********** */
/* INITIATION */
/* ********** */

// C++ dependencies
#include <Rcpp.h>

using namespace Rcpp;


/* ********** */
/* ITERATIONS */
/* ********** */

// [[Rcpp::export]]
List DEMCpp(List argList)
{

  // unwrap arguments
  Function dTarget = as<Function>(argList["dTarget"]);
  NumericVector Theta_local = as<NumericVector>(argList["Theta_0"]); int d = Theta_local.size();
  float gamma = 2.38/pow(2.0*d,0.5); // float gamma = as<float>(argList["gamma"]);
  float epsilon = as<float>(argList["epsilon"]);
  int nIt = as<int>(argList["nIt"]);

  // initiate chain
  NumericMatrix chain(nIt, d + 1);
  float target_local = as<float>(dTarget(Theta_local));
  NumericVector chain_local; chain_local = Theta_local; chain_local.push_front(target_local);
  chain(0,_) = chain_local;

  // iterations
  float r; int accepted = 0; int count = 0;
  for(int k = 0 ; k < nIt ; k++)
  {

    // set distance vector with lag
    // int lag = 100; if(lag > k){lag = k;}
    NumericVector delta = chain(as<int>(runif(1,k/2,k)),_) - chain(as<int>(runif(1,k/2,k)),_); delta.erase(0);
    // NumericVector delta = chain(as<int>(runif(1,0,k)),_) - chain(as<int>(runif(1,0,k)),_); delta.erase(0); // set distance vector (no lag)

    // update theta
    NumericVector Theta_newLocal = Theta_local + gamma*delta + epsilon*as<NumericVector>(runif(d,-1.0,1.0)); // uniform walker
    // Theta_newLocal = Theta_local + gamma*delta + epsilon*as<NumericVector>(rnorm(d,0,1)); // normal walker
    // (optional) Gelman's trick: uniform jumper as<float>(runif(1,0,2.0))*gamma*delta

    // update target
    float target_newLocal = as<float>(dTarget(Theta_newLocal));

    // metropolis ratio
    r = exp(target_newLocal - target_local);

    // metropolis test
    if(drand48() < r)
    {

      accepted += 1;
      target_local = target_newLocal;
      Theta_local = Theta_newLocal;

    }

    // update the chain
    chain_local = Theta_local; chain_local.push_front(target_local);
    chain(k,_) = chain_local;

    // message
    if(count == (int)(nIt/10))
    {
      Rcout << k << "/" << nIt << " | " << chain_local[0] << " | " << chain_local[1] << " | " << chain_local[2] << "\n";
      count = 0;
    }
    count +=1;


  }


  // return results
  Rcout << "p = " << (float)accepted/(float)nIt << "\n";
  return Rcpp::List::create(
    Rcpp::Named("dTarget") = dTarget(Theta_local),
    Rcpp::Named("Theta_0") = Theta_local,
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("chainList") = chain
    );

}
