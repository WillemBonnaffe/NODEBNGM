//* MCMCMH ALGORITHM *//

// goal: code an MCMC algorithm implemented with MH

// updates:
// @ v0.3 - 20/06/2019 - implemented storage of all steps again
// @ v0.4 - 21/06/2019 - removed the order of the samples taken
// @ v0.5 - 17/09/2019 - added automatic calculation of gamma


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
List DEMCOpp(List argList)
{

  // unwrap arguments
  Function dTarget = as<Function>(argList["dTarget"]);
  NumericVector Theta_local = as<NumericVector>(argList["Theta_0"]); int d = Theta_local.size();
  float gamma = 2.38/pow(2.0*d,0.5);
  float epsilon = as<float>(argList["epsilon"]);
  int lambda = as<int>(argList["lambda"]);
  int nIt = as<int>(argList["nIt"]);

  // initiate chain
  NumericMatrix chain(nIt, d + 1);
  float target_local = as<float>(dTarget(Theta_local));
  NumericVector chain_local; chain_local = Theta_local; chain_local.push_front(target_local);
  chain(0,_) = chain_local;

  // iterations
  int accepted = 0; int count = 0;
  for(int k = 0 ; k < nIt ; k++)
  {

    // set distance vector with lag
    int lag = lambda; if(lag > k){lag = k;}
    NumericVector delta = chain(as<int>(runif(1,k-lag,k)),_) - chain(as<int>(runif(1,k-lag,k)),_); delta.erase(0);
    // NumericVector delta = chain(as<int>(runif(1,0,k)),_) - chain(as<int>(runif(1,0,k)),_); delta.erase(0); // set distance vector (no lag)

    // update theta
    NumericVector Theta_newLocal = Theta_local + as<float>(runif(1,0,2.0))*gamma*delta + epsilon*as<NumericVector>(runif(d,-1.0,1.0)); // uniform walker
    // Theta_newLocal = Theta_local + gamma*delta + epsilon*as<NumericVector>(rnorm(d,0,1)); // normal walker
    // (optional) Gelman's trick: uniform jumper as<float>(runif(1,0,2.0))*gamma*delta

    // update target
    float target_newLocal = as<float>(dTarget(Theta_newLocal));

    // test
    if(target_newLocal > target_local)
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
    Rcpp::Named("chainList") = chain,
    Rcpp::Named("p") = (float)accepted/(float)nIt
    );

}
