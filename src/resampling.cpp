#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

// function to permute the ancestry vector so that A^i = i for as many i as possible
void permute(IntegerVector & a, const int & nparticles){
  int swap;
  for (int i = 0; i < nparticles; i++){
    if (a(i) != i && a(a(i)) != a(i)){
      swap = a(a(i));
      a(a(i)) = a(i);
      a(i) = swap;
      i--;
    }
  }
}

// [[Rcpp::export]]
IntegerVector systematic_resampling_(int nsamples, const NumericVector & weights){
  RNGScope scope;
  IntegerVector ancestors(nsamples);
  NumericVector u_vec = runif(1);
  double u = u_vec(0) / nsamples;
  int j = 0;
  double csw = weights(0);
  for(int k = 0; k < nsamples; k++){
    while(csw < u){
      j++;
      csw += weights(j);
    }
    ancestors(k) = j;
    u = u + 1./nsamples;
  }
  return ancestors + 1; // add 1 so that smallest is 1, and largest is N
}

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
IntegerVector multinomial_resampling_(int nsamples, const NumericVector & weights){
  RNGScope scope;
  IntegerVector ancestors(nsamples);
  int nparents = weights.size();
  NumericVector cumsumw = cumsum(weights);
  NumericVector uniforms = runif(nsamples);
  double sumw = cumsumw(nparents - 1);
  double lnMax = 0;
  int j = nparents;
  for (int i = nsamples; i > 0; i--){
    if (j > 0){
      lnMax += log(uniforms(i-1)) / i;
      uniforms(i-1) = sumw * exp(lnMax); // sorted uniforms
      while (j > 0 && uniforms(i-1) < cumsumw(j-1)){
        j --;
      }
      ancestors(i-1) = j;
    } else {
      ancestors(i-1) = 0;
    }
  }
  std::random_shuffle(ancestors.begin(), ancestors.end(), randWrapper);
  return ancestors + 1;
}
