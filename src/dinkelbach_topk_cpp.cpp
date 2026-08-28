#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List dinkelbach_topk_cpp(NumericVector x, NumericVector r, int k, int sgn,
                         double sum_x2, int max_iter = 50, double tol = 1e-9) {
  const int n = x.size();
  if (r.size() != n) stop("x and r must have same length");
  if (k < 1 || k > n) stop("invalid k");
  if (!(sgn == 1 || sgn == -1)) stop("sgn must be +/-1");
  
  std::vector<double> nval(n), dval(n), w(n);
  std::vector<int> ids(n);
  for (int i = 0; i < n; ++i) {
    nval[i] = sgn * x[i] * r[i];
    dval[i] = -x[i] * x[i];
    ids[i] = i;
  }
  
  double lambda = 0.0;
  int iterations = 0;
  std::vector<int> top(k);
  
  for (int it = 0; it < max_iter; ++it) {
    iterations = it + 1;
    for (int i = 0; i < n; ++i) w[i] = nval[i] - lambda * dval[i];
    
    std::iota(ids.begin(), ids.end(), 0);
    std::partial_sort(ids.begin(), ids.begin() + k, ids.end(),
                      [&](int a, int b) { return w[a] > w[b]; });
    
    double num = 0.0, den = sum_x2;
    for (int j = 0; j < k; ++j) {
      top[j] = ids[j];
      num += nval[top[j]];
      den += dval[top[j]];
    }
    if (std::abs(den) < 1e-15) break;
    const double new_lambda = num / den;
    if (std::abs(new_lambda - lambda) < tol) {
      lambda = new_lambda;
      break;
    }
    lambda = new_lambda;
  }
  
  IntegerVector out(k);
  for (int j = 0; j < k; ++j) out[j] = top[j] + 1; // R is 1-based
  return List::create(
    _["indices"] = out,
    _["dfbeta"] = sgn * lambda,
    _["lambda"] = lambda,
    _["iterations"] = iterations
  );
}