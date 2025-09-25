#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List SimpLinCpp(NumericVector x, NumericVector y, double t_val) {
  int n = x.size();

  // Compute means
  double x_mean = mean(x);
  double y_mean = mean(y);

  // Compute square sum
  double ss_xy = 0.0;
  double ss_xx = 0.0;

  for(int i = 0; i < n; ++i) {
    double x_diff = x[i] - x_mean;
    double y_diff = y[i] - y_mean;
    ss_xy += x_diff * y_diff;
    ss_xx += x_diff * x_diff;
  }

  // Compute regression coefficients
  double beta1 = ss_xy / ss_xx;
  double beta0 = y_mean - beta1 * x_mean;



  // Compute residual
  NumericVector y_hat(n);
  NumericVector residuals(n);

  for (int i = 0; i < n; ++i) {
    y_hat[i] = beta0 + beta1 * x[i];
    residuals[i] = y[i] - y_hat[i];
  }


  // Compute residual standard errors
  double rss = 0.0;
  for (int i = 0; i < n; ++i) {
    rss += residuals[i] * residuals[i];
  }

  double rse = sqrt(rss / (n - 2)); // freedom is 2

  // Compute standard error of coefficients
  double se_beta1 = rse / sqrt(ss_xx);
  double se_beta0 = rse * sqrt((1.0 / n) + (x_mean * x_mean) / ss_xx);

  // 95% IC with t_val in R
  // Compute t_val
  //double t_val = R::qt(0.975, n-2, 1, 0); // Parameters: p, df, lower_tail, log_p
  //ci_beta0 = [beta0 - t_val * se_beta0, beta0 + t_val * se_beta0];
  //ci_beta1 = [beta1 - t_val * se_beta1, beta1 + t_val * se_beta1];
  NumericVector ci_beta0 = NumericVector::create(
    beta0 - t_val * se_beta0,
    beta0 + t_val * se_beta0
  );

  NumericVector ci_beta1 = NumericVector::create(
    beta1 - t_val * se_beta1,
    beta1 + t_val * se_beta1
  );

  // Return the result
  return List::create(
    _["coefficients"] = NumericVector::create(_["(Intercept)"] = beta0, _["x"] = beta1),
    _["std_err"] = NumericVector::create(_["(Intercept)"] = se_beta0, _["x"] = se_beta1),
    _["ci_95"] = List::create(_["(Intercept)"] = ci_beta0, _["x"] = ci_beta1),
    _["residuals"] = residuals,
    _["predicted_values"] = y_hat
  );
}
