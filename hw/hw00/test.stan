data {
  int<lower = 1> n;
  int<lower = 1> p;
  vector[n] Y;
  matrix[n, p + 1] X;
}
parameters {
  vector[p + 1] beta;
  real<lower = 0> sigma;
}
model {
  Y ~ normal(X * beta, sigma);
}
