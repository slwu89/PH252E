#include <Rcpp.h>
#include <math.h> 
using namespace Rcpp;

/*
 * PH252E Causal Inference 2: Lab 1
 * Sean Wu 9/14/2016
 */


/*
 * Question 2
 * Given vector x of length n, write a function that returns n iid random variables 
 * uniformly distrubuted between 0 and x.
 */
// [[Rcpp::export]]
NumericVector uniform_cpp(NumericVector x){
  NumericVector out(x.size());
  for(int i=0; i<x.size(); i++){
    out[i] = R::runif(0.0,x[i]);
  }
  return(out);
}

/***R
set.seed(1)
uniform_cpp(c(10,100,1000))

#benchmark
library(microbenchmark)
uniform <- function(x){
  out <- sapply(x,function(xx){runif(n=1,min=0,max=xx)})
  return(out)
}
microbenchmark(
  uniform(c(1e1,1e2,1e3)),
  uniform_cpp(c(1e1,1e2,1e3)),
  times=1e3)
*/


/*
 * Question 3
 * Write a function called rexpit that takes input vector x and returns a vector of iid
 * Bernoulli random variables z, where P(z = 1) = expit(b).
 */
// [[Rcpp::export]]
double plogis_cpp(double x){
  double out;
  out = 1 / (1 + exp(-x));
  return(out);
}

// [[Rcpp::export]]
NumericVector rexpit_cpp(NumericVector x){
  NumericVector out(x.size());
  for(int i=0; i<x.size(); i++){
    out[i] = R::rbinom(1.0,plogis_cpp(x[i]));
  }
  return(out);
}

/*
 * Find mean(rexpit(rnorm(1e6, mean = 2, sd = 1)))
 */
/***R
set.seed(1)
mean(rexpit_cpp(rnorm(1e6, mean = 2, sd = 1)))

#benchmark
set.seed(1)
rexpit <- function(x) {
  return(rbinom(length(x), size = 1, prob = plogis(x)))
}
mean(rexpit(rnorm(1e6, mean = 2, sd = 1)))
# microbenchmark(
#   rexpit(rnorm(1e6, mean = 2, sd = 1)),
#   rexpit_cpp(rnorm(1e6, mean = 2, sd = 1)),
#   times=20)
*/


/*
 * Question 4
 * Now suppose instead that we want iid standard normal random variables, a and b and Bernoulli 
 * random variable z where (equation). Find E[z|a > −1].
 */
// [[Rcpp::export]]
double question4_cpp(int n){
  double out;
  NumericVector a = rnorm(n,0.0,1.0);
  NumericVector b = rnorm(n,0.0,1.0);
  NumericVector z(n,NA_REAL);
  //logical indexing with Rcpp sugar (won't work for STL objects)
  z[a > 0.7] = rexpit_cpp(b[a > 0.7]); //scalar version of rexpit_cpp
  z[a > -1.0 & a <= 0.7] = 1.0;
  NumericVector out_vector = z[a > -1.0];
  out = mean(out_vector);
  return(out);
}

/***R
set.seed(2)
n <- 1e6
a <- rnorm(n)
b <- rnorm(n)
z2 <- rep(NA, n)
z2[a > 0.7] <- rexpit(b[a > 0.7]) #logical indexing 
z2[a > -1 & a <= 0.7] <- 1
answer2 <- mean(z2[a > -1])
answer2

set.seed(2)
n <- 1e6
question4_cpp(n)
*/


/*
 * Question 5
 * Consider the data generating process O = (L(0), A(0), L(1), A(1), ..., L(29), A(29), L(30))
 * Find E[L(30)]
 * Hint: Use a loop!
 */
// [[Rcpp::export]]
double question5_cpp(int n, int t_max){
  NumericVector prevA(n,0.0);
  NumericVector prevL(n,0.0);
  NumericVector L;
  NumericVector A;
  for(int i=0; i<t_max; i++){
    //none of this will work in standard C++; relies on Rcpp syntatic sugar
    L = prevA + prevL + rnorm(n,0.0,1.0);
    A = rexpit_cpp(-1.5 * L + prevA);
    prevL = L;
    prevA = A;
  }
  double out;
  out = mean(L);
  return(out);
}

/***R
set.seed(1)
n <- 1e6
prev.A <- 0 # you can also do this without using prev.A and prev.L but I think
prev.L <- 0 # this is easier to read, especially when equations become more complicated 
for (t in 1:30) {
  L <- prev.A + prev.L + rnorm(n)
  A <- rexpit(-1.5 * L + prev.A)
  prev.L <- L
  prev.A <- A 
}
mean(L)

set.seed(1)
question5_cpp(n=1e6,t_max=30)
*/


/*
 * 3: Simulating the data
 * Data Structure and running example for this course
 */