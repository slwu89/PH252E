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
n <- 1e3
a <- rnorm(n)
b <- rnorm(n)
z2 <- rep(NA, n)
z2[a > 0.7] <- rexpit(b[a > 0.7]) #logical indexing 
z2[a > -1 & a <= 0.7] <- 1
answer2 <- mean(z2[a > -1])
answer2

set.seed(2)
n <- 1e3
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
n <- 1e3
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
question5_cpp(n=1e3,t_max=30)
*/


/*
 * 3: Simulating the data
 * Data Structure and running example for this course.
 * In reality we might observe n iid copies of O from the observed data distribution, P0 ∈ M, a statistical model in which the true model, P0, lives. 
 * In this lab we will create the true underlying distribution out of which P0 arises. We will therefore be able to use our knowledge of P0 
 * (i.e. the true underlying data generating functions) to create a very large hypothetical target population. 
 * This will allow us to evaluate the true value of parameters of both P0 and causal parameters under specific interventions on P0.
 */


/*
 * 3.1: Generating the observed data
 * We will now look at specific structural equations for the longitudinal data generating process described above.
 */
// [[Rcpp::export]]
NumericMatrix GenerateData_cpp(int n){

  NumericVector prevD(n,0.0);
  NumericVector prevA(n,0.0);
  NumericVector prevY(n,0.0);
  NumericVector Y0 = prevY;

  NumericMatrix out(n,13);
  out(_,0) = Y0;

  int i = 1; //define output iterator
  
  for(int t=0; t<4; t++){

    LogicalVector alive;
    alive = prevY == 0.0;
    
    //Rcpp's API doesn't allow subsetting a vector using a logical if said vector has NAs; find better workaround later
    NumericVector D(n,9999.0);
    NumericVector A(n,9999.0);
    NumericVector Y(n,9999.0);

    NumericVector u = rnorm(n,0.0,1.0);

    //updateD
    NumericVector prevA_alive = prevA[alive];
    NumericVector prevD_alive = prevD[alive];
    NumericVector u_alive = u[alive];
    NumericVector D_alive_update = 0.3 * prevA_alive + 0.9 * prevD_alive + 0.2 * u_alive;
    D[alive] = D_alive_update;

    //update A
    NumericVector D_alive_prevA = D[alive & prevA == 0];
    A[alive & prevA == 0.0] = rexpit_cpp(0.7 - 0.6 * D_alive_prevA);
    A[alive & prevA == 1] = 1.0;

    //update Y
    NumericVector A_alive = A[alive];
    NumericVector D_alive = D[alive];
    Y[alive] = rexpit_cpp(-0.2 - 0.5 * A_alive - 0.7 * D_alive);
    Y[!alive] = 1.0;

    //update previous values
    prevD = D;
    prevA = A;
    prevY = Y;
  
    out(_,i) = D;
    out(_,i+1) = A;
    out(_,i+2) = Y;
    
    i = i + 3;
  }
  
  //replace temporary null value with NA_REAL
  LogicalVector na_index = out == 9999.0;
  for(int i=0; i<na_index.size(); i++){
    if(na_index[i]==true){
      out[i] = NA_REAL;
    } else {
      continue;
    }
  }
  
  return(out);
}


/***R

GenerateData <- function(n) {
  prev.D <- prev.A <- prev.Y <- rep(0, n) 
  cum.df <- data.frame(Y0 = prev.Y) #bonus 
  for (t in 1:4) {
    alive <- prev.Y == 0
    D <- A <- Y <- rep(NA, n)
    u <- rnorm(n)
    D[alive] <- 0.3 * prev.A[alive] + 0.9 * prev.D[alive] + 0.2 * u[alive]
    A[alive & prev.A == 0] <- rexpit(0.7 - 0.6 * D[alive & prev.A == 0])
    A[alive & prev.A == 1] <- 1
    Y[alive] <- rexpit(-0.2 - 0.5 * A[alive] - 0.7 * D[alive])
    Y[!alive] <- 1
    prev.D <- D
    prev.A <- A
    prev.Y <- Y
    df <- data.frame(D, A, Y) #bonus
    names(df) <- paste0(c("D", "A", "Y"), c(t - 1, t - 1, t)) #bonus 
    cum.df <- data.frame(cum.df, df) #bonus
  }
  print(head(cum.df, 10), digits = 3) #bonus
  return(mean(Y)) 
}

set.seed(1)
q5_r <- GenerateData(1e2)

set.seed(1)
q5_cpp <- GenerateData_cpp(1e2)
q5_cpp <- mean(q5_cpp[,13])
*/


/*
 * 4: Generating the True Values of Target Causal Parameters
 * 4.1: Generating counterfactual data under a given regime.
 * Now let’s intervene on treatment and set the regime to a ̄ = (1, 1, 1, 1) by enforcing a switch time at time 0.
 * This will generate counterfactual observed data: L ̄(1,1,1,1)(t), i.e. the time-varying values of the non- intervention covariates (CD4 count and survival)
 * under an intervention to swtich treatment immediately, a ̄ = (1,1,1,1).
 */
// [[Rcpp::export]]
NumericMatrix GenerateData_abar_cpp(int n, NumericVector abar){
  
  NumericVector prevD(n,0.0);
  NumericVector prevA(n,0.0);
  NumericVector prevY(n,0.0);
  NumericVector Y0 = prevY;
  
  NumericMatrix out(n,13);
  out(_,0) = Y0;
  
  int i = 1; //define output iterator
  
  for(int t=0; t<4; t++){
    
    LogicalVector alive;
    alive = prevY == 0.0;
    
    //Rcpp's API doesn't allow subsetting a vector using a logical if said vector has NAs; find better workaround later
    NumericVector D(n,9999.0);
    NumericVector A(n,9999.0);
    NumericVector Y(n,9999.0);
    
    NumericVector u = rnorm(n,0.0,1.0);
    
    //updateD
    NumericVector prevA_alive = prevA[alive];
    NumericVector prevD_alive = prevD[alive];
    NumericVector u_alive = u[alive];
    NumericVector D_alive_update = 0.3 * prevA_alive + 0.9 * prevD_alive + 0.2 * u_alive;
    D[alive] = D_alive_update;
    
    //update A
    A[alive] = abar[t];
    
    //update Y
    NumericVector A_alive = A[alive];
    NumericVector D_alive = D[alive];
    Y[alive] = rexpit_cpp(-0.2 - 0.5 * A_alive - 0.7 * D_alive);
    Y[!alive] = 1.0;
    
    //update previous values
    prevD = D;
    prevA = A;
    prevY = Y;
    
    out(_,i) = D;
    out(_,i+1) = A;
    out(_,i+2) = Y;
    
    i = i + 3;
  }
  
  //replace temporary null value with NA_REAL
  LogicalVector na_index = out == 9999.0;
  for(int i=0; i<na_index.size(); i++){
    if(na_index[i]==true){
      out[i] = NA_REAL;
    } else {
      continue;
    }
  }
  
  return(out);
}

/***R
GenerateData <- function(n, abar) {
  prev.D <- prev.A <- prev.Y <- rep(0, n) 
  cum.df <- data.frame(Y0 = prev.Y) #bonus 
  for (t in 1:4) {
  alive <- prev.Y == 0
  D <- A <- Y <- rep(NA, n)
  u <- rnorm(n)
  D[alive] <- 0.3 * prev.A[alive] + 0.9 * prev.D[alive] + 0.2 * u[alive]
  A[alive] <- abar[t]
  Y[alive] <- rexpit(-0.2 - 0.5 * A[alive] - 0.7 * D[alive])
  Y[!alive] <- 1
  prev.D <- D
  prev.A <- A
    prev.Y <- Y
    df <- data.frame(D, A, Y) #bonus
    names(df) <- paste0(c("D", "A", "Y"), c(t - 1, t - 1, t)) #bonus 
    cum.df <- data.frame(cum.df, df) #bonus
  }
  print(head(cum.df, 10), digits = 3) #bonus 
  return(mean(Y))
}

set.seed(1)
EY_1111 <- GenerateData(n = 1e6, abar = c(1, 1, 1, 1))

set.seed(1)
EY_1111_cpp <- GenerateData_abar_cpp(1e6,c(1,1,1,1))
EY_1111_cpp <- mean(EY_1111_cpp[,13])
*/


/*
 * 4.2: Target Causal Parameters
 * Repeat what we just did, except set A(t) = 0 for all t and compute the parameters defined below.
 * This will allow us to generate counterfactual data under an intervention to never switch treatments.
 */

/***R
set.seed(1)
EY_0000 <- GenerateData(1e6,c(0,0,0,0))

set.seed(1)
EY_0000_cpp <- GenerateData_abar_cpp(1e6,c(0,0,0,0))
EY_0000_cpp <- mean(EY_0000_cpp[,13])

risk.difference.r <- EY_1111 - EY_0000
risk.difference.cpp <- EY_1111_cpp - EY_0000_cpp

risk.difference.r
risk.difference.cpp

risk.ratio.r <- EY_1111/EY_0000
risk.ratio.cpp <- EY_1111_cpp/EY_0000_cpp

risk.ratio.r
risk.ratio.cpp
*/


/*
 * 5: Computing true target parameter of a marginal structural model (MSM)
 * We suspect the quicker one intervenes to swtich to second line treatment, the lower the risk of mortality. 
 * We will call the time of intervention, the switch time, which can take values from 0 to 4, where 4 would be never switching.
 */

/*
 * 5.1: Generating counterfactual means for all switch times, s, and times, t
 * First, we notice we have 5 regimes, each corresponding to a switch time, s = 0, 1, 2, 3, 4.
 * For each regime, we also have a counterfactual mean for Y (t) at each time, t = 1, 2, 3, 4.
 * Find E[Y (t)a ̄] for each a ̄ corresponding to s = 0, 1, 2, 3, 4, t = 1, 2, 3, 4
 */
// [[Rcpp::export]]
double GenerateData_abar_tmax_cpp(int n, int t_max, NumericVector abar){
  
  NumericVector prevD(n,0.0);
  NumericVector prevA(n,0.0);
  NumericVector prevY(n,0.0);
  NumericVector Y0 = prevY;
  
  NumericMatrix out(n,(3*t_max)+1);
  out(_,0) = Y0;
  
  int i = 1; //define output iterator
  
  for(int t=0; t<t_max; t++){
    
    LogicalVector alive;
    alive = prevY == 0.0;
    
    //Rcpp's API doesn't allow subsetting a vector using a logical if said vector has NAs; find better workaround later
    NumericVector D(n,9999.0);
    NumericVector A(n,9999.0);
    NumericVector Y(n,9999.0);
    
    NumericVector u = rnorm(n,0.0,1.0);
    
    //updateD
    NumericVector prevA_alive = prevA[alive];
    NumericVector prevD_alive = prevD[alive];
    NumericVector u_alive = u[alive];
    NumericVector D_alive_update = 0.3 * prevA_alive + 0.9 * prevD_alive + 0.2 * u_alive;
    D[alive] = D_alive_update;
    
    //update A
    A[alive] = abar[t];
    
    //update Y
    NumericVector A_alive = A[alive];
    NumericVector D_alive = D[alive];
    Y[alive] = rexpit_cpp(-0.2 - 0.5 * A_alive - 0.7 * D_alive);
    Y[!alive] = 1.0;
    
    //update previous values
    prevD = D;
    prevA = A;
    prevY = Y;
    
    out(_,i) = D;
    out(_,i+1) = A;
    out(_,i+2) = Y;
    
    i = i + 3;
  }
  
  //replace temporary null value with NA_REAL
  LogicalVector na_index = out == 9999.0;
  for(int i=0; i<na_index.size(); i++){
    if(na_index[i]==true){
      out[i] = NA_REAL;
    } else {
      continue;
    }
  }
  
  return(mean(out(_,out.ncol()-1)));
}

// [[Rcpp::export]]
NumericVector make_abar(int s){
  NumericVector abar(s+(4-s));
  for(int i=0; i<abar.length(); i++){
    if(i < s){
      abar[i] = 0.0;
    } else {
      abar[i] = 1.0;
    }
  }
  return(abar);
}

// sweep through values of s and t
// [[Rcpp::export]]
DataFrame st_sweep(){
  NumericVector s_vec;
  NumericVector t_vec;
  NumericVector EYa_vec;
  for(int s=0; s<5; s++){
    for(int t=1; t<5; t++){
      Rcout << "on iteration s: " << s << ", t: " << t << std::endl;
      //make abar
      NumericVector abar = make_abar(s);
      double EYa = GenerateData_abar_tmax_cpp(1E6,t,abar);
      s_vec.push_back(s);
      t_vec.push_back(t);
      EYa_vec.push_back(EYa);
    }
  }
  
  return(DataFrame::create(_["t"]=t_vec,_["s"]=s_vec,_["EYa"]=EYa_vec));
}



/***R
GenerateData <- function(n, max.time, abar) { 
  prev.D <- rep(0, n)
  prev.A <- rep(0, n)
  prev.Y <- rep(0, n)
  for (t in 1:max.time) {
    alive <- prev.Y == 0
    D <- A <- Y <- rep(NA, n)
    u <- rnorm(n)
    D[alive] <- 0.3 * prev.A[alive] + 0.9 * prev.D[alive] + .2*u[alive] 
    A[alive] <- abar[t]
    Y[alive] <- rexpit(-0.2 - 0.5 * A[alive] - 0.7 * D[alive]) 
    Y[!alive] <- 1
    prev.D <- D
    prev.A <- A
    prev.Y <- Y
  }
  return(mean(Y)) 
}

set.seed(1)
df <- NULL
for (s in 0:4) {
  for (t in 1:4) {
    print(paste0("on iteration s: ",s,", t: ",t))
    abar <- c(rep(0, s), rep(1, 4 - s)) 
    EYa <- GenerateData(n = 1e6, t, abar) 
    df <- rbind(df, data.frame(t, s, EYa))
  }
}

set.seed(1)
system.time(df_cpp <- st_sweep())
*/


/*
 * 5.2: Fitting the Survival MSM
 * 5.2.1: Discrete Hazard Computation
 * Compute H(s,t) = E[Y(t)a ̄(s)|Y(t 1)a ̄(s) = 0] for s = 0,1,2,3,4, t = 1,2,3,4. 
 * One way is to modify the loops you wrote in the previous question and add H to your data.frame.
 */

// sweep through values of s and t and calculate hazard
// [[Rcpp::export]]
DataFrame st_sweep_hazard(){
  NumericVector s_vec;
  NumericVector t_vec;
  NumericVector EYa_vec;
  for(int s=0; s<5; s++){
    for(int t=1; t<5; t++){
      Rcout << "on iteration s: " << s << ", t: " << t << std::endl;
      //make abar
      NumericVector abar = make_abar(s);
      double EYa = GenerateData_abar_tmax_cpp(1E6,t,abar);
      s_vec.push_back(s);
      t_vec.push_back(t);
      EYa_vec.push_back(EYa);
    }
  }
  
  return(DataFrame::create(_["t"]=t_vec,_["s"]=s_vec,_["EYa"]=EYa_vec));
}