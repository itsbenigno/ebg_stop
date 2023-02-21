#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <array>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <atomic>
#include <random>
#include <cstdlib>

using namespace std;

#define N_THREADS 8

/**
 * Compute the sign function
 * @param a real number
 * @return -1 if the number is negative, 1 if the number is positive
 */
long double get_sign(long double x)
{
  double p = 1.0;
  double n = -1.0;
  if (x >= 0)
    return p;
  else
    return n;
}

/**
 * Calculate the logarithm of a using base b
 * @param number a, for which compute the logarithm, b the base
 */
long double log_a_base_b(long double a, long double b)
{
  return log2(a) / log2(b);
}

/**
 * A class for storing and computing the mean of the variable we are trying to estimate
 * The constructor needs lower and upper bound of the variable
 */
class RandomVariable
{
public:
  long double lower_bound;
  long double upper_bound;
  long double mean = 0;
  long double mean_old = 0;
  long double variance = 0;

  RandomVariable(long double lb, long double ub)
  {
    lower_bound = lb;
    upper_bound = ub;
  }

public:
  void update_mean(long double x, int t)
  {
    mean_old = mean;
    mean = mean_old + ((x - mean_old) / t);
  }

public:
  void update_variance(long double x, int t)
  {
    variance = variance + ((x - mean_old) * (x - mean));
  }
};

/**
 * Execute a simulation of the model
 */
long double simulate_model()
{
  return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

/**
 * Calculate the variable c_t as defined in the EBGstop algorithm
 * @ param the randomVariable for its mean and variance, number of trial t and the sample x
 * @return the value of c_t
 */
long double calculate_ct(RandomVariable randomVariable, int t, long double x)
{
  long double empirical_std_deviation_t = sqrt(randomVariable.variance / t);
  long double R = randomVariable.upper_bound - randomVariable.lower_bound;
  long double c_t = (empirical_std_deviation_t * sqrt((2 * x) / t)) + ((3 * R * x) / t);
  return c_t;
}


long double ebgstop(long double epsilon, long double delta, long double p, long double beta, long double lower_bound, long double upper_bound)
{

  RandomVariable randomVariable(lower_bound, upper_bound);

  long double UB = std::numeric_limits<long double>::max(); //+infinity
  long double LB = 0;
  int t = 0;
  long double k = 0;
  long double x;

  long double sample = simulate_model();

  // check if the termination gave an output
  if (!isnan(sample))
  { 
    t += 1;
    randomVariable.update_mean(sample, t);
    randomVariable.update_variance(sample, t);
  }

  while (LB * (1 + epsilon) < UB * (1 - epsilon))
  {

    double results[N_THREADS];

    // for i in THREAD sample with the ith sample
    #pragma omp parallel for
    for (int i = 0; i < N_THREADS; i++)
    {
      sample = simulate_model();
      results[i] = sample;
    }

    // using previously computed samples to update the mean//
    for (int i = 0; i < N_THREADS; i++)
    {
      sample = results[i];
      if (isnan(sample))
      {
        continue;
      }
      t = t + 1;
      randomVariable.update_mean(sample, t);
      if (t > floor(pow(beta, k)))
      {
        k = k + 1;
        long double alpha = floor(pow(beta, k)) / floor(pow(beta, k - 1));
        long double c = delta * (p - 1) / p;
        long double d_k = c / pow(log_a_base_b(t, beta), p);
        x = -alpha * log(d_k / 3);
      }
      randomVariable.update_variance(sample, t);
      long double c_t = calculate_ct(randomVariable, t, x);
      LB = max(LB, abs(randomVariable.mean) - c_t);
      UB = min(UB, abs(randomVariable.mean) + c_t);

      if (!(LB * (1 + epsilon) < UB * (1 - epsilon)))
      {
        break;
      }
      // cout << "t " << t << " UB " << UB << " Sample " << results[i] << endl << flush;
    }
  }

  long double expected_value = get_sign(randomVariable.mean) * 0.5 * ((1 + epsilon) * LB + (1 - epsilon) * UB);
  /*array<long double, 2> output_mean_t;
  output_mean_t[0] = expected_value;
  output_mean_t[1] = t;*/

  return expected_value;
}

int main(int argc, char **argv)
{

  long double epsilon = 0.001;
  long double delta = 0.001;

  long double lower_bound = 0;
  long double upper_bound = 1;

  long double p = 1.1; //DO NOT CHANGE
  long double beta = 1.1; //DO NOT CHANGE

  long double exp_val = ebgstop(epsilon, delta, p, beta, lower_bound, upper_bound);

  string return_expected_value = to_string(exp_val);
  cout << "expected value: " << return_expected_value << endl;

  return 0;
}
