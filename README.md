# Empirical Bernstein stopping
This is a C++ implementation of the Empirical Bernstein stopping algorithm (in particular EBGstop, a faster version of it).

Mnih, Volodymyr & Szepesvári, Csaba & Audibert, Jean-Yves. (2008). [Empirical Bernstein stopping](https://dl.acm.org/doi/10.1145/1390156.1390241)

## Prerequisites
Make sure that your C++ compiler supports OpenMP

##Run 
N_THREADS defines how many parallel simulation of your random variable to run.
Modify the simulate_model() function with the random variable you want to estimate, then open a terminal inside the project directory and write 
```
make
```
to compile and run the code.


