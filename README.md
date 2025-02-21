# DistributionLRS
Code for Variance-Covariance and Distribution of LRS for structured population

In this repository, one will find the code related to the method for the computation of the distribution and variance of LRS

1- the code corresponding the the computations of the (expectation and) variance of the (abundance vectors and) LRS of the illustration can be found
        - for R, in [RVariance.R](https://github.com/ChristopheCoste/DistributionLRS/blob/main/RVariance.R)
        - for Matlab, in https://github.com/ChristopheCoste/DistributionLRS/blob/main/MatlabVariance.m
   This code can be applied to obtain the variance-covariance of LRS for any structured population model

2- the code yielding the pgf of LRS and therefore its distribution, using the symbolic package of Matlab can be found in   https://github.com/ChristopheCoste/DistributionLRS/blob/main/Distribution_Symbolic.m
   There it is applied to the illustration and yields Fig 1 of the main text.
   
   We provide the same code, for illustrative purposes, for a model explicitely embedding different reproduction probability distributions (Bernoulli, Poisson,..) for the different types, in            https://github.com/ChristopheCoste/DistributionLRS/blob/main/Distribution_Symbolic2Poisson.m
   
   Figs 1c and 1d can be obtained without symbolic computation, via the algorithm of Extension I (see below)
   
   
3-  The code for the algorithm of Extension I can be found 
      - for R, in https://github.com/ChristopheCoste/DistributionLRS/blob/main/RVariance.R
      - for Matlab, in https://github.com/ChristopheCoste/DistributionLRS/blob/main/Matlab_Distribution_Numeric.m
    This code can be applied to obtain the variance-covariance of LRS for any structured population model, and is there applied to the illustration, yielding Figs 1c and 1d

4 - The specific computation of the distribution of LRS for Tsuga can be found in https://github.com/ChristopheCoste/DistributionLRS/blob/main/Tsuga_Distribution_Numeric.m, yielding figure 2.

5- The computation of the variance-covaraince of LRS when survival and reproduction covary, as per the illustration of SM2, can be found in https://github.com/ChristopheCoste/DistributionLRS/blob/main/Moments_covaryingFS.m
