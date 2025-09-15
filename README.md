This repositery contains the code related to the various computations in the article "A Branching Process Approach to Lifetime Reproductive Success of Structured Populations: Variance-Covariance and Distribution" published in Methods in Ecology and Evolution (2025) by Christophe F.D. Coste
It contains code in both R and Matlab. 

**(I)** the code corresponding the the computations of the (expectation and) **variance of** the (abundance vectors and) **LRS** (Lifetime reproductive success) of the illustration can be found (a) for R, in RVariance.R and (b) for Matlab, in [MatlabVariance.m](https://github.com/ChristopheCoste/DistributionLRS/blob/main/MatlabVariance.m)
- this code can be applied to obtain the **variance-covariance of LRS for any structured population model**

**(IIa)** the code yielding the **pgf** (probability generating function) **of LRS** for the illustration and, from it, the **joint distribution of LRS**, using the symbolic package of Matlab can be found in [Distribution_Symbolic.m](https://github.com/ChristopheCoste/DistributionLRS/blob/main/Distribution_Symbolic.m)
- It yields Fig 1 of the main text.
  
**(IIb)**  We provide a similar code, applicable to any projection model, providing the probability of a given LRS vector (reproduction can have any distribution, such as Bernoulli, Poisson, etc ...) , in [Distribution_Symbolic_General.m](https://github.com/ChristopheCoste/DistributionLRS/blob/main/Distribution_Symbolic_General.m)
- It yields Fig 1 of the main text (Note that Figs 1c and 1d can be obtained without symbolic computation, via the algorithm of Extension I,see below)

**(III)** the code for the algorithm of Extension I can be found, for R as  [Rdistribution.R](https://github.com/ChristopheCoste/DistributionLRS/blob/main/Rdistribution.R) and, for Matlab, as [Matlab_Distribution_Numeric.m](https://github.com/ChristopheCoste/DistributionLRS/blob/main/Matlab_Distribution_Numeric.m) 
- This code can be applied to obtain the **distribution of total LRS for any structured population model**, and is there applied to the illustration, yielding Figs 1c and 1d.


**(IV)** the specific computation of the distribution of LRS for Tsuga (extension II) can be found in [Tsuga_Distribution_Numeric.m](https://github.com/ChristopheCoste/DistributionLRS/blob/main/Tsuga_Distribution_Numeric.m) , yielding figure 2.

**(V)** the computation of the variance-covaraince of LRS when survival and reproduction covary, as per the illustration of SM2, can be found in  [Moments_covaryingFS.m](https://github.com/ChristopheCoste/DistributionLRS/blob/main/Moments_covaryingFS.m)
