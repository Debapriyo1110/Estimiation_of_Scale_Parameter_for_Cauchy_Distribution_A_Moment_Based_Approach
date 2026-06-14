# Estimation of Scale Parameter for Cauchy Distribution

This project proposes a **fractional moment-based estimator** for the scale parameter of the Cauchy distribution, a challenging problem due to the absence of finite integer moments.

## Overview

* Derived a new estimator using fractional moments of order α ∈ (-1,1).
* Studied its limiting behaviour as α → 0.
* Compared the proposed estimator with:

  * Maximum Likelihood Estimator (MLE)
  * Quartile Deviation (QD)

## Methodology

* Simulated 1000 Cauchy-distributed samples in **R**
* Evaluated estimator performance using:

  * Bias
  * Mean Squared Error (MSE)

## Results

The proposed estimator showed:

* Lower bias than classical estimators for optimal α
* Lower MSE in specific parameter ranges
* Greater robustness against outliers and heavy-tailed behaviour

## Tools Used

* R Programming
* Overleaf (LaTeX)
* Packages: psych, cauchypca

## Future Scope

* Extension to multivariate Cauchy distributions
* Application to other heavy-tailed distributions
* Real-world use in finance, signal processing, and robust data analysis
