Estimation of Scale Parameter for Cauchy Distribution

📌 Overview

This project proposes a fractional moment-based estimator for the scale parameter (σ) of the Cauchy distribution — a heavy-tailed distribution with undefined mean and variance.
The method is compared with Quartile Deviation (QD) and Maximum Likelihood Estimator (MLE) via R simulations.

🎯 Key Points

Uses fractional moments (α ∈ (-1,0) ∪ (0,1)) for robust estimation.

Outperforms QD and MLE in terms of bias and MSE.

Applicable to other heavy-tailed distributions.

⚙️ Method

We have Derived an Estimator using method of moments estimation
 
Simulate Cauchy data (µ=5, σ=2) in R.

Compare bias & MSE for proposed, QD, and MLE estimators.

📊 Results

Best α ≈ 0.06 for minimum bias.

Proposed estimator: narrow, symmetric distribution around true σ.

QD: higher MSE; MLE: sensitive to outliers.

Author: Debapriyo Bhar, Arkabrata Mondal

Guide: Prof. Sabyasachi Bhattacharya, Professor, ISI Kolkata
