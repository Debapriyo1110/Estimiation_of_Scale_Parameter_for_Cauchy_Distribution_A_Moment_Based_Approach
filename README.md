<p align="center">
<img src="https://capsule-render.vercel.app/api?type=waving&height=260&color=0:1A2980,50:26D0CE,100:00C9FF&text=Estimation%20of%20Scale%20Parameter%20for%20the%20Cauchy%20Distribution&fontColor=ffffff&fontSize=30&animation=fadeIn"/>
</p>

<h1 align="center">
📊 Estimation of Scale Parameter for the Cauchy Distribution
</h1>

<h3 align="center">
A Fractional Moment-Based Approach
</h3>

<p align="center">

<b>Research Internship Project</b><br>

Applied Statistics & Econometrics Research Unit (AERU)<br>

Indian Statistical Institute (ISI), Kolkata

</p>

<p align="center">

<img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white"/>

<img src="https://img.shields.io/badge/Statistics-Research-blue?style=for-the-badge"/>

<img src="https://img.shields.io/badge/Monte%20Carlo-Simulation-success?style=for-the-badge"/>

<img src="https://img.shields.io/badge/Heavy--Tailed-Distributions-orange?style=for-the-badge"/>

<img src="https://img.shields.io/badge/LaTeX-008080?style=for-the-badge&logo=latex"/>

<img src="https://img.shields.io/badge/License-MIT-green?style=for-the-badge"/>

</p>

---

# 🌟 Project Highlights

- 🔹 Novel **Fractional Moment-Based Estimator** for the Cauchy scale parameter
- 🔹 Research on **Heavy-Tailed Statistical Distributions**
- 🔹 Extensive **Monte Carlo Simulation Study**
- 🔹 Comparative analysis with **Maximum Likelihood Estimation (MLE)** and **Quartile Deviation (QD)**
- 🔹 Performance evaluation using **Bias** and **Mean Squared Error (MSE)**
- 🔹 Implemented entirely in **R Programming**

---

# 📖 Project Overview

The **Cauchy distribution** is a classical heavy-tailed probability distribution whose mean and variance are undefined. As a result, conventional moment-based estimation techniques cannot be applied directly.

This research proposes a **fractional moment-based estimator** for estimating the **scale parameter** of the Cauchy distribution. The proposed approach exploits the existence of fractional moments over the interval

```math
-1 < \alpha < 1,
```

allowing efficient estimation where ordinary moments fail.

The estimator is evaluated through an extensive Monte Carlo simulation study and compared with well-known classical estimators, including **Maximum Likelihood Estimation (MLE)** and **Quartile Deviation (QD)**.

---

# 🎯 Objectives

- Investigate the estimation problem for the Cauchy scale parameter.
- Develop a fractional moment-based estimator.
- Study the limiting behaviour as **α → 0**.
- Compare the proposed estimator with MLE and Quartile Deviation.
- Evaluate estimator performance using Monte Carlo simulations.
- Analyse statistical performance through Bias and Mean Squared Error (MSE).

---

# 🧮 Mathematical Background

For

```math
X \sim C(\mu,\sigma),
```

the probability density function is

```math
f(x;\mu,\sigma)=
\frac{1}{\pi\sigma}
\frac{1}
{1+\left(\frac{x-\mu}{\sigma}\right)^2},
\qquad x\in\mathbb{R}.
```

Although the Cauchy distribution has **no finite mean or variance**, its fractional moments satisfy

```math
E|X-\mu|^{\alpha},
\qquad -1<\alpha<1,
```

which forms the theoretical foundation of the proposed estimator.

---

# 🔬 Research Workflow

```mermaid
flowchart LR

A[Cauchy Distribution]

-->B[Fractional Moments]

B-->C[Proposed Estimator]

C-->D[Monte Carlo Simulation]

D-->E[Bias Analysis]

D-->F[MSE Analysis]

E-->G[Performance Comparison]

F-->G

G-->H[Results & Conclusions]

```

---

# ⚙️ Simulation Configuration

| Parameter | Value |
|-----------|--------|
| Sample Size | **n = 14** |
| Number of Simulations | **1000** |
| Location Parameter | **μ = 5** |
| Scale Parameter | **σ = 2** |
| Fractional Order | **−0.99 ≤ α ≤ 0.99** |

---

# 📈 Proposed Estimator

The proposed estimator for the scale parameter is given by

```math
\hat{\sigma}_{\alpha}
=
\left(
\frac{1}{g(\alpha)n}
\sum_{i=1}^{n}
|X_i-\mu|^{\alpha}
\right)^{1/\alpha},
\qquad
\alpha\in(-1,0)\cup(0,1),
```

where

```math
g(\alpha)
=
\frac{
\Gamma\left(\frac{1+\alpha}{2}\right)
\Gamma\left(\frac{1-\alpha}{2}\right)
}{\pi}
=
\frac{1}{\cos(\alpha\pi/2)}.
```

As

```math
\alpha \rightarrow 0,
```

the estimator converges to the geometric mean estimator, providing a natural connection between fractional moment estimation and classical robust estimation techniques.

---

# 📊 Performance Evaluation

The proposed estimator was assessed using extensive Monte Carlo simulations and compared with two widely used classical estimators.

### Compared Estimators

- Maximum Likelihood Estimator (MLE)
- Quartile Deviation (QD)
- Proposed Fractional Moment Estimator

### Performance Measures

- **Bias**

```math
Bias(\hat{\sigma})=E(\hat{\sigma})-\sigma
```

- **Mean Squared Error (MSE)**

```math
MSE(\hat{\sigma})
=
E[(\hat{\sigma}-\sigma)^2]
```

---

# 🏆 Key Findings

The simulation study demonstrates that the proposed estimator:

- ✅ Achieves lower bias for suitable values of **α**.
- ✅ Produces lower Mean Squared Error in several parameter settings.
- ✅ Is robust against the heavy-tailed behaviour of the Cauchy distribution.
- ✅ Performs competitively with Maximum Likelihood Estimation (MLE).
- ✅ Outperforms the Quartile Deviation estimator in terms of efficiency for many cases.
- ✅ Converges smoothly to the geometric mean estimator as **α → 0**.

---

# 📊 Included Visualisations

The repository contains graphical analyses including:

- 📈 Absolute Bias vs. α
- 📈 Mean Squared Error vs. α
- 📊 Histogram of the Proposed Estimator
- 📊 Histogram of the Maximum Likelihood Estimator
- 📊 Histogram of the Quartile Deviation Estimator

---

# 🛠 Software & Technologies

### Programming Language

- **R**

### R Packages

- psych
- cauchypca
- graphics

### Documentation

- LaTeX (Overleaf)

### Statistical Techniques

- Fractional Moment Estimation
- Monte Carlo Simulation
- Bias Analysis
- Mean Squared Error Analysis
- Robust Statistical Estimation

---

# 📂 Repository Structure

```text
📦 Estimation-of-Scale-Parameter-for-Cauchy-Distribution
│
├── 📄 Certificate.pdf
├── 📘 Final Report.pdf
├── 📊 Internship_new.R
├── 📜 LICENSE
└── 📖 README.md
```

---

# 🔬 Research Contribution

This work introduces a fractional moment-based estimator for estimating the scale parameter of the Cauchy distribution by exploiting the existence of finite fractional moments.

A comprehensive Monte Carlo simulation study was conducted to investigate its statistical properties and compare its performance with established estimation methods using Bias and Mean Squared Error (MSE) as evaluation criteria.

---

# 🚀 Future Scope

Potential directions for future research include:

- Extension to multivariate Cauchy distributions.
- Generalization to other heavy-tailed probability distributions.
- Applications in financial risk modelling.
- Robust parameter estimation in signal processing.
- Bayesian approaches for heavy-tailed distributions.
- Integration with machine learning methods for robust statistical inference.

---

# 👨‍🎓 Author

**Debapriyo Bhar**

**B.Sc. Major in Statistics**

Ramakrishna Mission Residential College (Autonomous), Narendrapur

📧 **Email:**  
debapriyobhar@gmail.com

💼 **LinkedIn:**  
https://www.linkedin.com/in/debapriyo-bhar-5074a6303

---

# 🎓 Research Supervisor

**Prof. Sabyasachi Bhattacharya**

Professor

Applied Statistics & Econometrics Research Unit (AERU)

Indian Statistical Institute (ISI), Kolkata

---

# 🙏 Acknowledgements

I sincerely express my gratitude to **Prof. Sabyasachi Bhattacharya** for his invaluable guidance, continuous encouragement, and insightful suggestions throughout this research internship. I also acknowledge the **Applied Statistics & Econometrics Research Unit (AERU), Indian Statistical Institute (ISI), Kolkata**, for providing an excellent academic environment to conduct this work.

---

<p align="center">
<img src="https://capsule-render.vercel.app/api?type=waving&height=120&color=0:00C9FF,50:26D0CE,100:1A2980&section=footer"/>
</p>
