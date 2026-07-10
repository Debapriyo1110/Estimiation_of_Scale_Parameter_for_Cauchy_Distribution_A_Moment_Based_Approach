<p align="center">
<img src="https://capsule-render.vercel.app/api?type=waving&color=0:4F46E5,100:06B6D4&height=250&section=header&text=Estimation%20of%20Scale%20Parameter%20for%20Cauchy%20Distribution&fontSize=32&fontColor=ffffff"/>
</p>

<h1 align="center">
📊 Estimation of Scale Parameter for Cauchy Distribution
</h1>

<h3 align="center">
A Fractional Moment-Based Approach
</h3>

<p align="center">

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Statistics](https://img.shields.io/badge/Statistical%20Inference-Research-blue?style=for-the-badge)
![Monte Carlo](https://img.shields.io/badge/Monte%20Carlo-Simulation-green?style=for-the-badge)
![Heavy Tailed](https://img.shields.io/badge/Heavy--Tailed-Distributions-orange?style=for-the-badge)
![LaTeX](https://img.shields.io/badge/LaTeX-008080?style=for-the-badge&logo=latex)

</p>

---

# 📖 Project Overview

The **Cauchy distribution** is a fundamental heavy-tailed probability distribution whose mean and variance do not exist. Consequently, classical moment-based estimation techniques fail for parameter estimation.

This project proposes a **novel fractional moment-based estimator** for estimating the **scale parameter** of the Cauchy distribution.

The estimator exploits the fact that although integer moments do not exist, **fractional moments of order**

```math
\alpha \in (-1,1)
```

are finite and can therefore be used to construct efficient estimators.

---

# 🎯 Objectives

- Study the statistical challenges associated with the Cauchy distribution.
- Derive a fractional moment-based estimator for the scale parameter.
- Investigate the limiting behaviour as

```math
\alpha \rightarrow 0.
```

- Compare the proposed estimator with classical estimators:
  - Maximum Likelihood Estimator (MLE)
  - Quartile Deviation (QD)

- Evaluate estimator performance through simulation studies.

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
\qquad x\in\mathbb{R}
```

Although integer moments do not exist,

```math
E|X-\mu|^{\alpha}
```

exists for

```math
-1<\alpha<1.
```

This motivates the proposed estimator.

---

# ✨ Proposed Estimator

The proposed estimator is

```math
\hat{\sigma}_{\alpha}
=
\left(
\frac{1}{g(\alpha)n}
\sum_{i=1}^{n}
|X_i-\mu|^{\alpha}
\right)^{1/\alpha},
\qquad
\alpha\in(-1,0)\cup(0,1)
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

the estimator converges to the geometric mean estimator.

---

# 🔬 Methodology

```mermaid
flowchart LR

A[Cauchy Distribution]
-->B[Fractional Moments]

B
-->C[Proposed Estimator]

C
-->D[Monte Carlo Simulation]

D
-->E[Bias Analysis]

D
-->F[MSE Analysis]

E
-->G[Performance Comparison]

F
-->G
```

---

# ⚙️ Simulation Study

| Parameter | Value |
|-----------|--------|
| Sample Size | n = 14 |
| Simulations | 1000 |
| Location Parameter | μ = 5 |
| Scale Parameter | σ = 2 |
| Alpha Values | -0.99 to 0.99 |

---

# 📈 Performance Measures

The estimators were evaluated using:

### 📌 Bias

```math
Bias(\hat{\sigma})
=
E(\hat{\sigma})-\sigma
```

### 📌 Mean Squared Error

```math
MSE(\hat{\sigma})
=
E[(\hat{\sigma}-\sigma)^2]
```

---

# 🏆 Results

The proposed estimator demonstrated:

✅ Lower bias than classical estimators for optimal values of α.

✅ Lower MSE in specific parameter regions.

✅ Greater robustness against outliers and heavy-tailed behaviour.

✅ More concentrated sampling distribution compared to QD and MLE.

---

# 📊 Visualisations

The project includes:

📈 Absolute Bias vs α

📈 MSE vs α

📊 Histograms of QD Estimator

📊 Histograms of MLE Estimator

📊 Histograms of Proposed Estimator

---

# 💻 Software and Tools

### Programming Language

- R Programming

### Packages

- psych
- cauchypca
- graphics

### Documentation

- Overleaf (LaTeX)

---

# 📂 Repository Structure

```text
├── Report/
│   └── Internship_Report.pdf
│
├── R_Code/
│   └── Simulation_Code.R
│
├── Figures/
│   ├── Bias_vs_Alpha.png
│   ├── MSE_vs_Alpha.png
│   └── Histograms.png
│
└── README.md
```

---

# 🚀 Future Scope

- Extension to multivariate Cauchy distributions.
- Application to other heavy-tailed distributions.
- Robust estimation in finance and signal processing.
- Development of generalized fractional moment estimators.

---

# 👨‍🎓 Author

## Debapriyo Bhar

B.Sc. Major in Statistics

Ramakrishna Mission Residential College (Autonomous), Narendrapur

📧 Email:
debapriyobhar@gmail.com

💼 LinkedIn:
https://www.linkedin.com/in/debapriyo-bhar-5074a6303

---

# 🎓 Project Guide

## Prof. Sabyasachi Bhattacharya

Professor

Applied Statistics and Econometrics Research Unit (AERU)

Indian Statistical Institute (ISI), Kolkata

---

<p align="center">

⭐ If you found this project interesting, consider giving it a star!

</p>

<p align="center">
<img src="https://capsule-render.vercel.app/api?type=waving&color=0:06B6D4,100:4F46E5&height=120&section=footer"/>
</p>
