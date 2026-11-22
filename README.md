# Bounds of Cycle Integrals of Modular Forms

This repository contains computation scripts, numerical experiments, and theoretical frameworks for analyzing cycle integrals of modular forms evaluated at real quadratic irrationalities. 

This project was conducted as part of the **[UW–Madison Experimental Mathematics Lab (MXM Lab) - Spring 2026](https://mxm.math.wisc.edu/past-semesters/spring-2026/#bounds-of-cycle-integrals-of-modular-forms)** under the supervision of Prof. Yingkun Li and graduate student mentor Yiwen Bai.

## 📌 Background & Motivation

In 2009, Masanobu Kaneko introduced a method to evaluate the Klein $j$-invariant at real quadratic irrationalities $w$ using normalized hyperbolic cycle integrals:

$$val_j(w) = \frac{1}{2 \log \epsilon_w} \int_{\tau}^{A\tau} \frac{j(z) \sqrt{D}}{Q(z, 1)} \, dz$$

Kaneko conjectured—and Bengoechea, Herrero, and Imamoğlu recently proved—that for all real quadratic irrationals $w$, the real part is bounded below by its value at the golden ratio $\phi = \frac{1+\sqrt{5}}{2}$:

$$\mathrm{Re}(val_j(w)) \ge \mathrm{Re}(val_j(\phi))$$

### Core Research Question
Can Kaneko's lower bound be generalized to other weakly holomorphic modular forms of negative weight $-2n$? Specifically, what happens when we replace $j(\tau)$ with forms such as $f = \frac{E_4^2}{\Delta}$ (weight $0$) or $f = \frac{E_4 E_6}{\Delta}$ (weight $-2$)?


## 🔬 Key Experimental Findings

1. **Weight $-2$ ($f = \frac{E_4 E_6}{\Delta}$):**
   - $\mathrm{Re}(val_f(\phi)) = 0$ and $\mathrm{Re}(val_f(\overline{[2,2]})) = 0$ (proved theoretically; numerical artifacts yield $\approx 10^{-14}$).
   - Other values yield both positive and negative results (e.g., $\overline{[1,2]} \mapsto -26.26$, $\overline{[4,1]} \mapsto 34.72$).
   - **Conclusion:** No global minimum or maximum exists for weight $-2$.

2. **Weight $-4$ ($f = \frac{E_4^2}{\Delta}$):**
   - Testing across thousands of hyperbolic matrices shows that the golden ratio $\phi$ acts as an **upper bound** rather than a lower bound:
     - $\phi = \overline{[1]} \implies \approx 99.93$
     - $\overline{[1,2]} \implies \approx 62.47$
     - $\overline{[4,2]} \implies \approx 5.49$
     - $\overline{[100,1]} \implies \approx -19.94$


## 🔍 Analytical Approach Toward the Conjecture

Beyond numerical experiments, we explored an analytical approach to proving the conjectured upper bound for the weight $-4$ form

$$
f(\tau)=\frac{E_4(\tau)^2}{\Delta(\tau)}.
$$

The main idea was to rewrite the difference between two cycle values as an integral of real-valued functions and then investigate conditions under which this integral is nonnegative.

### 1. Reducing the Conjecture to an Integral Inequality

For a modular form $f$ of weight $-2n$ with even $n$ and real Fourier coefficients, define

$$
g_f(t)=e^{-nit}f(e^{it}),
\qquad
t\in I=\left[\frac{\pi}{3},\frac{2\pi}{3}\right].
$$

This function is real-valued and symmetric:

$$
g_f(\pi-t)=g_f(t).
$$

Using the finite-sum formula for cycle integrals, we write

$$
\operatorname{Re}(\operatorname{val}_f(w))
=
\int_I g_f(t)\widehat{S}_w^{(n)}(t)\,dt,
$$

where

$$
\widehat{S}_w^{(n)}(t)
=
\frac{1}{2\log\epsilon_w}
\operatorname{Re}
\left(
i\sum_{j=1}^{N}K_j(t)^{n+1}
\right).
$$

Here, the functions $K_j(t)$ are constructed from the continued fraction representation of the real quadratic irrationality $w$.

Since $g_f$ is symmetric, we can replace $\widehat{S}_w^{(n)}$ with its symmetrization

$$
S_w^{(n)}(t)
=
\frac{1}{2}
\left(
\widehat{S}_w^{(n)}(t)
+
\widehat{S}_w^{(n)}(\pi-t)
\right)
$$

without changing the integral.

We then define the difference function

$$
D_{\phi,w}^{(n)}(t)
=
S_\phi^{(n)}(t)-S_w^{(n)}(t),
$$

where $\phi=(1+\sqrt{5})/2$ is the golden ratio. This gives

$$
\boxed{
\begin{aligned}
&\operatorname{Re}(\operatorname{val}_f(\phi))
-\operatorname{Re}(\operatorname{val}_f(w))\\
&\qquad =
\int_I g_f(t)D_{\phi,w}^{(n)}(t)\,dt.
\end{aligned}
}
$$

Consequently, proving that

$$
\int_I g_f(t)D_{\phi,w}^{(n)}(t)\,dt\geq 0
$$

for every real quadratic irrationality $w$ would establish the conjectured upper bound.

## 🛠 Project Structure
```text
Bounds_of_Cycle_Integrals_of_Modular_forms/
│
├── data/
│   └── src/
│   │   ├── Integration of the j-function.sage   # code used to integrate the j-function
│   ├── data_E_4^2 over delta.csv            # E_4^2/Δ values
│   └── data_E_4^2 over delta.xlsx           # same data, spreadsheet form
│
├── docs/
│   ├── references/                              # papers and reference materials for the project
│   │   ├── Basic Computation with Modular Forms.pdf
│   │   ├── Kaneko_values_j_function.pdf
│   │   ├── REAL PART OF CYCLE INTEGRALS AND CONJECTURES OF KANEKO.pdf
│   │   └── The 1-2-3 of Modular Forms.pdf
│   ├── final_poster.pdf                         # poster presented at the MXM open house
│   ├── final_report.pdf                         # final report
│   └── mid_presentation.pdf                     # slides from the MXM Midsemester Meeting
│
├── .gitignore
├── LICENSE
└── README.md
```
