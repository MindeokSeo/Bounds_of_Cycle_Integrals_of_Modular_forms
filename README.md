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
\operatorname{Re}(\operatorname{val}_f(\phi))
-
\operatorname{Re}(\operatorname{val}_f(w))
=
\int_I g_f(t)D_{\phi,w}^{(n)}(t)\,dt.
}
$$

Consequently, proving that

$$
\int_I g_f(t)D_{\phi,w}^{(n)}(t)\,dt\geq 0
$$

for every real quadratic irrationality $w$ would establish the conjectured upper bound.

### 2. Investigating the Sign and Monotonicity of $g_f$

For our weight $-4$ example, $n=2$, so

$$
g(t)
=
e^{-2it}
\frac{E_4(e^{it})^2}{\Delta(e^{it})}.
$$

Numerical plots suggest that:

- $g(t)$ is nonpositive on $I$ and vanishes at the endpoints.
- $g(t)$ is symmetric about $t=\pi/2$.
- $g(t)$ is nonincreasing on $J=[\pi/3,\pi/2]$.

To investigate these properties analytically, we used Ramanujan's differential identities for Eisenstein series:

$$
\frac{1}{2\pi i}\frac{dE_4}{d\tau}
=
\frac{E_2E_4-E_6}{3},
\qquad
\frac{1}{2\pi i}\frac{d\Delta}{d\tau}
=
E_2\Delta.
$$

Differentiating $g(t)$ gives

$$
\frac{g'(t)}{g(t)}
=
-2i+
\frac{2\pi}{3}e^{it}
\left(
E_2(e^{it})
+
2\frac{E_6(e^{it})}{E_4(e^{it})}
\right).
$$

Since $g(t)$ is real, taking real parts yields

$$
\frac{g'(t)}{g(t)}
=
\frac{2\pi}{3}
\operatorname{Re}
\left[
e^{it}
\left(
E_2(e^{it})
+
2\frac{E_6(e^{it})}{E_4(e^{it})}
\right)
\right].
$$

Assuming $g(t)<0$ in the interior of $I$, it would therefore be sufficient to establish

$$
\operatorname{Re}
\left[
e^{it}
\left(
E_2(e^{it})
+
2\frac{E_6(e^{it})}{E_4(e^{it})}
\right)
\right]
\geq 0
$$

for $\pi/3<t<\pi/2$ in order to show that $g$ is nonincreasing on $J$.

The numerical graph supports these properties, but the required inequalities were not proved.

### 3. Exploring an Integral Comparison Lemma

We also investigated whether a comparison principle could establish the desired sign of the integral.

We considered the following form of Chebyshev's integral inequality. Suppose that two continuous functions $g,D:I\to\mathbb{R}$ satisfy

$$
g(\pi-t)=g(t),
\qquad
D(\pi-t)=D(t),
$$

and that $g$ is nondecreasing while $D$ is nonincreasing on $J$. If

$$
\int_I D(t)\,dt=0,
$$

then

$$
\int_I g(t)D(t)\,dt\leq 0.
$$

However, this lemma cannot be applied directly to our weight $-4$ case.

Indeed, when $n=0$, applying the cycle-integral formula to the constant function $f=1$ gives

$$
\operatorname{val}_1(w)=1.
$$

This implies

$$
\int_I S_w^{(0)}(t)\,dt=1
$$

for every $w$, and consequently

$$
\int_I D_{\phi,w}^{(0)}(t)\,dt=0.
$$

For $n\neq 0$, this mean-zero identity does not generally hold. Thus, the comparison lemma's hypothesis fails for $D_{\phi,w}^{(2)}$.

As an alternative, we explored a sign-based approach, suggested in our final poster: if both $g(t)$ and $D_{\phi,w}^{(2)}(t)$ could be shown to be nonpositive on $J$, their symmetry would yield a nonnegative integrand on $I$. However, establishing the required sign of $D_{\phi,w}^{(2)}$ for arbitrary $w$ remains an open step in this approach.

### 4. Remaining Challenges

Our analytical investigation reduced the conjecture to a sign problem involving two functions, $g(t)$ and $D_{\phi,w}^{(2)}(t)$. However, several questions remain unresolved:

1. **Proving the properties of $g$:** Establish rigorously that $g(t)\leq 0$ on $I$ and that it is nonincreasing on $J$, potentially by proving the inequality derived from Ramanujan's differential identities.

2. **Controlling the difference function:** Find suitable conditions on $D_{\phi,w}^{(2)}(t)$ that guarantee the desired sign of the integral, despite the failure of the mean-zero condition.

3. **Extending the method:** Investigate whether this analytical framework can be adapted to modular forms of other negative even weights and whether a refined version can handle odd $n$.

Our work provides a numerical conjecture and an analytical reduction toward its proof, but does not establish the conjectured upper bound.


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
