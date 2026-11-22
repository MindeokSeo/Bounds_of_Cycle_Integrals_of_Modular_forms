# Bounds of Cycle Integrals of Modular Forms

This repository contains computation scripts, numerical experiments, and theoretical frameworks for analyzing cycle integrals of modular forms evaluated at real quadratic irrationalities. 

This project was conducted as part of the **[UW–Madison Experimental Mathematics Lab (MXM Lab) - Spring 2026](https://mxm.math.wisc.edu/past-semesters/spring-2026/#bounds-of-cycle-integrals-of-modular-forms)** under the supervision of Prof. Yingkun Li and graduate student mentor Yiwen Bai.

## 📌 Background & Motivation

In 2009, Masanobu Kaneko introduced a method to evaluate the Klein $j$-invariant at real quadratic irrationalities $w$ using normalized hyperbolic cycle integrals:

$$\mathrm{val}_j(w) = \frac{1}{2 \log \epsilon_w} \int_{\tau}^{A\tau} \frac{j(z) \sqrt{D}}{Q(z, 1)} \, dz$$

Kaneko conjectured—and Bengoechea, Herrero, and Imamoğlu recently proved—that for all real quadratic irrationals $w$, the real part is bounded below by its value at the golden ratio $\phi = \frac{1+\sqrt{5}}{2}$:

$$\mathrm{Re}\big(\mathrm{val}_j(w)\big) \ge \mathrm{Re}\big(\mathrm{val}_j(\phi)\big)$$

### Core Research Question
Can Kaneko's lower bound be generalized to other weakly holomorphic modular forms of negative weight $-2n$? Specifically, what happens when we replace $j(\tau)$ with forms such as $f = \frac{E_4^2}{\Delta}$ (weight $-4$) or $f = \frac{E_4 E_6}{\Delta}$ (weight $-2$)?


## 🔬 Key Experimental Findings

1. **Weight $-2$ ($f = \frac{E_4 E_6}{\Delta}$):**
   - $\mathrm{Re}(\mathrm{val}_f(\phi)) = 0$ and $\mathrm{Re}(\mathrm{val}_f(\overline{[2,2]})) = 0$ (proved theoretically; numerical artifacts yield $\approx 10^{-14}$).
   - Other values yield both positive and negative results (e.g., $\overline{[1,2]} \mapsto -26.26$, $\overline{[4,1]} \mapsto 34.72$).
   - **Conclusion:** No global minimum or maximum exists for weight $-2$.

2. **Weight $-4$ ($f = \frac{E_4^2}{\Delta}$):**
   - Testing across thousands of hyperbolic matrices shows that the golden ratio $\phi$ acts as an **upper bound** rather than a lower bound:
     - $\phi = \overline{[1]}$ gives $\approx 99.93$
     - $\overline{[1,2]}$ gives $\approx 62.47$
     - $\overline{[4,2]}$ gives $\approx 5.49$
     - $\overline{[100,1]}$ gives $\approx -19.94$


## 🧩 Attempted Proof of the Weight $-4$ Conjecture

> **Status: not proven.** This section records the proof strategy we pursued (Section 5 of [`docs/final_report.pdf`](docs/final_report.pdf); the *A Useful Lemma* and *Extending the Lemma* panels of [`docs/final_poster.pdf`](docs/final_poster.pdf)), which ingredients we managed to establish, and exactly where the argument gets stuck. The project provides numerical evidence and an analytical *reduction* of the conjecture, but it does not establish the bound.

### The conjecture

**Conjecture 4.1.** For $f = E_4^2/\Delta$ (weight $-4$, so $n = 2$) and every real quadratic irrational $w$,

$$\mathrm{Re}\big(\mathrm{val}_f(w)\big) \le \mathrm{Re}\big(\mathrm{val}_f(\phi)\big) = 99.92938911748521\ldots$$

**Normalization.** All numerical values use

$$\Delta(\tau) = \frac{E_4(\tau)^3 - E_6(\tau)^2}{1728} = q\prod_{m\ge 1}(1-q^m)^{24}, \qquad q = e^{2\pi i \tau},$$

i.e. *without* a $(2\pi)^{12}$ factor, with the normalized Eisenstein series

$$E_{2k}(\tau) = \frac{1}{2\zeta(2k)}\sum_{(m,n)\in\mathbb{Z}^2\setminus\{(0,0)\}} \frac{1}{(m+n\tau)^{2k}}.$$

### Strategy at a glance

1. Rewrite $\mathrm{val}_f(w)$ as an integral over the arc $e^{it}$, $\frac{\pi}{3}\le t\le\frac{2\pi}{3}$, against an explicit finite sum built from the continued fraction of $w$.
2. For even $n$, split the integrand into a modular-form part $g$ (independent of $w$) and a combinatorial part $S_w^n$ (depends on $w$), and reduce the conjecture to the **sign of one integral** $\int g\,D_{\phi,w}^n\,dt$.
3. Try to prove that sign with a **Chebyshev-type comparison lemma**.
4. Check the lemma's hypotheses on $g$ (via Ramanujan's derivative formulas) and on $D_{\phi,w}^n$.

### Step 1 — A finite-sum integral formula

Every real quadratic irrational has an $SL_2(\mathbb{Z})$-equivalent with a purely periodic continued fraction $w = \overline{[a_1,\dots,a_{2\ell}]}$, and $\mathrm{val}_f$ is $SL_2(\mathbb{Z})$-invariant, so we may work with it. Put

$$T = \begin{pmatrix}1&1\\0&1\end{pmatrix}, \qquad V = \begin{pmatrix}1&0\\1&1\end{pmatrix},$$

and

$$B = T^{a_1}V^{a_2}\cdots T^{a_{2\ell-1}}V^{a_{2\ell}} = L_1L_2\cdots L_N, \qquad L_j\in\{T,V\}, \quad N = a_1+\cdots+a_{2\ell}.$$

With $P_j = L_1\cdots L_j$ ($P_0 = I$), set $u_j = P_{j-1}^{-1}w$, $v_j = P_{j-1}^{-1}w'$ ($w'$ the conjugate of $w$) and

$$K_j(t) = \frac{e^{it}}{e^{it}-u_j} - \frac{e^{it}}{e^{it}-v_j}.$$

Then, for a weakly holomorphic form $f$ of weight $-2n$,

$$\mathrm{val}_f(w) = \frac{1}{2\log\epsilon_w}\int_{\pi/3}^{2\pi/3} f(e^{it})\, i\, e^{-nit}\sum_{j=1}^{N} K_j(t)^{n+1}\,dt.$$

This is the formula of Bengoechea–Herrero–Imamoğlu, which extends directly from weight $0$ to weight $-2n$, and it is what the SageMath code computes.

### Step 2 — Reduction to a sign problem

Let $I = [\frac{\pi}{3},\frac{2\pi}{3}]$ and $J = [\frac{\pi}{3},\frac{\pi}{2}]$. For **$n$ even** and $f$ with real Fourier coefficients (true for $E_4^2/\Delta$), define

$$g_f(t) = e^{-nit}f(e^{it}),$$

which is real on $I$ and satisfies $g_f(\pi - t) = g_f(t)$. Define

$$\widehat{S}_w^n(t) = \frac{1}{2\log\epsilon_w}\mathrm{Re}\Big(i\sum_{j=1}^{N}K_j(t)^{n+1}\Big),$$

so that

$$\mathrm{Re}\big(\mathrm{val}_f(w)\big) = \int_I g_f(t)\,\widehat{S}_w^n(t)\,dt.$$

Since $g_f$ is symmetric about $\pi/2$, we may replace $\widehat{S}_w^n$ by its symmetrization

$$S_w^n(t) = \frac{1}{2}\Big(\widehat{S}_w^n(t) + \widehat{S}_w^n(\pi - t)\Big)$$

without changing the integral, so that $S_w^n(\pi - t) = S_w^n(t)$. With $D_{\phi,w}^n = S_\phi^n - S_w^n$ (here $n$ is a superscript, not an exponent),

$$\mathrm{Re}\big(\mathrm{val}_f(\phi)\big) - \mathrm{Re}\big(\mathrm{val}_f(w)\big) = \int_I g_f(t)\,D_{\phi,w}^n(t)\,dt.$$

So the conjecture follows if this integral is $\ge 0$ for **every** real quadratic irrational $w$.

### Step 3 — A comparison lemma

The classical Chebyshev-type argument needs the mean-zero hypothesis $\int_I D\,dt = 0$. This holds for $n = 0$ (because $\mathrm{val}_1(w) = 1$ for every $w$) but is not available for $n\ne 0$ (final report, Sec. 5.2), so the original lemma cannot be applied to $D_{\phi,w}^2$ directly. For our case the poster uses the following variant instead:

> **Variant (final poster).** If $g$ and $D$ are both nonpositive and decreasing on $J = [\frac{\pi}{3},\frac{\pi}{2}]$ (and symmetric about $\frac{\pi}{2}$), then $\int_I g\,D\,dt \ge 0$.

This variant does not use the mean-zero hypothesis. The final report states only the original lemma (Lemma 5.1) and lists the replacement condition as an open question (Question 2), so the variant above should be read as the poster's proposed replacement.

### Step 4 — The function $g$ for $f = E_4^2/\Delta$

For $n = 2$,

$$g(t) = e^{-2it}\,\frac{E_4(e^{it})^2}{\Delta(e^{it})}.$$

**Real and symmetric.** $g$ is real on $I$ and $g(\pi - t) = g(t)$.

**Nonpositive (numerical).** The numerical plot (Figure 1 in the report) strongly suggests $g\le 0$ on $I$, vanishing at the endpoints (where $E_4(e^{it}) = 0$) and attaining its minimum at $t = \pi/2$.

**Monotonicity via Ramanujan's formulas.** Ramanujan's formulas read

$$\frac{1}{2\pi i}\frac{dE_4}{d\tau} = \frac{E_2E_4 - E_6}{3}, \qquad \frac{1}{2\pi i}\frac{d\Delta}{d\tau} = E_2\,\Delta,$$

where $E_2 = 1 - 24\sum_{m\ge 1}\sigma_1(m)q^m$ is the quasimodular Eisenstein series of weight two. They give

$$\frac{f'(\tau)}{f(\tau)} = -\frac{2\pi i}{3}\left(E_2(\tau) + 2\,\frac{E_6(\tau)}{E_4(\tau)}\right),$$

and therefore, for $\frac{\pi}{3} < t < \frac{2\pi}{3}$,

$$\frac{g'(t)}{g(t)} = -2i + \frac{2\pi}{3}\,e^{it}\left(E_2(e^{it}) + 2\,\frac{E_6(e^{it})}{E_4(e^{it})}\right).$$

Since $g$ is real, the left side is real, so we may take real parts. If $g<0$ in the interior of $I$, then $g'\le 0$ on $(\frac{\pi}{3},\frac{\pi}{2})$ as soon as

$$\mathrm{Re}\left(e^{it}\left[E_2(e^{it}) + 2\,\frac{E_6(e^{it})}{E_4(e^{it})}\right]\right) \ge 0 \quad \text{for } \tfrac{\pi}{3} < t < \tfrac{\pi}{2}. \qquad (\ast)$$

### Where things stand

The reduction in Step 2 is exact. Applying the poster's variant to $g$ and $D_{\phi,w}^2$ requires the following:

| Ingredient | Status |
|---|---|
| $g$ real on $I$ and $g(\pi-t) = g(t)$ | shown in the report ($n$ even, real coefficients) |
| $D_{\phi,w}^2(\pi-t) = D_{\phi,w}^2(t)$ | by construction of $S_w^n$ |
| $g\le 0$ on $I$ | numerical evidence only (Figure 1) |
| $g$ nonincreasing on $J$ | reduced to inequality $(\ast)$, which is unproven |
| $D_{\phi,w}^2$ nonpositive and nonincreasing on $J$ | not written up in the report |

### Open questions

1. What techniques give a rigorous proof that $g\le 0$ on $I$ and that $g$ is nonincreasing on $J$, for instance by proving $(\ast)$?
2. Which conditions on $D_{\phi,w}^2$ can replace the mean-zero hypothesis in the comparison lemma? (The poster's variant in Step 3 is one proposal; its hypotheses on $D_{\phi,w}^2$ are not verified.)
3. How can the conjecture be refined for odd $n$?
4. Does this method work for other modular forms of weight $-2n$ with $n$ even?

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

## 📚 References

1. P. Bengoechea, S. Herrero, Ö. İmamoğlu. *A Lyapunov exponent attached to modular functions.* 2025. [arXiv:2503.16343](https://arxiv.org/abs/2503.16343).
2. P. Bengoechea, S. Herrero, Ö. İmamoğlu. *Real part of cycle integrals and conjectures of Kaneko.* 2025. [arXiv:2505.14500](https://arxiv.org/abs/2505.14500).
3. M. Kaneko. *Observations on the "values" of the elliptic modular function $j(\tau)$ at real quadratics.* Kyushu J. Math. 63(2) (2009), 353–364.
