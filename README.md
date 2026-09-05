# Bounds of Cycle Integrals of Modular Forms

This repository contains computation scripts, numerical experiments, and theoretical frameworks for analyzing cycle integrals of modular forms evaluated at real quadratic irrationalities. 

This project was conducted as part of the **UW–Madison Experimental Mathematics Lab (MXM Lab) - Spring 2026** under the supervision of Prof. Yingkun Li and mentor Yiwen Bai.

---

## 📌 Background & Motivation

In 2009, Masanobu Kaneko introduced a method to evaluate the Klein $j$-invariant at real quadratic irrationalities $w$ using normalized hyperbolic cycle integrals:

$$val_j(w) = \frac{1}{2 \log \epsilon_w} \int_{\tau}^{A\tau} \frac{j(z) \sqrt{D}}{Q(z, 1)} \, dz$$

Kaneko conjectured—and Bengoechea, Herrero, and Imamoğlu recently proved—that for all real quadratic irrationals $w$, the real part is bounded below by its value at the golden ratio $\phi = \frac{1+\sqrt{5}}{2}$:

$$\mathrm{Re}(val_j(w)) \ge \mathrm{Re}(val_j(\phi))$$

### Core Research Question
Can Kaneko's lower bound be generalized to other weakly holomorphic modular forms of negative weight $-2n$? Specifically, what happens when we replace $j(\tau)$ with forms such as $f = \frac{E_4^2}{\Delta}$ (weight $0$) or $f = \frac{E_4 E_6}{\Delta}$ (weight $-2$)?

---

## 🔬 Key Experimental Findings

1. **Weight $-2$ ($f = \frac{E_4 E_6}{\Delta}$):**
   - $\mathrm{Re}(val_f(\phi)) = 0$ and $\mathrm{Re}(val_f(\overline{[2,2]})) = 0$ (proved theoretically; numerical artifacts yield $\approx 10^{-14}$).
   - Other values yield both positive and negative results (e.g., $\overline{[1,2]} \mapsto -26.26$, $\overline{[4,1]} \mapsto 34.72$).
   - **Conclusion:** No global minimum or maximum exists for weight $-2$.

2. **Weight $0$ ($f = \frac{E_4^2}{\Delta}$):**
   - Testing across thousands of hyperbolic matrices shows that the golden ratio $\phi$ acts as an **upper bound** rather than a lower bound:
     - $\phi = \overline{[1]} \implies \approx 99.93$
     - $\overline{[1,2]} \implies \approx 62.47$
     - $\overline{[4,2]} \implies \approx 5.49$
     - $\overline{[100,1]} \implies \approx -19.94$

## 🛠 Project Structure
---
```text
.
├── src/
│   ├── continued_fractions.py  # Matrix representation (T^a V^b) and sequence generators
│   └── cycle_integrals.py      # Core SageMath integration algorithm along u=e^{it}
├── notebooks/
│   └── experimental_data.ipynb # Notebook executing experiments across ~15,000 matrices
├── requirements.txt
└── README.md
