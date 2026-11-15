import cmath
import math

# 1. Calculate the j-invariant coefficients once and store them in memory 
# (Optimization for faster execution)
j_series = j_invariant_qexp(15)
coeffs = [float(c) for c in j_series.list()]
val = int(j_series.valuation())

print("Loading UI...")

# 2. SageMath-specific UI generator (@interact)
@interact
def calculate_j_integral(
    a = input_box(default=5, label='Matrix a (Top Left):', width=10),
    b = input_box(default=2, label='Matrix b (Top Right):', width=10),
    c = input_box(default=2, label='Matrix c (Bottom Left):', width=10),
    d = input_box(default=1, label='Matrix d (Bottom Right):', width=10),
    re_z0 = input_box(default=-0.5, label='z0 Real Part:', width=10),
    im_z0 = input_box(default=0.5, label='z0 Imaginary Part:', width=10)
):
    # --- Logic that automatically calculates based on input values ---
    
    # Set the input starting point z0
    z0 = complex(re_z0, im_z0)
    
    # Check if the denominator becomes 0 (Safety mechanism)
    if c * z0 + d == 0:
        print("Error: c*z0 + d is 0, so z1 (endpoint) cannot be calculated. Please enter a different z0.")
        return

    # Calculate endpoint z1 and dz_dt
    z1 = (a * z0 + b) / (c * z0 + d)
    dz_dt = z1 - z0
    
    # Calculate discriminant D
    D = (d - a)**2 + 4 * b * c
    
    # Display configuration output
    print("-" * 45)
    print(f"▶ Configured starting point z0: {z0}")
    print(f"▶ Calculated endpoint z1:     {z1}")
    print(f"▶ Discriminant D:             {D}")
    print("-" * 45)
    
    # Handle complex root if D < 0, real root if D > 0
    if D < 0:
        print("Warning: Discriminant D is negative. The result may be a complex number.")
        sqrt_D = cmath.sqrt(D)
    elif D == 0:
        print("Warning: Discriminant D is 0. The denominator becomes 0, and the integral may diverge.")
        return
    else:
        sqrt_D = math.sqrt(D)

    print("Calculating numerical integral... (Please wait)\n")

    # 3. Declare integrand functions (reflecting new a, b, c, d)
    def integrand_real(t):
        z_t = (1.0 - t) * z0 + t * z1
        q_t = cmath.exp(2j * cmath.pi * z_t)
        j_t = sum(coeffs[k] * (q_t**(val + k)) for k in range(len(coeffs)))
        
        Q_val = c*(z_t**2) + (d - a)*z_t - b
        return (j_t * (sqrt_D / Q_val) * dz_dt).real

    def integrand_imag(t):
        z_t = (1.0 - t) * z0 + t * z1
        q_t = cmath.exp(2j * cmath.pi * z_t)
        j_t = sum(coeffs[k] * (q_t**(val + k)) for k in range(len(coeffs)))
        
        Q_val = c*(z_t**2) + (d - a)*z_t - b
        return (j_t * (sqrt_D / Q_val) * dz_dt).imag

    # 4. Execute numerical integration
    re_part, re_err = numerical_integral(integrand_real, 0, 1)
    im_part, im_err = numerical_integral(integrand_imag, 0, 1)

    # 5. Output final results
    print("--- Final Integration Result ---")
    print(f"Real part:      {re_part}")
    print(f"Imaginary part: {im_part}")
    
    # Add a comment if the imaginary part is practically 0 (< 1e-8)
    if abs(im_part) < 1e-8:
        print("\n(Note: The imaginary part is extremely small due to machine precision. Mathematically, this result is a pure real number.)")
