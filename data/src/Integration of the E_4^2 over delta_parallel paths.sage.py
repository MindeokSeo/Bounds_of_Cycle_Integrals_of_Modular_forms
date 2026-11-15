import cmath
import math
import time

# ==========================================
# 1. Initialization: f = E4^2 / Delta (weight -4) and j-function (weight 0)
# ==========================================
print("Initializing Modular Forms and Parameters...")

# Initialize f function (weight -4)
e4 = EisensteinForms(1, 4).q_echelon_basis(60)[0]
e6 = EisensteinForms(1, 6).q_echelon_basis(60)[0]
delta = (e4**3 - e6**2) / 1728
f_func = e4**2 / delta
f_coeffs = [float(k) for k in f_func.list()]
f_val = int(f_func.valuation())

# Initialize j function (weight 0)
j_series = j_invariant_qexp(60)
j_coeffs = [float(c) for c in j_series.list()]
j_val = int(j_series.valuation())

# --- SageMath compatible fundamental domain mapping function (Shared) ---
def get_fundamental_domain_matrix(z):
    p, q, r, s = 1, 0, 0, 1
    curr_z = z
    for _ in range(100):
        rz = float(real_part(curr_z))
        if abs(rz) > 0.5:
            k = round(rz)
            p, q = p - k * r, q - k * s
            curr_z = curr_z - k
        if abs(curr_z) < 0.999999:
            p, q, r, s = -r, -s, p, q
            curr_z = -1.0 / curr_z
        else:
            if abs(float(real_part(curr_z))) <= 0.500001: break
    return p, q, r, s, curr_z

# --- Line segment integration function for f(z) (weight -4) ---
def integral_f_on_segment(start_z, end_z, a, b, c, d, D_factor):
    dz = end_z - start_z
    if abs(dz) < 1e-12: return 0
    
    def segment_integrand(t):
        zt = start_z + t * dz
        p, q, r, s, wt = get_fundamental_domain_matrix(zt)
        q_wt = cmath.exp(2j * cmath.pi * wt)
        f_wt = sum(f_coeffs[k] * (q_wt**(f_val + k)) for k in range(len(f_coeffs)))
        
        # Since f has weight -4, multiply by (r*zt + s)^4 to restore
        f_zt = ((r * zt + s)**4) * f_wt
        Q_val = c*(zt**2) + (d - a)*zt - b
        
        res_complex = f_zt * (D_factor / (Q_val**3)) * dz
        return float(real_part(res_complex))
    
    res, _ = numerical_integral(segment_integrand, 0, 1)
    return res

# --- Line segment integration function for j(z) (weight 0, for validation) ---
def integral_j_on_segment(start_z, end_z, a, b, c, d, sqrt_D):
    dz = end_z - start_z
    if abs(dz) < 1e-12: return 0
    
    def segment_integrand(t):
        zt = start_z + t * dz
        _, _, _, _, wt = get_fundamental_domain_matrix(zt)
        q_wt = cmath.exp(2j * cmath.pi * wt)
        j_wt = sum(j_coeffs[k] * (q_wt**(j_val + k)) for k in range(len(j_coeffs)))
        
        # Since j-function has weight 0, the (r*zt + s) factor is not multiplied (j_zt = j_wt)
        Q_val = c*(zt**2) + (d - a)*zt - b
        
        # Integrand for j-function: (sqrt(D) / Q) * dz
        res_complex = j_wt * (sqrt_D / Q_val) * dz
        return float(real_part(res_complex))
    
    res, _ = numerical_integral(segment_integrand, 0, 1)
    return res

# ==========================================
# 2. Execution setup and exhaustive search for f(z)
# ==========================================
results = []
start_time = time.time()
time_limit = 1.95 
M = 2.0  # Safe altitude (Steeple Path)
max_search = 12

print(f"\nStarting f_func calculation using Steeple Paths (Piecewise Linear)...")

for a in range(1, max_search):
    if time.time() - start_time > time_limit: break
    for b in range(1, max_search):
        for c in range(1, max_search):
            if (1 + b * c) % a == 0:
                d = (1 + b * c) // a
                
                z0 = complex(-d/c, 1/c)
                z1 = (a * z0 + b) / (c * z0 + d)
                
                D = (d - a)**2 + 4 * b * c
                D_factor = D * math.sqrt(D)
                lam = (a + d + math.sqrt(D)) / 2
                log_factor = 2 * math.log(lam)
                
                try:
                    i_raw = (integral_f_on_segment(z0, z0 + M*1j, a, b, c, d, D_factor) +
                             integral_f_on_segment(z0 + M*1j, z1 + M*1j, a, b, c, d, D_factor) +
                             integral_f_on_segment(z1 + M*1j, z1, a, b, c, d, D_factor))
                    
                    val_res = i_raw / log_factor
                    results.append((a, b, c, d, i_raw, val_res))
                except: continue
                
                if time.time() - start_time > time_limit: break

# Print results for f(z)
print("\n" + "=" * 85)
print(f"{'Matrix [a, b, c, d]':<18} | {'I_f(A) (Raw)':<18} | {'val_f(A)':<20}")
print("-" * 85)
for a, b, c, d, i_raw, val_res in results:
    print(f"[{a:2d}, {b:2d}, {c:2d}, {d:2d}]        | {i_raw:>18.4f} | {val_res:>18.6f}")
print("=" * 85)

# ==========================================
# 3. Professor's suggestion: Add j-function validation logic
# ==========================================
print("\nStarting Validation: j-function integrals for [2,1; 1,1] and [5,2; 2,1]")
validation_matrices = [
    (2, 1, 1, 1),
    (5, 2, 2, 1)
]

for a, b, c, d in validation_matrices:
    z0 = complex(-d/c, 1/c)
    z1 = (a * z0 + b) / (c * z0 + d)
    
    D = (d - a)**2 + 4 * b * c
    sqrt_D = math.sqrt(D)
    
    try:
        # Use Steeple Path: Same algorithm as before, same safe altitude (M=2.0)
        j_raw = (integral_j_on_segment(z0, z0 + M*1j, a, b, c, d, sqrt_D) +
                 integral_j_on_segment(z0 + M*1j, z1 + M*1j, a, b, c, d, sqrt_D) +
                 integral_j_on_segment(z1 + M*1j, z1, a, b, c, d, sqrt_D))
        
        print(f"j-function integral for Matrix [{a}, {b}, {c}, {d}]: {j_raw:.6f}")
    except Exception as e:
        print(f"Error calculating j-function integral for Matrix [{a}, {b}, {c}, {d}]: {e}")
print("=" * 85)
