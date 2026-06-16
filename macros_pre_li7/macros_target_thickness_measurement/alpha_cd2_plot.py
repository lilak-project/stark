"""
Plot: 4He (alpha) 3.271 MeV remaining energy vs CD2 thickness.
Reads thickness_and_alpha_energy_remaing.dat if present,
otherwise recomputes from stopping power tables.
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os

# ── Stopping power (NIST ASTAR, Bragg additivity) ────────────────────────────
nist_C_E = np.array([
    0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.30,
    0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
    0.80, 0.85, 0.90, 0.95, 1.00, 1.25, 1.50, 1.75, 2.00,
    2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50,
    6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00,
])
nist_C_S = np.array([
    4.918, 5.051, 5.097, 5.094, 4.999, 4.870, 4.727, 4.577, 4.426,
    4.130, 3.858, 3.617, 3.405, 3.218, 3.052, 2.904, 2.771, 2.650,
    2.541, 2.442, 2.352, 2.269, 2.194, 1.898, 1.682, 1.517, 1.389,
    1.286, 1.202, 1.131, 1.072, 0.975, 0.899, 0.838, 0.789, 0.747,
    0.712, 0.681, 0.654, 0.629, 0.608, 0.588, 0.570, 0.554, 0.540,
])
nist_H_E = nist_C_E.copy()
nist_H_S = np.array([
    15.91, 16.18, 15.96, 15.51, 14.93, 14.29, 13.63, 12.97, 12.34,
    11.16, 10.10, 9.181, 8.386, 7.699, 7.104, 6.586, 6.132, 5.732,
    5.378, 5.063, 4.781, 4.528, 4.298, 3.373, 2.749, 2.309, 1.990,
    1.748, 1.559, 1.408, 1.284, 1.086, 0.942, 0.832, 0.746, 0.678,
    0.622, 0.576, 0.537, 0.503, 0.473, 0.448, 0.425, 0.405, 0.387,
])

def extend_low(E_t, S_t, E_min=0.001):
    slope = np.log(S_t[1]/S_t[0]) / np.log(E_t[1]/E_t[0])
    E_l = np.geomspace(E_min, E_t[0]*0.99, 30)
    return np.concatenate([E_l, E_t]), np.concatenate([S_t[0]*(E_l/E_t[0])**slope, S_t])

EC, SC = extend_low(nist_C_E, nist_C_S)
EH, SH = extend_low(nist_H_E, nist_H_S)
iC = interp1d(np.log(EC), np.log(SC), kind='cubic', fill_value='extrapolate')
iH = interp1d(np.log(EH), np.log(SH), kind='cubic', fill_value='extrapolate')

def S_CD2(E):
    E = np.clip(np.atleast_1d(float(E)), 0.001, 10.0)
    return 0.75 * np.exp(iC(np.log(E))) + 0.25 * np.exp(iH(np.log(E)))

rho = 1.06
E0  = 3.271

def dEdx(x, E):
    return [0.0] if E[0] <= 0.001 else [-S_CD2(E[0])[0] * rho * 0.1]

def stop(x, E): return E[0] - 0.01
stop.terminal = True; stop.direction = -1

sol = solve_ivp(dEdx, [0, 25], [E0], events=stop,
                max_step=0.02, dense_output=True, rtol=1e-9)
x_stop = float(sol.t_events[0][0])

# ── Full curve ────────────────────────────────────────────────────────────────
x_full = np.linspace(0, x_stop, 800)
E_full = np.maximum(sol.sol(x_full)[0], 0)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(x_full, E_full, color='royalblue', lw=2.5, label='Remaining energy')
ax.fill_between(x_full, E_full, alpha=0.12, color='royalblue')

# 10 um marker
E_10 = float(sol.sol(10)[0])
ax.axvline(10, color='gray', lw=1, ls='--', alpha=0.6)
ax.plot(10, E_10, 'o', color='orange', ms=8, zorder=5)
ax.annotate(f'10 µm → {E_10:.3f} MeV',
            xy=(10, E_10), xytext=(10.6, E_10 + 0.25), fontsize=10,
            arrowprops=dict(arrowstyle='->', color='orange'))

# range marker
ax.annotate(f'Range = {x_stop:.1f} µm',
            xy=(x_stop, 0), xytext=(x_stop - 3.8, 0.4), fontsize=10,
            color='firebrick',
            arrowprops=dict(arrowstyle='->', color='firebrick'))

ax.set_xlabel('CD₂ Thickness (µm)', fontsize=13)
ax.set_ylabel('Remaining Energy (MeV)', fontsize=13)
ax.set_title(r'$^4$He (α)  3.271 MeV  in CD₂  (ρ = 1.06 g/cm³)', fontsize=13)
ax.set_xlim(0, x_stop * 1.05)
ax.set_ylim(-0.1, 3.5)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

plt.tight_layout()
import os; os.makedirs("figures", exist_ok=True)
outpng = "figures/alpha_cd2_energy_vs_thickness.png"
plt.savefig(outpng, dpi=150)
print(f"Saved: {outpng}")
plt.show()
