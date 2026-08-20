"""
Alpha (4He) 3.271 MeV energy remaining after CD2 target.
Prints table from 0.6 to 10 um in 0.2 um steps.

Stopping power: NIST ASTAR, Bragg additivity (C 75% + D 25% by mass)
CD2 density: 1.06 g/cm3
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# NIST ASTAR: alpha stopping power in Carbon [MeV cm2/mg]
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

# NIST ASTAR: alpha stopping power in Hydrogen [MeV cm2/mg]
# D (deuterium) electronic stopping ~ H per unit mass
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
    S_l = S_t[0] * (E_l/E_t[0])**slope
    return np.concatenate([E_l, E_t]), np.concatenate([S_l, S_t])

EC, SC = extend_low(nist_C_E, nist_C_S)
EH, SH = extend_low(nist_H_E, nist_H_S)
iC = interp1d(np.log(EC), np.log(SC), kind='cubic', fill_value='extrapolate')
iH = interp1d(np.log(EH), np.log(SH), kind='cubic', fill_value='extrapolate')

def S_CD2(E):
    E = np.clip(np.atleast_1d(float(E)), 0.001, 10.0)
    return 0.75 * np.exp(iC(np.log(E))) + 0.25 * np.exp(iH(np.log(E)))

rho = 1.06   # g/cm3
E0  = 3.271  # MeV

def dEdx(x, E):
    return [0.0] if E[0] <= 0.001 else [-S_CD2(E[0])[0] * rho * 0.1]

def stop(x, E): return E[0] - 0.01
stop.terminal = True
stop.direction = -1

sol = solve_ivp(dEdx, [0, 25], [E0], events=stop,
                max_step=0.02, dense_output=True, rtol=1e-9)

x_stop = float(sol.t_events[0][0])
print(f"Range in CD2: {x_stop:.2f} um\n")

print(f"{'Thickness (um)':>16} {'mg/cm2':>10} {'E remain (MeV)':>16} {'E loss (MeV)':>14}")
print("-" * 60)

out_lines = ["um energy_remaining energy_loss mg_per_cm2"]
xs = np.arange(0.6, 10.01, 0.2)
for x in xs:
    e = max(float(sol.sol(x)[0]), 0)
    t = x * rho * 0.1
    loss = E0 - e
    print(f"{x:>16.1f} {t:>10.4f} {e:>16.4f} {loss:>14.4f}")
    out_lines.append(f"{x:.1f} {e:.6f} {loss:.6f} {t:.6f}")

outfile = "thickness_and_alpha_energy_remaing.dat"
with open(outfile, "w") as f:
    f.write("\n".join(out_lines) + "\n")
print(f"\nSaved: {outfile}")
