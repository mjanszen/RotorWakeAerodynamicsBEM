import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Task ?? : Plot Circulation
# =============================================================================

# Load Results from file
df_results = pd.read_csv("../data/BEM_results.dat", delimiter=',')
R = 50
NBlades = 3

# Operational Conditions
TSR_ref = [6, 8, 10]
pitch_angle = -2
v0 = 10

# Plot Non-Dimensional Circulation
fig, ax = plt.subplots(figsize=(12, 6))
for tsr in TSR_ref:
    temp = df_results.loc[(df_results['tsr'] == tsr)
                          & (df_results['pitch'] == pitch_angle)
                          & (df_results['v0'] == v0)]
    ax.plot(temp['r_centre'] / R, temp['circulation'] / (np.pi * R * temp["v0"] / (tsr * NBlades)), label=f'tsr:{tsr}')

ax.set_title('Non-dimensionalized circulation')
ax.set_xlabel(r'$r/R$')
ax.set_ylabel(r'$\frac{\Gamma }{(\pi U_{\infty}^2 / \Omega B)}$')
ax.grid()
ax.legend()
plt.show()
