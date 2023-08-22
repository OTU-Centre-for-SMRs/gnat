import numpy as np
import scipy.special as sc
from matplotlib import pyplot as plt
import matplotlib as mat
import pandas as pd

plt.rcParams.update({'font.size': 36})
plt.rcParams.update({'text.usetex': True})

# The analytical solution for the scalar flux in a slab reactor.
# q is the intensity of the isotropic neutron source [cm^{-3} s^{-1}].
q = 1000.0
# a is the width of the system, which varies from -a / 2 to a / 2 [cm].
a = 100.0
# src_pos is the position of the source within the bounds of -a / 2 to a / 2.
src_pos = 0.0
# sigma_t is the macroscopic total cross-section [cm^{-1}].
sigma_t = 0.1

# The scalar flux in the system.
def analytical_scalar_flux(x):
  optical_thickness = np.abs(src_pos - x) * sigma_t
  return sc.exp1(optical_thickness) * q / 2.0

# Parse the CSV files.
df_25a_25s = pd.read_csv('data/25a_25s.csv')
df_25a_51s = pd.read_csv('data/25a_51s.csv')
df_25a_101s = pd.read_csv('data/25a_101s.csv')
df_25a_251s = pd.read_csv('data/25a_251s.csv')
df_25a_501s = pd.read_csv('data/25a_501s.csv')
df_25a_1001s = pd.read_csv('data/25a_1001s.csv')

x_25a_25s = df_25a_25s['Points:0']
flux_25a_25s = df_25a_25s['flux_moment_1_0_0']
x_25a_51s = df_25a_51s['Points:0']
flux_25a_51s = df_25a_51s['flux_moment_1_0_0']
x_25a_101s = df_25a_101s['Points:0']
flux_25a_101s = df_25a_101s['flux_moment_1_0_0']
x_25a_251s = df_25a_251s['Points:0']
flux_25a_251s = df_25a_251s['flux_moment_1_0_0']
x_25a_501s = df_25a_501s['Points:0']
flux_25a_501s = df_25a_501s['flux_moment_1_0_0']
x_25a_1001s = df_25a_1001s['Points:0']
flux_25a_1001s = df_25a_1001s['flux_moment_1_0_0']

df_50a_25s = pd.read_csv('data/50a_25s.csv')
df_50a_51s = pd.read_csv('data/50a_51s.csv')
df_50a_101s = pd.read_csv('data/50a_101s.csv')
df_50a_251s = pd.read_csv('data/50a_251s.csv')
df_50a_501s = pd.read_csv('data/50a_501s.csv')
df_50a_1001s = pd.read_csv('data/50a_1001s.csv')

x_50a_25s = df_50a_25s['Points:0']
flux_50a_25s = df_50a_25s['flux_moment_1_0_0']
x_50a_51s = df_50a_51s['Points:0']
flux_50a_51s = df_50a_51s['flux_moment_1_0_0']
x_50a_101s = df_50a_101s['Points:0']
flux_50a_101s = df_50a_101s['flux_moment_1_0_0']
x_50a_251s = df_50a_251s['Points:0']
flux_50a_251s = df_50a_251s['flux_moment_1_0_0']
x_50a_501s = df_50a_501s['Points:0']
flux_50a_501s = df_50a_501s['flux_moment_1_0_0']
x_50a_1001s = df_50a_1001s['Points:0']
flux_50a_1001s = df_50a_1001s['flux_moment_1_0_0']

df_75a_1001s = pd.read_csv('data/75a_1001s.csv')

x_75a_1001s = df_75a_1001s['Points:0']
flux_75a_1001s = df_75a_1001s['flux_moment_1_0_0']

df_100a_25s = pd.read_csv('data/100a_25s.csv')
df_100a_51s = pd.read_csv('data/100a_51s.csv')
df_100a_101s = pd.read_csv('data/100a_101s.csv')
df_100a_251s = pd.read_csv('data/100a_251s.csv')
df_100a_501s = pd.read_csv('data/100a_501s.csv')
df_100a_1001s = pd.read_csv('data/100a_1001s.csv')

x_100a_25s = df_100a_25s['Points:0']
flux_100a_25s = df_100a_25s['flux_moment_1_0_0']
x_100a_51s = df_100a_51s['Points:0']
flux_100a_51s = df_100a_51s['flux_moment_1_0_0']
x_100a_101s = df_100a_101s['Points:0']
flux_100a_101s = df_100a_101s['flux_moment_1_0_0']
x_100a_251s = df_100a_251s['Points:0']
flux_100a_251s = df_100a_251s['flux_moment_1_0_0']
x_100a_501s = df_100a_501s['Points:0']
flux_100a_501s = df_100a_501s['flux_moment_1_0_0']
x_100a_1001s = df_100a_1001s['Points:0']
flux_100a_1001s = df_100a_1001s['flux_moment_1_0_0']

ratio_25a_25s   = flux_25a_25s   / analytical_scalar_flux(x_25a_25s - (a / 2.0))
ratio_25a_51s   = flux_25a_51s   / analytical_scalar_flux(x_25a_51s - (a / 2.0))
ratio_25a_101s  = flux_25a_101s  / analytical_scalar_flux(x_25a_101s - (a / 2.0))
ratio_25a_251s  = flux_25a_251s  / analytical_scalar_flux(x_25a_251s - (a / 2.0))
ratio_25a_501s  = flux_25a_501s  / analytical_scalar_flux(x_25a_501s - (a / 2.0))
ratio_25a_1001s = flux_25a_1001s / analytical_scalar_flux(x_25a_1001s - (a / 2.0))

ratio_50a_25s   = flux_50a_25s   / analytical_scalar_flux(x_50a_25s - (a / 2.0))
ratio_50a_51s   = flux_50a_51s   / analytical_scalar_flux(x_50a_51s - (a / 2.0))
ratio_50a_101s  = flux_50a_101s  / analytical_scalar_flux(x_50a_101s - (a / 2.0))
ratio_50a_251s  = flux_50a_251s  / analytical_scalar_flux(x_50a_251s - (a / 2.0))
ratio_50a_501s  = flux_50a_501s  / analytical_scalar_flux(x_50a_501s - (a / 2.0))
ratio_50a_1001s = flux_50a_1001s / analytical_scalar_flux(x_50a_1001s - (a / 2.0))

ratio_75a_1001s = flux_75a_1001s / analytical_scalar_flux(x_75a_1001s - (a / 2.0))

ratio_100a_25s   = flux_100a_25s   / analytical_scalar_flux(x_100a_25s - (a / 2.0))
ratio_100a_51s   = flux_100a_51s   / analytical_scalar_flux(x_100a_51s - (a / 2.0))
ratio_100a_101s  = flux_100a_101s  / analytical_scalar_flux(x_100a_101s - (a / 2.0))
ratio_100a_251s  = flux_100a_251s  / analytical_scalar_flux(x_100a_251s - (a / 2.0))
ratio_100a_501s  = flux_100a_501s  / analytical_scalar_flux(x_100a_501s - (a / 2.0))
ratio_100a_1001s = flux_100a_1001s / analytical_scalar_flux(x_100a_1001s - (a / 2.0))

WIDTH = 16
HEIGHT = 10
PAD = 0.5
LINEWIDTH = 3

plt.figure(figsize=(WIDTH,HEIGHT))
plt.loglog(x_25a_1001s - (a / 2.0), ratio_25a_1001s, label='$N_{L} = 25$, $N_{e} = 1001$', linestyle="--", linewidth=LINEWIDTH)
plt.loglog(x_50a_1001s - (a / 2.0), ratio_50a_1001s, label='$N_{L} = 50$, $N_{e} = 1001$', linestyle="-.", linewidth=LINEWIDTH)
plt.loglog(x_75a_1001s - (a / 2.0), ratio_75a_1001s, label='$N_{L} = 75$, $N_{e} = 1001$', linestyle=":", linewidth=LINEWIDTH)
plt.loglog(x_100a_1001s - (a / 2.0), ratio_100a_1001s, label='$N_{L} = 100$, $N_{e} = 1001$', linestyle="-", linewidth=LINEWIDTH)
plt.legend()
plt.xlabel("Distance from the Center of the Source [cm]")
plt.ylabel("$\Phi / \Phi_{a}$")
plt.tight_layout(pad=PAD)
plt.savefig("flux_angular.svg", format='svg')
plt.savefig("flux_angular.png", format='png')
plt.show()

plt.figure(figsize=(WIDTH,HEIGHT))
plt.loglog(x_100a_25s - (a / 2.0), ratio_100a_25s, label='$N_{L} = 100$, $N_{e} = 25$', linestyle="--", linewidth=LINEWIDTH, markersize=5)
plt.loglog(x_100a_101s - (a / 2.0), ratio_100a_101s, label='$N_{L} = 100$, $N_{e} = 101$', linestyle="-.", linewidth=LINEWIDTH, markersize=5)
plt.loglog(x_100a_501s - (a / 2.0), ratio_100a_501s, label='$N_{L} = 100$, $N_{e} = 501$', linestyle=":", linewidth=LINEWIDTH, markersize=5)
plt.loglog(x_100a_1001s - (a / 2.0), ratio_100a_1001s, label='$N_{L} = 100$, $N_{e} = 1001$', linestyle="-", linewidth=LINEWIDTH, markersize=5)
plt.legend()
plt.xlabel("Distance from the Center of the Source [cm]")
plt.ylabel("$\Phi / \Phi_{a}$")
plt.tight_layout(pad=PAD)
plt.savefig("flux_spatial.svg", format='svg')
plt.savefig("flux_spatial.png", format='png')
plt.show()
