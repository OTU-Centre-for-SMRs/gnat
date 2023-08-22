import numpy as np
import scipy.special as sc
from matplotlib import pyplot as plt
import pandas as pd

plt.rcParams.update({'font.size': 36})
plt.rcParams.update({'text.usetex': True})

# a is the width of the system, which varies from 0 to a [cm].
a = 10.0
vel = 1.0
D = 1.0
decay = 1.0

r_1 = (-vel + np.sqrt(np.square(vel) + 4.0 * D * decay)) / (-2.0 * D)
r_2 = (-vel - np.sqrt(np.square(vel) + 4.0 * D * decay)) / (-2.0 * D)

def analytical_solution(x):
  B = (r_1 * (np.exp(r_1 * a))) / (r_1 * (np.exp(r_1 * a)) - r_2 * (np.exp(r_2 * a)))
  return (1.0 - B) * np.exp(r_1 * x) + B * np.exp(r_2 * x)

s25 = pd.read_csv('data/25s.csv')
s100 = pd.read_csv('data/100s.csv')
s500 = pd.read_csv('data/500s.csv')
s1000 = pd.read_csv('data/100s.csv')

x_25s = s25['Points:0']
u_25s = s25['u']
x_100s = s100['Points:0']
u_100s = s100['u']
x_500s = s500['Points:0']
u_500s = s500['u']
x_1000s = s1000['Points:0']
u_1000s = s1000['u']

ratio_25s   = u_25s   / analytical_solution(x_25s)
ratio_100s  = u_100s  / analytical_solution(x_100s)
ratio_500s  = u_500s  / analytical_solution(x_500s)
ratio_1000s = u_1000s / analytical_solution(x_1000s)

WIDTH = 16
HEIGHT = 10
PAD = 0.5
LINEWIDTH = 3

plt.figure(figsize=(WIDTH,HEIGHT))
plt.loglog(x_25s,   ratio_25s, label='$N_{e} = 25$', linestyle="-", linewidth=LINEWIDTH)
plt.loglog(x_100s,  ratio_100s, label='$N_{e} = 100$', linestyle="-.", linewidth=LINEWIDTH)
plt.loglog(x_500s,  ratio_500s, label='$N_{e} = 500$', linestyle=":", linewidth=LINEWIDTH)
plt.loglog(x_1000s, ratio_1000s, label='$N_{e} = 1000$', linestyle="--", linewidth=LINEWIDTH)
plt.legend()
plt.xlabel("Distance from the Left Boundary [cm]")
plt.ylabel("$N / N_{a}$")
plt.tight_layout(pad=PAD)
plt.savefig("nuclide_spatial.svg", format='svg')
plt.savefig("nuclide_spatial.png", format='png')
plt.show()
