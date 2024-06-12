import numpy as np
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'font.weight': 'black'})
plt.rcParams.update({'text.usetex': True})

WIDTH = 16
HEIGHT = 10
PAD = 0.5
LINEWIDTH = 2.0

# Quad meshes.
h_quads = np.array([4.695357e-01, 2.394210e-01, 1.211106e-01])
sigma_0_quads = np.array([3.053544e-03, 5.889862e-04, 1.130405e-04])
sigma_1_quads = np.array([2.994832e-03, 6.179738e-04, 1.284101e-04])
sigma_2_quads = np.array([2.898093e-02, 9.006073e-03, 2.016663e-03])

sigma_0_quads_coef = np.polyfit(np.log10(h_quads), np.log10(sigma_0_quads), 1)
sigma_1_quads_coef = np.polyfit(np.log10(h_quads), np.log10(sigma_1_quads), 1)
sigma_2_quads_coef = np.polyfit(np.log10(h_quads), np.log10(sigma_2_quads), 1)

# Tri meshes.
h_tris = np.array([4.751901e-01, 2.395790e-01, 1.202166e-01])
sigma_0_tris = np.array([5.014959e-03, 1.458384e-03, 3.494491e-04])
sigma_1_tris = np.array([3.370195e-03, 9.216808e-04, 2.240193e-04])
sigma_2_tris = np.array([2.672866e-02, 7.341917e-03, 1.644150e-03])

sigma_0_tris_coef = np.polyfit(np.log10(h_tris), np.log10(sigma_0_tris), 1)
sigma_1_tris_coef = np.polyfit(np.log10(h_tris), np.log10(sigma_1_tris), 1)
sigma_2_tris_coef = np.polyfit(np.log10(h_tris), np.log10(sigma_2_tris), 1)

fig, ax = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax.plot(h_quads, sigma_0_quads, label='Quads, $\sigma_{t,a} = 0.0$: ' + '{:.{precision}f}'.format(sigma_0_quads_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'o')
ax.plot(h_quads, sigma_1_quads, label='Quads, $\sigma_{t,a} = 1.0$: ' + '{:.{precision}f}'.format(sigma_1_quads_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'v')
ax.plot(h_quads, sigma_2_quads, label='Quads, $\sigma_{t,a} = 10.0$: ' + '{:.{precision}f}'.format(sigma_2_quads_coef[0], precision=3) , linewidth=LINEWIDTH, marker = '^')

ax.plot(h_tris, sigma_0_tris, label='Tris, $\sigma_{t,a} = 0.0$: ' + '{:.{precision}f}'.format(sigma_0_tris_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'o')
ax.plot(h_tris, sigma_1_tris, label='Tris, $\sigma_{t,a} = 1.0$: ' + '{:.{precision}f}'.format(sigma_1_tris_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'v')
ax.plot(h_tris, sigma_2_tris, label='Tris, $\sigma_{t,a} = 10.0$: ' + '{:.{precision}f}'.format(sigma_2_tris_coef[0], precision=3) , linewidth=LINEWIDTH, marker = '^')

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('Average Element Size ($h$)')
ax.set_ylabel('$L_{2}$ Error')
ax.legend(loc='lower right', ncol = 2)
ax.grid(True, which='both', color=[0.8]*3)

fig.tight_layout(pad=PAD)
fig.set_size_inches(WIDTH, HEIGHT)
plt.savefig("all.png", format='png')
plt.show()
