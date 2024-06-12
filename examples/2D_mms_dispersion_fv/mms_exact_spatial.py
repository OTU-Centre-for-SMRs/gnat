#!/usr/bin/env python3
import mms
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'text.usetex': True})
PAD = 0.5

df1 = mms.run_spatial('mms_spatial.i', 4, 'FVKernels/Advection/advected_interp_method=upwind', console=False, executable='../../')
df2 = mms.run_spatial('mms_spatial.i', 4, 'FVKernels/Advection/advected_interp_method=vanLeer', console=False, executable='../../')

fig = mms.ConvergencePlot(xlabel='Element Size ($h$)', ylabel='$L_2$ Error')
fig.plot(df1, label='1st Order (Constant Upwinding)', marker='o', markersize=8)
fig.plot(df2, label='2nd Order (van Leer)', marker='o', markersize=8)
plt.tight_layout(pad=PAD)
fig.save('nuclide_mms_spatial.png')
