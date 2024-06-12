#!/usr/bin/env python3
import mms
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'text.usetex': True})
PAD = 0.5

df1 = mms.run_temporal('mms_temporal.i', 4, console=False, executable='../../')
df2 = mms.run_temporal('mms_temporal.i', 4, 'Executioner/scheme=bdf2', console=False,
                       executable='../../')
fig = mms.ConvergencePlot(xlabel='$\Delta$t', ylabel='$L_2$ Error')
fig.plot(df1, label='1st Order (Implicit Euler)', marker='o', markersize=8)
fig.plot(df2, label='2nd Order (Backward Difference)', marker='o', markersize=8)
plt.tight_layout(pad=PAD)
fig.save('nuclide_mms_temp.png')
