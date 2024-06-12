#!/usr/bin/env python3
import mms
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'text.usetex': True})
PAD = 0.5

# A near source region of 10.0%.
fig = mms.ConvergencePlot(xlabel='Element Size ($h$)', ylabel='$L_2$ Error')

df = mms.run_spatial('nsr_10_percent_1_quad.i', 4, 'xs=0.0', console=False, executable='../../../')
fig.plot(df, label='1 Quadrant', marker='o', markersize=8)
df = mms.run_spatial('nsr_10_percent_2_quad.i', 4, 'xs=0.0', console=False, executable='../../../')
fig.plot(df, label='2 Quadrants', marker='o', markersize=8)
df = mms.run_spatial('nsr_10_percent_3_quad.i', 4, 'xs=0.0', console=False, executable='../../../')
fig.plot(df, label='3 Quadrants', marker='o', markersize=8)
df = mms.run_spatial('nsr_10_percent_4_quad.i', 4, 'xs=0.0', console=False, executable='../../../')
fig.plot(df, label='4 Quadrants', marker='o', markersize=8)

plt.tight_layout(pad=PAD)
fig.save('10_percent_nsr_mult_q.png')
plt.show()
