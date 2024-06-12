import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 32})
plt.rcParams.update({'text.usetex': True})

WIDTH = 16
HEIGHT = 10
PAD = 0.5
LINEWIDTH = 3

anal_df = pd.read_csv('test_out_LineAnal_0001.csv', sep=',')
anal = anal_df[anal_df.columns[0]].to_numpy()

dims = anal_df[anal_df.columns[2]].to_numpy()

num_1_df = pd.read_csv('test_out_LineNum_1_0001.csv', sep=',')
num_1 = num_1_df[num_1_df.columns[1]].to_numpy()

num_2_df = pd.read_csv('test_out_LineNum_2_0001.csv', sep=',')
num_2 = num_2_df[num_2_df.columns[1]].to_numpy()

num_3_df = pd.read_csv('test_out_LineNum_3_0001.csv', sep=',')
num_3 = num_3_df[num_3_df.columns[1]].to_numpy()

err_1_df = pd.read_csv('test_out_LineError_1_0001.csv', sep=',')
err_1 = err_1_df[err_1_df.columns[0]].to_numpy()
err_2_df = pd.read_csv('test_out_LineError_2_0001.csv', sep=',')
err_2 = err_2_df[err_2_df.columns[0]].to_numpy()
err_3_df = pd.read_csv('test_out_LineError_3_0001.csv', sep=',')
err_3 = err_3_df[err_3_df.columns[0]].to_numpy()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("")
ax1.set_xlabel("$x$ ($cm$)")
ax1.set_ylabel("Scalar Flux ($cm^{-2}s^{-1}$)")
ax1.set_yscale('log')
ax1.plot(dims, num_1, label='Gnat, $N = 20$', linestyle="--", linewidth=LINEWIDTH)
ax1.plot(dims, num_2, label='Gnat, $N = 40$', linestyle="--", linewidth=LINEWIDTH)
ax1.plot(dims, num_3, label='Gnat, $N = 80$', linestyle="--", linewidth=LINEWIDTH)
ax1.plot(dims, anal,  label='Analytical',     linestyle="-", linewidth=LINEWIDTH)

fig.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95))
fig.tight_layout()
plt.savefig("./plots/1D_analytical_flux.png", format='png')
plt.show()

fig, ax1 = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax1.set_title("")
ax1.set_xlabel("$x$ ($cm$)")
ax1.set_ylabel("Scalar Flux Error ($\%$)")
ax1.set_yscale('linear')
ax1.plot(dims, err_1, label='Gnat, $N = 20$', linestyle="--", linewidth=LINEWIDTH)
ax1.plot(dims, err_2, label='Gnat, $N = 40$', linestyle="--", linewidth=LINEWIDTH)
ax1.plot(dims, err_3, label='Gnat, $N = 80$', linestyle="--", linewidth=LINEWIDTH)

fig.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95))
fig.tight_layout()
plt.savefig("./plots/1D_analytical_flux_error.png", format='png')
plt.show()
