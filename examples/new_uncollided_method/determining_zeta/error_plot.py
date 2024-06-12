import numpy as np
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'font.weight': 'black'})
plt.rcParams.update({'text.usetex': True})

WIDTH = 16
HEIGHT = 10
PAD = 0.5
LINEWIDTH = 2.0

xi = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])

# For void regions.
l2_error_norm_0_void = np.array([2.807142e-02, 2.472129e-02, 2.219298e-02, 2.026398e-02, 1.879529e-02, 1.769118e-02, 1.688019e-02, 1.630570e-02, 1.592131e-02, 1.568850e-02,
                                 1.557526e-02, 1.555525e-02, 1.560702e-02, 1.571333e-02, 1.586042e-02, 1.603746e-02, 1.623595e-02, 1.644927e-02, 1.667231e-02, 1.690112e-02])
l2_error_norm_1_void = np.array([1.086233e-02, 9.436198e-03, 8.447699e-03, 7.737110e-03, 7.219234e-03, 6.841113e-03, 6.566786e-03, 6.370539e-03, 6.233422e-03, 6.141222e-03,
                                 6.083175e-03, 6.051084e-03, 6.038687e-03, 6.041185e-03, 6.054900e-03, 6.077001e-03, 6.105312e-03, 6.138155e-03, 6.174233e-03, 6.212548e-03])
l2_error_norm_2_void = np.array([2.811875e-03, 2.349936e-03, 2.098142e-03, 1.960367e-03, 1.889282e-03, 1.858138e-03, 1.851023e-03, 1.858280e-03, 1.873981e-03, 1.894465e-03,
                                 1.917457e-03, 1.941541e-03, 1.965838e-03, 1.989808e-03, 2.013128e-03, 2.035611e-03, 2.057162e-03, 2.077739e-03, 2.097338e-03, 2.115979e-03])
l2_error_norm_3_void = np.array([6.148280e-04, 5.122187e-04, 4.831211e-04, 4.811285e-04, 4.891436e-04, 5.004684e-04, 5.124235e-04, 5.239478e-04, 5.346540e-04, 5.444406e-04,
                                 5.533268e-04, 5.613791e-04, 5.686788e-04, 5.753074e-04, 5.813409e-04, 5.868477e-04, 5.918879e-04, 5.965142e-04, 6.007724e-04, 6.047026e-04])
l2_error_norm_4_void = np.array([1.296648e-04, 1.214298e-04, 1.261382e-04, 1.319995e-04, 1.371560e-04, 1.414397e-04, 1.449768e-04, 1.479191e-04, 1.503936e-04, 1.524983e-04,
                                 1.543080e-04, 1.558793e-04, 1.572559e-04, 1.584714e-04, 1.595523e-04, 1.605198e-04, 1.613906e-04, 1.621786e-04, 1.628948e-04, 1.635486e-04])
l2_error_norm_5_void = np.array([3.043312e-05, 3.306416e-05, 3.542887e-05, 3.705548e-05, 3.820661e-05, 3.905722e-05, 3.970955e-05, 4.022507e-05, 4.064249e-05, 4.098724e-05,
                                 4.127673e-05, 4.152321e-05, 4.173558e-05, 4.192044e-05, 4.208281e-05, 4.222656e-05, 4.235470e-05, 4.246964e-05, 4.257333e-05, 4.266733e-05])

fig, ax = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax.plot(xi, l2_error_norm_0_void, label='Refinement 0: 396 elements',     linewidth=LINEWIDTH, marker = 'o')
ax.plot(xi, l2_error_norm_1_void, label='Refinement 1: 1,584 elements',   linewidth=LINEWIDTH, marker = 'v')
ax.plot(xi, l2_error_norm_2_void, label='Refinement 2: 6,336 elements',   linewidth=LINEWIDTH, marker = '^')
ax.plot(xi, l2_error_norm_3_void, label='Refinement 3: 25,344 elements',  linewidth=LINEWIDTH, marker = 's')
ax.plot(xi, l2_error_norm_4_void, label='Refinement 4: 101,376 elements', linewidth=LINEWIDTH, marker = 'D')
ax.plot(xi, l2_error_norm_5_void, label='Refinement 5: 405,504 elements', linewidth=LINEWIDTH, marker = 'p')
ax.set_yscale('log')
ax.set_xlabel('$\\zeta$')
ax.set_ylabel('$L_{2}$ Error')
ax.set_xticks(np.arange(min(xi), max(xi) + 0.1, 0.1))
ax.set_ylim(bottom = 2e-5, top = 1.3e-1)
ax.legend(loc='upper right', ncol = 2)
ax.grid(True, which='both', color=[0.8]*3)

fig.tight_layout(pad=PAD)
fig.set_size_inches(WIDTH, HEIGHT)
plt.savefig("zeta_void.png", format='png')
plt.show()

# For a cross-section of 10^{-3}
l2_error_norm_0_nonvoid_1 = np.array([3.166948e-02, 3.149458e-02, 3.108717e-02, 3.055893e-02, 2.997805e-02, 2.938356e-02, 2.879730e-02, 2.823106e-02, 2.769085e-02, 2.717928e-02,
                                      2.669698e-02, 2.624343e-02, 2.581749e-02, 2.541770e-02, 2.504244e-02, 2.469007e-02, 2.435896e-02, 2.404758e-02, 2.375448e-02, 2.347829e-02])
l2_error_norm_1_nonvoid_1 = np.array([2.106758e-02, 1.462423e-02, 1.160444e-02, 9.800734e-03, 8.579431e-03, 7.687080e-03, 7.001143e-03, 6.454510e-03, 6.007019e-03, 5.632998e-03,
                                      5.315173e-03, 5.041440e-03, 4.803028e-03, 4.593409e-03, 4.407608e-03, 4.241761e-03, 4.092812e-03, 3.958313e-03, 3.836272e-03, 3.725055e-03])
l2_error_norm_2_nonvoid_1 = np.array([4.430640e-03, 2.648531e-03, 1.973561e-03, 1.602587e-03, 1.362328e-03, 1.191745e-03, 1.063330e-03, 9.626644e-04, 8.813706e-04, 8.142098e-04,
                                      7.577163e-04, 7.094948e-04, 6.678311e-04, 6.314621e-04, 5.994351e-04, 5.710162e-04, 5.456303e-04, 5.228199e-04, 5.022158e-04, 4.835171e-04])
l2_error_norm_3_nonvoid_1 = np.array([7.014674e-04, 4.146160e-04, 3.046148e-04, 2.438404e-04, 2.045614e-04, 1.768322e-04, 1.561082e-04, 1.399867e-04, 1.270662e-04, 1.164692e-04,
                                      1.076161e-04, 1.001073e-04, 9.365790e-05, 8.805894e-05, 8.315346e-05, 7.882117e-05, 7.496828e-05, 7.152051e-05, 6.841823e-05, 6.561302e-05])
l2_error_norm_4_nonvoid_1 = np.array([1.063650e-04, 6.208477e-05, 4.481466e-05, 3.532370e-05, 2.925608e-05, 2.502256e-05, 2.189352e-05, 1.948385e-05, 1.757000e-05, 1.601289e-05,
                                      1.472130e-05, 1.363279e-05, 1.270317e-05, 1.190025e-05, 1.120001e-05, 1.058415e-05, 1.003850e-05, 9.551892e-06, 9.115405e-06, 8.721837e-06])
l2_error_norm_5_nonvoid_1 = np.array([1.589453e-05, 9.011778e-06, 6.367521e-06, 4.944512e-06, 4.050831e-06, 3.436181e-06, 2.987123e-06, 2.644556e-06, 2.374593e-06, 2.156380e-06,
                                      1.976370e-06, 1.825374e-06, 1.696939e-06, 1.586393e-06, 1.490276e-06, 1.405965e-06, 1.331438e-06, 1.265111e-06, 1.205723e-06, 1.152259e-06])

fig, ax = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax.plot(xi, l2_error_norm_0_nonvoid_1, label='Refinement 0: 396 elements',     linewidth=LINEWIDTH, marker = 'o')
ax.plot(xi, l2_error_norm_1_nonvoid_1, label='Refinement 1: 1,584 elements',   linewidth=LINEWIDTH, marker = 'v')
ax.plot(xi, l2_error_norm_2_nonvoid_1, label='Refinement 2: 6,336 elements',   linewidth=LINEWIDTH, marker = '^')
ax.plot(xi, l2_error_norm_3_nonvoid_1, label='Refinement 3: 25,344 elements',  linewidth=LINEWIDTH, marker = 's')
ax.plot(xi, l2_error_norm_4_nonvoid_1, label='Refinement 4: 101,376 elements', linewidth=LINEWIDTH, marker = 'D')
ax.plot(xi, l2_error_norm_5_nonvoid_1, label='Refinement 5: 405,504 elements', linewidth=LINEWIDTH, marker = 'p')
ax.set_yscale('log')
ax.set_xlabel('$\\zeta$')
ax.set_ylabel('$L_{2}$ Error')
ax.set_xticks(np.arange(min(xi), max(xi) + 0.1, 0.1))
ax.set_ylim(bottom = 1e-6, top = 4.8e-1)
ax.legend(loc='upper right', ncol = 2)
ax.grid(True, which='both', color=[0.8]*3)

fig.tight_layout(pad=PAD)
fig.set_size_inches(WIDTH, HEIGHT)
plt.savefig("zeta_nonvoid_1.png", format='png')
plt.show()

# For a cross-section of 10^{0}
l2_error_norm_0_nonvoid_2 = np.array([1.171509e-02, 1.171509e-02, 1.171509e-02, 1.171509e-02, 1.171509e-02, 1.171509e-02, 1.171509e-02, 1.171509e-02, 1.171509e-02, 1.171509e-02,
                                      1.129589e-02, 1.091297e-02, 1.056231e-02, 1.024041e-02, 9.944196e-03, 9.671005e-03, 9.418490e-03, 9.184596e-03, 8.967514e-03, 8.765651e-03])
l2_error_norm_1_nonvoid_2 = np.array([4.013621e-03, 4.013621e-03, 4.013621e-03, 4.013621e-03, 4.013620e-03, 3.655618e-03, 3.368311e-03, 3.131116e-03, 2.931217e-03, 2.760079e-03,
                                      2.611729e-03, 2.481821e-03, 2.367096e-03, 2.265049e-03, 2.173715e-03, 2.091527e-03, 2.017214e-03, 1.949738e-03, 1.888236e-03, 1.831986e-03])
l2_error_norm_2_nonvoid_2 = np.array([1.092135e-03, 1.092135e-03, 9.907391e-04, 8.445297e-04, 7.415365e-04, 6.638056e-04, 6.025570e-04, 5.528340e-04, 5.115661e-04, 4.767237e-04,
                                      4.468988e-04, 4.210784e-04, 3.985111e-04, 3.786266e-04, 3.609825e-04, 3.452303e-04, 3.310907e-04, 3.183377e-04, 3.067855e-04, 2.962803e-04])
l2_error_norm_3_nonvoid_2 = np.array([2.850308e-04, 2.191759e-04, 1.717905e-04, 1.427064e-04, 1.226352e-04, 1.078265e-04, 9.640691e-05, 8.731623e-05, 7.990258e-05, 7.374031e-05,
                                      6.853845e-05, 6.409072e-05, 6.024653e-05, 5.689318e-05, 5.394454e-05, 5.133362e-05, 4.900747e-05, 4.692368e-05, 4.504785e-05, 4.335180e-05])
l2_error_norm_4_nonvoid_2 = np.array([5.588236e-05, 3.623510e-05, 2.730181e-05, 2.202946e-05, 1.851848e-05, 1.600416e-05, 1.411256e-05, 1.263743e-05, 1.145513e-05, 1.048681e-05,
                                      9.679723e-06, 8.997206e-06, 8.412956e-06, 7.907603e-06, 7.466558e-06, 7.078624e-06, 6.735054e-06, 6.428922e-06, 6.154666e-06, 5.907770e-06])
l2_error_norm_5_nonvoid_2 = np.array([9.182571e-06, 5.564067e-06, 4.029746e-06, 3.169501e-06, 2.617179e-06, 2.232183e-06, 1.948474e-06, 1.730815e-06, 1.558651e-06, 1.419169e-06,
                                      1.303962e-06, 1.207280e-06, 1.125059e-06, 1.054341e-06, 9.929248e-07, 9.391349e-07, 8.916754e-07, 8.495275e-07, 8.118790e-07, 7.780748e-07])

fig, ax = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax.plot(xi, l2_error_norm_0_nonvoid_2, label='Refinement 0: 396 elements',     linewidth=LINEWIDTH, marker = 'o')
ax.plot(xi, l2_error_norm_1_nonvoid_2, label='Refinement 1: 1,584 elements',   linewidth=LINEWIDTH, marker = 'v')
ax.plot(xi, l2_error_norm_2_nonvoid_2, label='Refinement 2: 6,336 elements',   linewidth=LINEWIDTH, marker = '^')
ax.plot(xi, l2_error_norm_3_nonvoid_2, label='Refinement 3: 25,344 elements',  linewidth=LINEWIDTH, marker = 's')
ax.plot(xi, l2_error_norm_4_nonvoid_2, label='Refinement 4: 101,376 elements', linewidth=LINEWIDTH, marker = 'D')
ax.plot(xi, l2_error_norm_5_nonvoid_2, label='Refinement 5: 405,504 elements', linewidth=LINEWIDTH, marker = 'p')
ax.set_yscale('log')
ax.set_xlabel('$\\zeta$')
ax.set_ylabel('$L_{2}$ Error')
ax.set_xticks(np.arange(min(xi), max(xi) + 0.1, 0.1))
ax.set_ylim(bottom = 5e-7, top = 2.5e-1)
ax.legend(loc='upper right', ncol = 2)
ax.grid(True, which='both', color=[0.8]*3)

fig.tight_layout(pad=PAD)
fig.set_size_inches(WIDTH, HEIGHT)
plt.savefig("zeta_nonvoid_2.png", format='png')
plt.show()

# For fun, I guess. \zeta = 0.0.
h = np.array([1.414214e+00, 7.071069e-01, 3.535534e-01, 1.767767e-01, 8.838836e-02, 4.419418e-02])
sigma_0 = np.array([3.148120e-02, 4.820672e-02, 5.376586e-02, 4.991569e-02, 3.602374e-02, 1.697540e-02]) # sigma_{t} = 1e-3
sigma_1 = np.array([3.127685e-02, 4.408470e-02, 3.982541e-02, 2.243415e-02, 8.107307e-03, 2.339496e-03]) # sigma_{t} = 1e-2
sigma_2 = np.array([2.893509e-02, 2.323818e-02, 1.001239e-02, 3.040344e-03, 8.260108e-04, 2.162838e-04]) # sigma_{t} = 1e-1
sigma_3 = np.array([1.171509e-02, 4.013621e-03, 1.092135e-03, 2.850308e-04, 7.303575e-05, 1.869793e-05]) # sigma_{t} = 1e-0

sigma_0_coef = np.polyfit(np.log10(h), np.log10(sigma_0), 1)
sigma_1_coef = np.polyfit(np.log10(h), np.log10(sigma_1), 1)
sigma_2_coef = np.polyfit(np.log10(h), np.log10(sigma_2), 1)
sigma_3_coef = np.polyfit(np.log10(h), np.log10(sigma_3), 1)

fig, ax = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax.plot(h, sigma_0, label='$\sigma_{t} = 10^{-3}$: ' + '{:.{precision}f}'.format(sigma_0_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'o')
ax.plot(h, sigma_1, label='$\sigma_{t} = 10^{-2}$: ' + '{:.{precision}f}'.format(sigma_1_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'v')
ax.plot(h, sigma_2, label='$\sigma_{t} = 10^{-1}$: ' + '{:.{precision}f}'.format(sigma_2_coef[0], precision=3) , linewidth=LINEWIDTH, marker = '^')
ax.plot(h, sigma_3, label='$\sigma_{t} = 10^{0}$: ' + '{:.{precision}f}'.format(sigma_3_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 's')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('$h$')
ax.set_ylabel('$L_{2}$ Error')
ax.legend(loc='lower right', ncol = 2)
ax.grid(True, which='both', color=[0.8]*3)

fig.tight_layout(pad=PAD)
fig.set_size_inches(WIDTH, HEIGHT)
plt.savefig("zeta_eq_zero.png", format='png')
plt.show()

# For fun, I guess. \zeta = 1.0.
sigma_0 = np.array([2.717928e-02, 5.632998e-03, 8.142098e-04, 1.164692e-04, 1.601289e-05, 2.156380e-06]) # sigma_{t} = 1e-3
sigma_1 = np.array([2.663223e-02, 5.530749e-03, 8.019523e-04, 1.148967e-04, 1.580411e-05, 2.127090e-06]) # sigma_{t} = 1e-2
sigma_2 = np.array([2.275496e-02, 4.809748e-03, 7.187487e-04, 1.043953e-04, 1.442690e-05, 1.936637e-06]) # sigma_{t} = 1e-1
sigma_3 = np.array([1.171509e-02, 2.760079e-03, 4.767237e-04, 7.374031e-05, 1.048681e-05, 1.419169e-06]) # sigma_{t} = 1e-0

sigma_0_coef = np.polyfit(np.log10(h), np.log10(sigma_0), 1)
sigma_1_coef = np.polyfit(np.log10(h), np.log10(sigma_1), 1)
sigma_2_coef = np.polyfit(np.log10(h), np.log10(sigma_2), 1)
sigma_3_coef = np.polyfit(np.log10(h), np.log10(sigma_3), 1)

fig, ax = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax.plot(h, sigma_0, label='$\sigma_{t} = 10^{-3}$: ' + '{:.{precision}f}'.format(sigma_0_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'o')
ax.plot(h, sigma_1, label='$\sigma_{t} = 10^{-2}$: ' + '{:.{precision}f}'.format(sigma_1_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'v')
ax.plot(h, sigma_2, label='$\sigma_{t} = 10^{-1}$: ' + '{:.{precision}f}'.format(sigma_2_coef[0], precision=3) , linewidth=LINEWIDTH, marker = '^')
ax.plot(h, sigma_3, label='$\sigma_{t} = 10^{0}$: ' + '{:.{precision}f}'.format(sigma_3_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 's')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('$h$')
ax.set_ylabel('$L_{2}$ Error')
ax.legend(loc='lower right', ncol = 2)
ax.grid(True, which='both', color=[0.8]*3)

fig.tight_layout(pad=PAD)
fig.set_size_inches(WIDTH, HEIGHT)
plt.savefig("zeta_eq_one.png", format='png')
plt.show()

# For fun, I guess. \zeta = 2.0.
sigma_0 = np.array([2.347829e-02, 3.725055e-03, 4.835171e-04, 6.561302e-05, 8.721837e-06, 1.152259e-06]) # sigma_{t} = 1e-3
sigma_1 = np.array([2.292951e-02, 3.656042e-03, 4.766707e-04, 6.479120e-05, 8.614404e-06, 1.137097e-06]) # sigma_{t} = 1e-2
sigma_2 = np.array([1.903218e-02, 3.168910e-03, 4.303418e-04, 5.932128e-05, 7.909434e-06, 1.039195e-06]) # sigma_{t} = 1e-1
sigma_3 = np.array([8.765651e-03, 1.831986e-03, 2.962803e-04, 4.335180e-05, 5.907770e-06, 7.780748e-07]) # sigma_{t} = 1e-0

sigma_0_coef = np.polyfit(np.log10(h), np.log10(sigma_0), 1)
sigma_1_coef = np.polyfit(np.log10(h), np.log10(sigma_1), 1)
sigma_2_coef = np.polyfit(np.log10(h), np.log10(sigma_2), 1)
sigma_3_coef = np.polyfit(np.log10(h), np.log10(sigma_3), 1)

fig, ax = plt.subplots()
fig.set_size_inches(w=WIDTH,h=HEIGHT)
ax.plot(h, sigma_0, label='$\sigma_{t} = 10^{-3}$: ' + '{:.{precision}f}'.format(sigma_0_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'o')
ax.plot(h, sigma_1, label='$\sigma_{t} = 10^{-2}$: ' + '{:.{precision}f}'.format(sigma_1_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 'v')
ax.plot(h, sigma_2, label='$\sigma_{t} = 10^{-1}$: ' + '{:.{precision}f}'.format(sigma_2_coef[0], precision=3) , linewidth=LINEWIDTH, marker = '^')
ax.plot(h, sigma_3, label='$\sigma_{t} = 10^{0}$: ' + '{:.{precision}f}'.format(sigma_3_coef[0], precision=3) , linewidth=LINEWIDTH, marker = 's')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('$h$')
ax.set_ylabel('$L_{2}$ Error')
ax.legend(loc='lower right', ncol = 2)
ax.grid(True, which='both', color=[0.8]*3)

fig.tight_layout(pad=PAD)
fig.set_size_inches(WIDTH, HEIGHT)
plt.savefig("zeta_eq_two.png", format='png')
plt.show()
