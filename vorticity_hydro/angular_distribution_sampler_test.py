# Python script to simulate sampling from the Λ -> pπ angular distribution
# and study how the sampled mean of cos(theta) approaches the true mean.
# Uses an inverse-CDF sampling formula to sample the angular distribution
# of the decay.
# Produces a plot of relative error vs sample size and a plot with the samples
# of the related distribution.

import numpy as np
import matplotlib.pyplot as plt
import argparse
# Parameters
alpha = 0.747           # alpha_Lambda
P_mag = 1e-2           # magnitude of polarization (test other variations of this, such as the 1e-2 that gets out of the plots by Vitor)
# Sofisticating with an argument parser:
parser = argparse.ArgumentParser()
parser.add_argument("value", type=float, help="A float value (e.g., 1e-1)")
args = parser.parse_args()
P_mag = args.value

# choose quantization axis aligned with P, so true mean <\cos\theta> = alpha/3 * P_mag
true_cosine_mean = alpha/3.0 * P_mag

# A parameter used in inverse CDF: A = alpha * |P|
A = alpha * P_mag

    # sample sizes to test:
max_N_samples = 1e7
sample_sizes_array = np.arange(1,int(max_N_samples))

print(sample_sizes_array)

# number of independent repeats for each sample size to estimate average error
repeats = 25

rel_errors = np.zeros((repeats, sample_sizes_array.size))
abs_errors = np.zeros((repeats, sample_sizes_array.size))

rng = np.random.default_rng(12345)

mean_cos_array = [[0] for n_repeats in range(repeats)] # Creates an array with 0's in the first entry
# print("\nStarter mean_cos_array:", mean_cos_array)

for j, N in enumerate(sample_sizes_array): # j is an index to acess the vectors, and N is the number of samples for this part of the loop
    if N % int(0.05 * max_N_samples) == 0:
        print("Now on cumulative sample number {} of {} ({}%)".format(N, max_N_samples, N * 1./max_N_samples * 100))
    for r in range(repeats):
        # u = rng.random(N)
        u = rng.random() # Non-vectorized form
        # inverse CDF sampling for cos(theta)
        if A < 1e-9:
            cos_t = 2.0 * u - 1.0
        else:
            # vectorized implementation of the provided formula
            sqrt_term = np.sqrt((1.0 - A)**2 + 4.0 * A * u)
            cos_t = (-1.0 + sqrt_term) / A
        
        current_mean = mean_cos_array[r][-1]*(N-1) + cos_t # Updates the last step's average with the current one's average
        # print("\t\tmean_cos_array[r][-1]", mean_cos_array[r][-1])
        # print("\t\tcos_t", cos_t)
        current_mean /= N # Averages by total number of samples that was taken up until this point

        mean_cos_array[r].append(current_mean)
        # print("\ncurrent_mean", current_mean)
        # print("\ntrue_cosine_mean", true_cosine_mean)
        abs_err = np.abs(current_mean - true_cosine_mean)
        abs_errors[r, j] = abs_err
        # relative error: handle mu ~ 0
        if np.abs(true_cosine_mean) < 1e-30:
            rel_errors[r, j] = np.nan
        else:
            rel_errors[r, j] = abs_err / np.abs(true_cosine_mean)

# summarize: median value and relative error across repeats:
median_mean_cos = np.median(np.array(mean_cos_array), axis=0)
median_rel_error = np.nanmedian(rel_errors, axis=0)
median_abs_error = np.median(abs_errors, axis=0)

# Plotting
fig, ax = plt.subplots(ncols=2, figsize=(8,4))
ax[0].loglog(sample_sizes_array, median_mean_cos[1:], linestyle='-')
ax[0].axhline(true_cosine_mean, linestyle='--', label='True value')
ax[0].set_xlabel('Number of samples N')
ax[0].set_ylabel('Median value of sample mean')
ax[0].grid(True, which='both', ls=':')
ax[0].set_ylim(true_cosine_mean/2, 1)
ax[0].legend()

ax[1].loglog(sample_sizes_array, median_rel_error, linestyle='-') # Skips the first value, which is a zero and is not needed.
ax[1].axhline(0.01, linestyle='--', label='1% relative error')
ax[1].set_xlabel('Number of samples N')
ax[1].set_ylabel('Median relative error of sample mean')
ax[1].grid(True, which='both', ls=':')
ax[1].legend()

fig.suptitle('Convergence of sampled <cos θ> to true mean (median over repeats) - P_mag = {:.1e}'.format(P_mag))
fig.tight_layout()
fig.savefig("/home/cicero/results/hydro_vorticity/smaller_tests/angular_distribution_sampler_test_PolIs{:.0e}_{:.0e}samples.png".format(P_mag, max_N_samples), dpi = 300)
# plt.show()