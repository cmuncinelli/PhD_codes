import numpy as np
import time

"""
================================================================================
Ratio uncertainty propagation via Monte Carlo resampling (bootstrap method)
================================================================================

This function computes the ratio R = N / D between two binned quantities
(numerator and denominator) and estimates *asymmetric uncertainties* on the
ratio using a Monte Carlo resampling (bootstrap) technique.

Why this approach?
------------------
The ratio of two quantities is a *non-linear* function. Even if both numerator
and denominator have symmetric (Gaussian) uncertainties, the resulting ratio
generally has an *asymmetric* uncertainty distribution. A simple linear error
propagation or a symmetric standard deviation can therefore be misleading.

Instead, we:
  1. Treat statistical uncertainties in the denominator as Gaussian.
  2. Treat systematic uncertainties in the denominator as asymmetric shifts.
  3. Resample the denominator many times (bootstrap). (Just the denominator though!)
  4. Compute the ratio for each resampled spectrum.
  5. Extract a central confidence interval using percentiles.

Choice of percentiles (16% and 84%)
----------------------------------
For a perfectly Gaussian distribution, the interval between the 16th and 84th
percentiles corresponds to a central 68% confidence interval (~ +/- 1 sigma).
Even though the ratio distribution is not guaranteed to be Gaussian, these
percentiles provide a robust, *distribution-agnostic* estimator of the central
uncertainty band. This is generally preferable to quoting a symmetric standard
deviation for an asymmetric distribution.

Assumptions
-----------
- Numerator and denominator are evaluated at the *same bin centers*.
  No interpolation or extrapolation is needed.
- Statistical uncertainties are uncorrelated between bins.
- Statistical uncertainties are Gaussian (reasonable for large statistics).
- Systematic uncertainties in the denominator are treated as bin-wise shifts.
- Numerator uncertainties are propagated linearly after the ratio is formed.

Inputs
------
numerator : array-like
    Bin contents of the numerator.
denominator : array-like
    Bin contents of the denominator.
num_err : array-like
    Total uncertainty of the numerator (stat and sys will be propagated similarly.
    The asymmetric part is the one being divided in the denominator).
den_stat_err : array-like
    Statistical uncertainty of the denominator (symmetric).
den_sys_err : array-like or float, optional
    Systematic uncertainty of the denominator. Can be:
      - a single scalar (symmetric),
      - a vector (symmetric),
      - or a tuple (minus_sys, plus_sys) for asymmetric systematics.
n_samples : int
    Number of samples to be taken in the boostrapping procedure

Returns
-------
ratio_uncertainty : ndarray, shape (2, N)
    Asymmetric uncertainties on the ratio:
      ratio_uncertainty[0] --> lower error
      ratio_uncertainty[1] --> upper error
ratio : ndarray
    Central value of the ratio N / D.
================================================================================
"""


def ratio_uncertainty_mc(
    numerator,
    denominator,
    num_err,
    den_stat_err,
    den_sys_err=0.0,
    n_samples=int(5e4),
):
    # -------------------------------------------------------------------------
    # Convert inputs explicitly to NumPy arrays
    # -------------------------------------------------------------------------
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    num_err = np.asarray(num_err, dtype=float)
    den_stat_err = np.asarray(den_stat_err, dtype=float)

    # Handle systematic uncertainties in the denominator
    # Allow for symmetric or asymmetric systematics
    if isinstance(den_sys_err, (tuple, list)): # Just checking if den_sys_err is actually a vector of vectors, containing the asymmetric parts
        den_sys_minus = np.asarray(den_sys_err[0], dtype=float)
        den_sys_plus = np.asarray(den_sys_err[1], dtype=float)
    else:
        den_sys_minus = np.asarray(den_sys_err, dtype=float)
        den_sys_plus = np.asarray(den_sys_err, dtype=float)

    # -------------------------------------------------------------------------
    # Central value of the ratio
    # -------------------------------------------------------------------------
    ratio = numerator / denominator

    # Numerator contribution propagated linearly:
    #   R = N / D  -->  sigma_R_N = sigma_N / D
    num_uncertainty_term = num_err / denominator

    # -------------------------------------------------------------------------
    # Bootstrap resampling of the denominator
    # -------------------------------------------------------------------------
    print(f"\nStarting bootstrapping procedure with {n_samples} samples.")
    start_time = time.time()

    ratio_samples = []

    for _ in range(n_samples):
        # Statistical fluctuations:
        # One Gaussian fluctuation per bin
        stat_fluct = np.random.normal(loc=0.0, scale=den_stat_err)

        # Shifted denominator including:
        # - statistical fluctuation
        # - systematic upward / downward shifts
        den_up = denominator + stat_fluct + den_sys_plus
        den_down = denominator + stat_fluct - den_sys_minus

        # Compute ratio for both systematic extremes
        ratio_samples.append(numerator / den_up)
        ratio_samples.append(numerator / den_down)

    ratio_samples = np.asarray(ratio_samples)

    # -------------------------------------------------------------------------
    # Extract central confidence interval from percentiles
    # -------------------------------------------------------------------------
    q16, q84 = np.percentile(ratio_samples, [16, 84], axis=0)

    # Asymmetric denominator-driven uncertainty
    bootstrapping_upper_error = np.abs(q84 - ratio)
    bootstrapping_lower_error = np.abs(q16 - ratio)

    # Combine numerator and denominator contributions in quadrature
    upper_error = np.sqrt(
        bootstrapping_upper_error**2 + num_uncertainty_term**2
    )
    lower_error = np.sqrt(
        bootstrapping_lower_error**2 + num_uncertainty_term**2
    )

    # -------------------------------------------------------------------------
    # Special case: numerator == denominator --> ratio = 1 with zero uncertainty
    # -------------------------------------------------------------------------
    # (Just for paranoia haha)
    if np.allclose(numerator, denominator, atol=1e-6):
        upper_error[:] = 0.0
        lower_error[:] = 0.0

    ratio_uncertainty = np.vstack([lower_error, upper_error])

    elapsed = time.time() - start_time
    print(
        f"Ending bootstrapping procedure. "
        f"Took {elapsed:.2f} s total "
        f"({elapsed / n_samples:.2e} s per sample)."
    )

    return ratio_uncertainty, ratio